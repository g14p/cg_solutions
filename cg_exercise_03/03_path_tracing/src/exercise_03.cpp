#include <bit>
#include <cglib/rt/renderer.h>
#include <cglib/rt/intersection_tests.h>
#include <cglib/rt/raytracing_context.h>
#include <cglib/rt/intersection.h>
#include <cglib/rt/ray.h>
#include <cglib/rt/scene.h>
#include <cglib/rt/light.h>
#include <cglib/rt/material.h>
#include <cglib/rt/render_data.h>
#include <cglib/core/thread_local_data.h>
#include <cglib/core/assert.h>
#include <glm/geometric.hpp>
#include <limits>
/*
 * Creates a random sample on a unit sphere
 *
 * Parameters:
 * - u1: random number [0,1)
 * - u2: random number [0,1)
 */
glm::vec3
uniform_sample_sphere(float u1, float u2)
{
	// TODO AmbientOcclusion/IndirectIllumination: 
	// implement uniform sampling on a unit sphere
        float h = 1.f - 2*u1;
        float r = sqrt(1.f-h*h); 
	glm::vec3 d = glm::vec3(
                r*cos(2.f * std::numbers::pi * u2 ),
                h,
                r*sin(2.f * std::numbers::pi * u2 )
        );
        return glm::normalize(d);
}

/*
 * Creates a random sample on a hemisphere with
 * normal direction N
 *
 * Parameters:
 * - data: RenderData for access to random number generator
 * - N: main direction of the hemisphere
 */
glm::vec3
uniform_sample_hemisphere(RenderData& data, glm::vec3 const& N)
{
	// TODO AmbientOcclusion/IndirectIllumination: 
	// implement uniform sampling on a unit hemisphere.
	// data.tld->rand() creates uniform [0, 1] random numbers
	// TIP: use uniform_sample_sphere 
        
        cg_assert(1.f - glm::length(N) < EPSILON);
        
        glm::vec3 d = -N;
        bool d_valid = false;
        do {
            float u2 = data.tld->rand();
            float u1 = data.tld->rand();
            d = uniform_sample_sphere(u1, u2);
            d_valid = glm::dot(N, d) > 0.f;
        } // reselect d until if it is located on the wrong hemisphere (opposed to N) 
        while (!d_valid);
        return d;
}

float evaluate_ambient_occlusion(
	RenderData& data,           // class containing raytracing information
	glm::vec3 const& P,         // world space position
	glm::vec3 const& N)         // normal at the position (already normalized)
{
	// TODO AmbientOcclusion: compute ambient occlusion
	float ambient_occlusion = 0.f;
	for (int i = 0; i < data.context.params.ao_rays; ++i)
	{
            glm::vec3 direction = uniform_sample_hemisphere(data, N);
            float dist = max_unobstructed_distance(data, P, direction, N); 
            float c = 1.f / (data.context.params.half_ao_radius * data.context.params.half_ao_radius);
            float visibility = (dist == std::numeric_limits<float>::max()) ? 1.f
               : 1 - 1.f / (1 + c * dist * dist);
            ambient_occlusion += visibility * glm::dot(N, direction);

	}
        ambient_occlusion *= 2.f / data.context.params.ao_rays;

	return ambient_occlusion;
}

glm::vec3 evaluate_illumination_from_light(
	RenderData& data,           // class containing raytracing information
	MaterialSample const& mat,  // the material at position
	Light const& light,         // the light source
	glm::vec3 const& LP,        // a point on the light source
	glm::vec3 const& P,         // world space position
	glm::vec3 const& N,         // normal at the position (already normalized)
	glm::vec3 const& V)         // view vector (already normalized)
{
	glm::vec3 L = LP - P;                       // direction to the light
	const float dist2 = glm::dot(L, L);         // compute squared distance to light point
	L /= sqrt(dist2);                           // normalize direction

	float visibility = 1.f;
	if (data.context.params.shadows)
	{
		if (!visible(data, P, LP))
		{
			visibility = 0.f;
		}
	}

	auto incomingLight = visibility * light.getEmission(-L) / dist2;
	return evaluate_phong_BRDF(data, mat, L, N, V) * incomingLight;
}

glm::vec3 evaluate_illumination(
	RenderData& data,           // class containing raytracing information
	MaterialSample const& mat,  // the material at position
	glm::vec3 const& P,         // world space position
	glm::vec3 const& N,         // normal at the position (already normalized)
	glm::vec3 const& V,         // view vector (already normalized)
	int depth)                  // the current recursion depth
{
	glm::vec3 direct_illumination(0.f);
	if (!data.context.params.disable_direct || depth > 1)
	{
		for (auto& light : data.context.get_active_scene()->lights)
		{
			const glm::vec3 LP = light->getPosition();
			direct_illumination += evaluate_illumination_from_light(
				data, mat, *light, LP, P, N, V);
		}
		direct_illumination /= data.context.get_active_scene()->lights.size();
	}

	glm::vec3 indirect_illumination(0.f);
	if (data.context.params.indirect)
	{
            // TODO IndirectIllumination: compute indirect illumination with russian roulette
            
            for( int i = 0; i < data.context.params.indirect_rays; i++) {
                

                for (auto& light : data.context.get_active_scene()->lights)
		{
                        
                        glm::vec3 direction = uniform_sample_hemisphere(data, N);
                        Ray ray(P, direction);
                        //float cos_theta = glm::dot(N, direction); 
			const glm::vec3 LP = light->getPosition();
                        const glm::vec3 L = glm::normalize(LP - P);
                        glm::vec3 brdf = evaluate_phong_BRDF(data, mat, L, N, V);
                        if (data.context.params.russian_roulette)
                        {
                            if (depth >= data.context.params.russian_roulette_depth_thresh)
                            {
                                bool ray_survives = data.tld->rand() <= data.context.params.russian_roulette_surv_prop; 
                                if(ray_survives)
                                    indirect_illumination += brdf * trace_recursive(data, ray, depth + 1);
                            }
                        }
                        else {
                                indirect_illumination += brdf * trace_recursive(data, ray, depth + 1);
                        }
                        //std::cout<<"depth is " << depth << std::endl;

                }
		indirect_illumination /= data.context.get_active_scene()->lights.size();

           }
            indirect_illumination *= 2.f / data.context.params.indirect_rays;
	}

	return direct_illumination + indirect_illumination;
}// CG_REVISION d49323c0b80887ae263ae8263b5cfb0cc5a956b8
