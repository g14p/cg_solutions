#include <cglib/rt/renderer.h>
#include <cglib/rt/intersection_tests.h>
#include <cglib/rt/raytracing_context.h>
#include <cglib/rt/intersection.h>
#include <cglib/rt/ray.h>
#include <cglib/rt/scene.h>
#include <cglib/rt/light.h>
#include <cglib/rt/material.h>
#include <cglib/rt/render_data.h>

/*
 * TODO: implement a ray-sphere intersection test here.
 * The sphere is defined by its center and radius.
 *
 * Return true, if (and only if) the ray intersects the sphere.
 * In this case, also fill the parameter t with the distance such that
 *    ray_origin + t * ray_direction
 * is the intersection point.
 */
bool intersect_sphere(
    glm::vec3 const& ray_origin,    // starting point of the ray
    glm::vec3 const& ray_direction, // direction of the ray
    glm::vec3 const& center,        // position of the sphere
    float radius,                   // radius of the sphere
    float* t)                       // output parameter which contains distance to the hit point
{
    cg_assert(t);
	cg_assert(std::fabs(glm::length(ray_direction) - 1.f) < EPSILON);

    //-------------------Begin Georg Solution -------------------------------
    //intersection leads to need of solving at**2+bt+c=0 for t
    // with coefficients a, b, c being:
    float a = glm::dot(ray_direction, ray_direction);
    float b = 2 * glm::dot(ray_direction, ray_origin - center);
    float c = glm::dot(ray_origin - center, ray_origin - center) - pow(radius, 2); 
    // analytically we get t= ( b +- sqrt( b**2 - 4ac )) / 2a
    float discriminant = b * b - 4*a*c;

    if (discriminant < 0) return false; // squareroot is complex -> no intersection
    else if (discriminant == 0) { // intersect in tangent style :)
        *t = b / (2*a);
        return true;
    }
    else{ //pierce the ball
        if (*t<0) return false; // we only consider halbgerade
        // case 1: squareroot counts negative --> ray shoots out of sphere --> t_small
        // case 2: squareroot counts positive --> ray shoots into sphere --> t_big
        // we prefer t_small for some reason. maybe the ray shoots to the user?!
        *t = (float)(-b  - sqrt(discriminant)) / (2 * a) ; 
        return true;
    }
    //-------------------End Georg Solution ---------------------------------
}

/*
 * emission characteristic of a spotlight
 */
glm::vec3 SpotLight::getEmission(
		glm::vec3 const& omega // world space direction
		) const
{
	cg_assert(std::fabs(glm::length(omega) - 1.f) < EPSILON);
 
	// TODO: implement a spotlight emitter as specified on the exercise sheet
	return glm::vec3(0.f);
}

glm::vec3 evaluate_phong(
	RenderData &data,			// class containing raytracing information
	MaterialSample const& mat,	// the material at position
	glm::vec3 const& P,			// world space position
	glm::vec3 const& N,			// normal at the position (already normalized)
	glm::vec3 const& V)			// view vector (already normalized)
{
	cg_assert(std::fabs(glm::length(N) - 1.f) < EPSILON);
	cg_assert(std::fabs(glm::length(V) - 1.f) < EPSILON);

	glm::vec3 contribution(0.f);
        glm::vec3 fake_contribution(0.3,0.1,0.1); // georgs humble faker

	// iterate over lights and sum up their contribution
	for (auto& light_uptr : data.context.get_active_scene()->lights) 
	{
		// TODO: calculate the (normalized) direction to the light
		const Light *light = light_uptr.get();
		glm::vec3 L(0.0f, 1.0f, 0.0f);
                // ----------------- Georg Begin Solution ------------------
                L = (light->getPosition() - P) / glm::length(light->getPosition() - P);
                float valid_light_angle = 1.0;
                // ----------------- Georg end  Solution ------------------


		float visibility = 1.f;
		if (data.context.params.shadows) {
			// TODO: check if light source is visible
                        // ----------------- Georg Begin Solution ------------------
                        visibility = visible(data, P, light->getPosition());
                        // ----------------- Georg end  Solution ------------------
		}

		glm::vec3 diffuse(0.f);
		if (data.context.params.diffuse) {
			// TODO: compute diffuse component of phong model
                        // ----------------- Georg Begin Solution ------------------
                        float cos_theta = glm::length(L) / glm::length(N); //select N as hypothenuse
                        cos_theta = glm::dot(L,N);
                        valid_light_angle = (cos_theta > 0) ? 1.f : 0.f;

                        diffuse = mat.k_d * std::max(0.f, cos_theta);

                        // ----------------- Georg end  Solution ------------------
		}

		glm::vec3 specular(0.f);
		if (data.context.params.specular) {
			// TODO: compute specular component of phong model
                        // ----------------- Georg Begin Solution ------------------
                        glm::vec3 R = -L + 2 * glm::dot(L,N)*N; // Reflectance Vector R
                        R = R / glm::length(R); //normalize R
                        float cos_psi = glm::length(R) / glm::length(V); // select V as hypothenuse
                        cos_psi = glm::dot(R, V);
                        specular = mat.k_s * pow(std::max(0.f, cos_psi),mat.n);
                        // ----------------- Georg end  Solution ------------------
		}

		glm::vec3 ambient = data.context.params.ambient ? mat.k_a : glm::vec3(0.0f);

		// TODO: modify this and implement the phong model as specified on the exercise sheet
                // ----------------- Georg Begin Solution ------------------
                float squared_dist = pow(glm::length(P - light->getPosition()), 2); //later weaken the lighting acc. to squared distance to light
                // squared_dist mit dot_product
                ambient = ambient / squared_dist;
                contribution += light->getEmission(-L)*visibility*valid_light_angle / squared_dist * (diffuse + specular);
                // ----------------- Georg end  Solution ------------------


		contribution += ambient * light->getPower();
	}
        
	return contribution;
}

glm::vec3 evaluate_reflection(
	RenderData &data,			// class containing raytracing information
	int depth,					// the current recursion depth
	glm::vec3 const& P,			// world space position
	glm::vec3 const& N,			// normal at the position (already normalized)
	glm::vec3 const& V)			// view vector (already normalized)
{
	// TODO: calculate reflective contribution by constructing and shooting a reflection ray.
	return glm::vec3(0.f);
}

glm::vec3 evaluate_transmission(
	RenderData &data,			// class containing raytracing information
	int depth,					// the current recursion depth
	glm::vec3 const& P,			// world space position
	glm::vec3 const& N,			// normal at the position (already normalized)
	glm::vec3 const& V,			// view vector (already normalized)
	float eta)					// the relative refraction index
{
	// TODO: calculate transmissive contribution by constructing and shooting a transmission ray.
	glm::vec3 contribution(0.f);
	return contribution;
}

glm::vec3 handle_transmissive_material_single_ior(
	RenderData &data,			// class containing raytracing information
	int depth,					// the current recursion depth
	glm::vec3 const& P,			// world space position
	glm::vec3 const& N,			// normal at the position (already normalized)
	glm::vec3 const& V,			// view vector (already normalized)
	float eta)					// the relative refraction index
{
	if (data.context.params.fresnel) {
		// TODO: replace with proper fresnel handling.
		return evaluate_transmission(data, depth, P, N, V, eta);
	}
	else {
		// just regular transmission
		return evaluate_transmission(data, depth, P, N, V, eta);
	}
}

glm::vec3 handle_transmissive_material(
	RenderData &data,					// class containing raytracing information
	int depth,							// the current recursion depth
	glm::vec3 const& P,					// world space position
	glm::vec3 const& N,					// normal at the position (already normalized)
	glm::vec3 const& V,					// view vector (already normalized)
	glm::vec3 const& eta_of_channel)	// relative refraction index of red, green and blue color channel
{
	if (data.context.params.dispersion && !(eta_of_channel[0] == eta_of_channel[1] && eta_of_channel[0] == eta_of_channel[2])) {
		// TODO: split ray into 3 rays (one for each color channel) and implement dispersion here
		glm::vec3 contribution(0.f);
		return contribution;
	}
	else {
		// dont handle transmission, take average refraction index instead.
		const float eta = 1.f/3.f*(eta_of_channel[0]+eta_of_channel[1]+eta_of_channel[2]);
		return handle_transmissive_material_single_ior(data, depth, P, N, V, eta);
	}
	return glm::vec3(0.f);
}
// CG_REVISION 5ec3fe57a58d0be911914e19912c1e83a3f0d67c
