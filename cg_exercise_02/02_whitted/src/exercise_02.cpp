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
#include <glm/detail/qualifier.hpp>
#include <glm/detail/type_vec.hpp>
#include <glm/geometric.hpp>

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
    double a = glm::dot(ray_direction, ray_direction);
    double b = 2 * glm::dot(ray_direction, ray_origin - center);
    double c = glm::dot(ray_origin - center, ray_origin - center) - std::pow(radius, 2); 
    //potentiall make synergies... TODO
    // analytically we get t= ( -b +- sqrt( b**2 - 4ac )) / 2a
    double discriminant = b * b - 4*a*c;
    if (discriminant < 0.0) return false; // squareroot definitly complex -> no intersection
    //else we pierce the ball (neglect assume intersect in tangent style and assume two intersections :)
    *t = (float)(-b - sqrt(discriminant)) / (2 * a) ; 
    if (*t<0.f) return false;
    // we only consider halbgerade
    // case 1: squareroot counts negative --> ray shoots out of sphere --> t_small
    // case 2: squareroot counts positive --> ray shoots into sphere --> t_big
    // we prefer t_small for some reason. maybe the ray shoots to the user?!
    return true;
//-------------------End Georg Solution ---------------------------------
/*
    cg_assert(t);

    const glm::vec3 e_c = ray_origin - center;
    const float c = glm::dot(e_c, e_c) - radius * radius;
    const float b = glm::dot(ray_direction, e_c);
    const float a = glm::dot(ray_direction, ray_direction);

    const float d = b * b - a * c;
    if (d >= 0.0f)
    {
        const float e = sqrt(d);
        const float f = 1.0f / a;
        const float t1 = (-b + e) * f;
        const float t2 = (-b - e) * f;

        const bool t1valid = t1 >= 0.0f;
        const bool t2valid = t2 >= 0.0f;
        *t = (t1valid && t2valid ? glm::min(t1, t2)
            : (t1valid ? t1
            : (t2valid ? t2 
            : -1)));

        if (*t >= 0)
        {
            return true;
        }
    }
    return false;
*/
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
        //Begin Georg Solution ----------------------------------------
	return getPower() * (falloff + 2) * std::pow(std::max(0.f, glm::dot(omega, direction)),falloff);
        //End Georg Solution ----------------------------------------
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

	// iterate over lights and sum up their contribution
	for (auto& light_uptr : data.context.get_active_scene()->lights) 
	{
                // calculate the (normalized) direction to the light
		const Light *light = light_uptr.get();
		glm::vec3 L(0.0f, 1.0f, 0.0f);
                L = glm::normalize(light->getPosition() - P); 
                // ignore the influence of light that doesnt come from a 'upper hemisphere' angle, called O in exercise :)
                float valid_light_angle = 1.0; // updated later 
                // ignore influence lights that are not visible on object, called S in exercise :-)
                float visibility = 1.f; // updated later 
                if (data.context.params.shadows) {
			// check if light source is visible 
                        visibility = visible(data, P, light->getPosition());
		}
		glm::vec3 diffuse(0.f);
		if (data.context.params.diffuse) {
			// compute diffuse component of phong model
                        // ----------------- Georg Begin Solution ------------------
                        float cos_theta = glm::dot(N,L);
                        valid_light_angle = (cos_theta > 0.f) ? 1.f : 0.f;
                        diffuse = mat.k_d * std::max(0.f, cos_theta);
                        // ----------------- Georg end  Solution ------------------
		}

		glm::vec3 specular(0.f);
		if (data.context.params.specular) {
			// compute specular component of phong model
                        // ----------------- Georg Begin Solution ------------------
                        glm::vec3 R = -L + 2 * glm::dot(L,N)*N; // Reflectance Vector R
                        R = glm::normalize(R);
                        float cos_psi = glm::dot(R, V);
                        specular = mat.k_s * pow(std::max(0.f, cos_psi),mat.n);
                        // ----------------- Georg end  Solution ------------------
		}

		glm::vec3 ambient = data.context.params.ambient ? mat.k_a : glm::vec3(0.0f);

                // We have now three parts of contribution> 
                //    * ambient
                //    * specular
                //    * diffuse
                // Each one is to understand as an intensity of light - as light has three colors (here only three) - its a 3D vector!



		// implement the phong model as specified on the exercise sheet
                // ----------------- Georg Begin Solution ------------------
                float dist = glm::length(P - light->getPosition()); 
                float squared_dist = dist * dist; 
		contribution += light->getPower() / squared_dist * ambient; 
                contribution += light->getEmission(-L)*visibility*valid_light_angle / squared_dist * diffuse;
                contribution += light->getEmission(-L)*visibility*valid_light_angle / squared_dist * specular;
                // ----------------- Georg end  Solution ------------------

	}
       	return contribution;
        /*
	cg_assert(std::fabs(glm::length(N) - 1.f) < EPSILON);
	cg_assert(std::fabs(glm::length(V) - 1.f) < EPSILON);

	glm::vec3 contribution(0.f);
	// iterate over lights and sum up their contribution
	for (auto& light : data.context.get_active_scene()->lights) {
		// TODO: calculate the (normalized) direction to the light
		const glm::vec3 L = glm::normalize(light->getPosition() - P);

		float visibility = 1.f;
		if (data.context.params.shadows) {
			// TODO: check if light source is visible
			if (!visible(data, P, light->getPosition())) {
				visibility = 0.f;
			}
		}

		glm::vec3 diffuse(0.f);
		if (data.context.params.diffuse) {
			// TODO: compute diffuse component of phong model
			if (visibility > 0.f) {
				diffuse = std::max(0.f, glm::dot(N, L)) * mat.k_d;
			}
		}

		glm::vec3 specular(0.f);
		if (data.context.params.specular) {
			// TODO: compute specular component of phong model
			if ((visibility > 0.f) && (glm::dot(L, N) > 0.f)) {
				const glm::vec3 R = reflect(L, N);
				specular = std::pow(std::max(0.f, glm::dot(R, V)), mat.n) * mat.k_s;
			}
		}

		glm::vec3 ambient = data.context.params.ambient ? mat.k_a : glm::vec3(0.0f);

		// TODO: modify this and implement the phong model as specified on the exercise sheet
		const float dist = glm::length(light->getPosition() - P);
		contribution += (visibility * (diffuse + specular) + ambient) * light->getEmission(-L) / (dist*dist);
	}

	return contribution;
        */

}

glm::vec3 evaluate_reflection(
	RenderData &data,			// class containing raytracing information
	int depth,					// the current recursion depth
	glm::vec3 const& P,			// world space position
	glm::vec3 const& N,			// normal at the position (already normalized)
	glm::vec3 const& V)			// view vector (already normalized)
{
	// TODO: calculate reflective contribution by constructing and shooting a reflection ray.
        float cos_theta = glm::dot(V, N);
        glm::vec3 R = 2.f * N - V;
        Ray ray(P+EPSILON*R,R); 
        glm::vec3 contribution = trace_recursive(data, ray, depth);

	return contribution;
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
