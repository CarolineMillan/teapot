/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the entry point function for the program you need to create for lab two.
 * You should not need to modify this code.
 * It creates a framebuffer, loads an triangle mesh object, calls the drawing function to render the object and then outputs the framebuffer as a ppm file.
 *
 * On linux.bath.ac.uk:
 *
 * Compile the code using g++ -o lab2executable main_lab2.cpp framebuffer.cpp linedrawer.cpp polymesh.cpp -lm
 *
 * Execute the code using ./lab2executable
 *
 * This will produce an image file called test.ppm. You can convert this a png file for viewing using
 *
 * pbmropng test.ppm > test.png
 *
 * You are expected to fill in the missing code in polymesh.cpp.
 */



//    g++ -o csgexecutable main_csg.cpp framebuffer.cpp polymesh.cpp directional_light.cpp sphere.cpp phong.cpp -lm


#include "framebuffer.h"
#include "ray.h"
#include "hit.h"
#include "polymesh.h"
#include "sphere.h"
#include "light.h"
#include "directional_light.h"
#include "material.h"
#include "phong.h"
#include "csg.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

void object_test(Ray ray, Object *objects, Hit &best_hit)
{
	Object *obj = objects;

	best_hit.flag = false;

	while(obj != 0)
	{
		Hit obj_hit;
		obj_hit.flag=false;

		obj->intersection(ray, obj_hit);

		if (obj_hit.flag)
		{
			if (obj_hit.t > 0.0f)
			{
				if (best_hit.flag == false)
				{
					best_hit = obj_hit;
				} else if (obj_hit.t < best_hit.t)
				{
					best_hit = obj_hit;
				}
			}
		}
		obj = obj->next;
	}
	return;
}

void raytrace(Ray ray, Object *objects, Light *lights, Colour &colour, float &depth, int ref_limit) //gives depth and colour of hit point
{
	// first step, find the closest primitive

	Hit shadow_hit;
	Hit best_hit;
	object_test(ray, objects, best_hit);

	// if we found a primitive then compute the colour we should see
	if(best_hit.flag)
	{
		if(best_hit.what->material->eta != 0)
				{
					std::cout << "hit point eta: " << best_hit.what->material->eta << std::endl;
				}
		best_hit.what->material->compute_base_colour(colour);
		depth = best_hit.t;
		Light *light = lights;
		int no_of_reflections = 0; //depth counter to stop infinite recursion

		while (light != (Light *)0) //there is a light source
		{
			Vector viewer;
			Vector ldir;

			viewer.x = -best_hit.position.x; //viewing ray position
			viewer.y = -best_hit.position.y;
			viewer.z = -best_hit.position.z;
			viewer.normalise();

			bool lit;
			lit = light->get_direction(best_hit.position, ldir);

			if(ldir.dot(best_hit.normal)>0)
			{
				std::cout << "shadow: " << best_hit.what->material->eta << std::endl;

				lit=false;//light is facing wrong way.
			}

			if(lit)
			{
				Ray shadow_ray;

				shadow_ray.direction.x = -ldir.x;
				shadow_ray.direction.y = -ldir.y;
				shadow_ray.direction.z = -ldir.z;
				shadow_ray.position.x = best_hit.position.x + (0.0001f * shadow_ray.direction.x);
				shadow_ray.position.y = best_hit.position.y + (0.0001f * shadow_ray.direction.y);
				shadow_ray.position.z = best_hit.position.z + (0.0001f * shadow_ray.direction.z);

				object_test(shadow_ray, objects, shadow_hit);

				if(shadow_hit.flag==true)
				{
					if (shadow_hit.t < 1000000000.0f)
					{
						std::cout << "shadow: " << best_hit.what->material->eta << std::endl;

						lit = false; //there's a shadow so no lighting, if realistically close
					}
				}
			}

			if (lit) //do colour
			{

				std::cout << "lit2: " << best_hit.what->material->eta << std::endl;

				Colour intensity;
				Colour scaling;

				light->get_intensity(best_hit.position, scaling);

				best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);

				intensity.scale(scaling); //scale by scaling to get intensity

				colour.add(intensity); //add intensity to the colour
			}

			light = light->next; //move on to next light
		}

		// TODO: compute refraction ray if material supports it.
		if(best_hit.what->material->transparent)
		{
			Ray T;
			Colour trans_colour;
			Vector n;
			float trans_depth, eta, cos_theta_t, cos_theta_i, rpar, rper, kr, kt;

			//get refractive index for visible object
			eta = best_hit.what->material->eta;

			n =  best_hit.normal;
			cos_theta_i = n.dot(-1*ray.direction); //cosine of the angle of incidence
			cos_theta_t = sqrt(1 - ( 1 / (pow(eta,2) ) ) * ( 1 - pow(cos_theta_i,2) ) ); //cosine of the angle of refraction

			//test for total internal reflection
			if(cos_theta_t >= 0)
			{
				//get the refracted ray
				T.direction = (1/eta)*ray.direction - (cos_theta_t - (1/eta)*cos_theta_i)*best_hit.normal;

				T.position.x = best_hit.position.x + 0.001*T.direction.x;
				T.position.y = best_hit.position.y + 0.001*T.direction.y;
				T.position.z = best_hit.position.z + 0.001*T.direction.z;

				//raytrace the refracted ray
				raytrace(T, objects, lights, trans_colour, trans_depth, ref_limit);

				//fresnel equations
				rpar = (eta*cos_theta_i - cos_theta_t)/(eta*cos_theta_i + cos_theta_t);
				rper = (cos_theta_i - eta*cos_theta_t)/(cos_theta_i + eta*cos_theta_t);
				kr = 0.5*( pow(rpar,2) + pow(rper,2) );
				kt = 1-kr;

				//assign the calculated values to the object - kr may be used in reflection.
				best_hit.what->material->kr = kr;
				best_hit.what->material->kt = kt;

				//scale trans_colour by kt
				trans_colour.r *= kt;
				trans_colour.g *= kt;
				trans_colour.b *= kt;

				//add the refracted colour to the colour
				colour.add(trans_colour);
			}

		}
		// TODO: compute reflection ray if material supports it.
		if(best_hit.what->material->reflective) //material supports reflection
		{
			Ray R;
			float kr, ref_depth;
			Hit rhit;
			Colour ref_colour;

			//calculate reflected eye ray
			R.direction = ray.direction - 2*((best_hit.normal).dot(ray.direction))*best_hit.normal; //ray = viewer ray,

			R.position.x = best_hit.position.x + 0.001*R.direction.x;
			R.position.y = best_hit.position.y + 0.001*R.direction.y;
			R.position.z = best_hit.position.z + 0.001*R.direction.z;

			//tests if there is anything to be reflected in the surface
			object_test(R, objects, rhit);

			ref_limit -=1;

			if(ref_limit<0)
			{
				return;
			}

			//if there is something to be reflected, calculate the new colour of the point.
			if(rhit.flag)
			{
				//get colour of reflected object
				raytrace(R, objects, lights, ref_colour, ref_depth, ref_limit);

				//get reflection coefficient
				kr = best_hit.what->material->kr;

				//scale by kr
				ref_colour.r *= kr;
				ref_colour.g *= kr;
				ref_colour.b *= kr;

				//add reflection colour to the main object's colour.
				colour.add(ref_colour);
			}
		}
	}
	else
	{
		depth = 7.0f;
		colour.r = 0.0f;
		colour.g = 0.0f;
		colour.b = 0.0f;
	}


	std::cout << "r: " << colour.r << std::endl;

}



int main(int argc, char *argv[])
{
	int width = 128;
	int height = 128;
	// Create a framebuffer
	FrameBuffer *fb = new FrameBuffer(width,height);

	// The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
	Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
			0.0f, 0.0f, 1.0f, -2.7f,
			0.0f, 1.0f, 0.0f, 5.0f,
			0.0f, 0.0f, 0.0f, 1.0f);

	//  Read in the teapot model.
	PolyMesh *pm = new PolyMesh((char *)"teapot_smaller.ply", transform);

	Transform *t = new Transform();

	PolyMesh *fl = new PolyMesh((char *)"square.ply", t);
	PolyMesh *ce = new PolyMesh((char *)"ceiling.ply", t);
	PolyMesh *w1 = new PolyMesh((char *)"wall1.ply", t);
	PolyMesh *w2 = new PolyMesh((char *)"wall2.ply", t);
	PolyMesh *w3 = new PolyMesh((char *)"wall3.ply", t);

	Vertex v;
	v.x = 0.0f;
	v.y = 0.0f;
	v.z = 5.0f;

	Sphere *sphere = new Sphere(v, 1.5f);

	Vertex v2;
	v2.x = 0.0f; //x axis goes across
	v2.y = 1.0f; //y axis goes up and down
	v2.z = 5.0f; //z axis is the depth

	Sphere *sphere2 = new Sphere(v2,3.0f); //put in another sphere

	Ray ray;

	ray.position.x = 0.0001f;
	ray.position.y = 0.0f;
	ray.position.z = 0.0f;

	DirectionalLight *dl = new DirectionalLight(Vector(1.01f, -1.0f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source

	Phong bp1; //creates Phong surface illumination model for polymesh

	bp1.ambient.r = 0.2f;
	bp1.ambient.g = 0.0f;
	bp1.ambient.b = 0.0f;
	bp1.diffuse.r = 0.2f;
	bp1.diffuse.g = 0.0f;
	bp1.diffuse.b = 0.0f;
	bp1.specular.r = 0.4f;
	bp1.specular.g = 0.4f;
	bp1.specular.b = 0.4f;
	bp1.power = 40.0f;

	pm->material = &bp1;

	pm->material->transparent = false;

	pm->material->reflective = false;

	pm->material->eta = 0;

	pm->material->kr = 0.2;

	Phong bp2;

	bp2.ambient.r = 0.0f;
	bp2.ambient.g = 0.0f;
	bp2.ambient.b = 0.2f;
	bp2.diffuse.r = 0.0f;
	bp2.diffuse.g = 0.0f;
	bp2.diffuse.b = 0.3f;
	bp2.specular.r = 0.4f;
	bp2.specular.g = 0.4f;
	bp2.specular.b = 0.4f;
	bp2.power = 40.0f;



	sphere->material = &bp1;

	sphere->material->transparent = true;

	sphere->material->reflective = true;

	sphere->material->eta = 1.33;

	//sphere->material->kr = 0;

	//sphere->material->kt = 0;

	sphere2->material = &bp1;

	sphere2->material->transparent = true;

	sphere2->material->reflective = true;

	sphere2->material->eta = 1.33;

	//sphere2->material->kr = 0;

	//sphere2->material->kt = 0;





	fl->material = &bp1;

	fl->material->transparent = false;

	fl->material->reflective = false;

	fl->material->eta = 0;

	fl->material->kr = 0.41;

	fl->material->kt = 0;


		ce->material = &bp1;

		ce->material->transparent = false;

		ce->material->reflective = false;

		ce->material->eta = 0;

		ce->material->kr = 0.42;

		ce->material->kt = 0;


		w1->material = &bp1;

		w1->material->transparent = false;

		w1->material->reflective = false;

		w1->material->eta = 0;

		w1->material->kr = 0.43;

		w1->material->kt = 0;



		w2->material = &bp1;

		w2->material->transparent = false;

		w2->material->reflective = false;

		w2->material->eta = 0;

		w2->material->kr = 0.44;

		w2->material->kt = 0;


		w3->material = &bp1;

		w3->material->transparent = false;

		w3->material->reflective = false;

		w3->material->eta = 0;

		w3->material->kr = 0.45;

		w3->material->kt = 0;



		CSG *c = new CSG(sphere->center,sphere->radius, sphere2->center, sphere2->radius, true, false, false); //creates a an object from the difference of two spheres

		c->material = &bp1;

		c->material->transparent = true;

		c->material->reflective = true;

		c->material->eta = 1.52;

		c->material->kr = 0.4;

		c->material->kt = 0;

	c->next = fl;
/*
	//pm->next = fl;

	fl->next = ce;

	ce->next = w1;

	w1->next = w2;

	//w2->next = w3;
*/

	for (int y = 0; y < height; y += 1)
	{
		for (int x = 0; x < width; x += 1) //for each pixel
		{
			float fx = (float)x/(float)width;
			float fy = (float)y/(float)height;

			Vector direction;

			ray.direction.x = (fx-0.5f);
			ray.direction.y = (0.5f-fy);
			ray.direction.z = 0.5f;
			ray.direction.normalise();

			Colour colour;
			float depth;

			raytrace(ray, c, dl, colour, depth, 4);

			fb->plotPixel(x, y, colour.r, colour.g, colour.b);
			fb->plotDepth(x,y, depth);

		}

		cerr << "*" << flush;
	}

	// Output the framebuffer.
	fb->writeRGBFile((char *)"test.ppm");
	//  fb->writeDepthFile((char *)"depth.ppm");
	return 0;

}
