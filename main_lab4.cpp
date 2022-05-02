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



//    g++ -o lab4executable main_lab4.cpp framebuffer.cpp polymesh.cpp directional_light.cpp sphere.cpp phong.cpp -lm


#include "framebuffer.h"
#include "ray.h"
#include "hit.h"
#include "polymesh.h"
#include "sphere.h"
#include "light.h"
#include "directional_light.h"
//#include "point_light.h"
#include "material.h"
#include "phong.h"

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
		best_hit.what->material->compute_base_colour(colour); //this gives 0?
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
			lit = light->get_direction(best_hit.position, ldir); //is it hit by the light?

			if(ldir.dot(best_hit.normal)>0)
			{
				//std::cout << "shadow: " << best_hit.what->material->kr << std::endl;

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
						//std::cout << "shadow: " << best_hit.what->material->kr << std::endl;

						lit = false; //there's a shadow so no lighting, if realistically close
					}
				}
			}

			if (lit) //do colour
			{
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
			float trans_depth, eta, cos_theta_t, cos_theta_i, rpar, rper, kr, kt;

			eta = best_hit.what->material->eta;

			Vector n;
			n =  best_hit.normal;
			cos_theta_i = n.dot(-1*ray.direction);
			cos_theta_t = sqrt(1 - ( 1 / (pow(eta,2) ) ) * ( 1 - pow(cos_theta_i,2) ) );

			//tests for total internal reflection
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
				kr = 0.5*( pow(rpar,2) + pow(rper,2) ); //do we use this in reflection?
				kt = 1-kr;

				best_hit.what->material->kr = kr; //this is the value needed for reflection in this surface at this point
				best_hit.what->material->kt = kt;

				//scale trans_colour by kt
				trans_colour.r *= kt; //use scale function?
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
			R.direction = ray.direction - 2*((best_hit.normal).dot(ray.direction))*best_hit.normal;

			R.position.x = best_hit.position.x + 0.001*R.direction.x;
			R.position.y = best_hit.position.y + 0.001*R.direction.y;
			R.position.z = best_hit.position.z + 0.001*R.direction.z;

			//tests if there is an object to be reflected in the surface
			object_test(R, objects, rhit); 

			ref_limit -=1;

			if(ref_limit<0)
			{
				return;
			}

			if(rhit.flag)
			{
				//get colour of reflection 
				raytrace(R, objects, lights, ref_colour, ref_depth, ref_limit); 

				kr = best_hit.what->material->kr;

				//scale by kr
				ref_colour.r *= kr; 
				ref_colour.g *= kr;
				ref_colour.b *= kr;

				//add reflective colour to the colour
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
}



int main(int argc, char *argv[])
{
	int width = 1280;
	int height = 1280;
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
	PolyMesh *w4 = new PolyMesh((char *)"wall4.ply", t);


	Vertex v;
	v.x = 0.0f;
	v.y = 0.0f;
	v.z = -3.0f;

	Sphere *sphere = new Sphere(v, 2.5f);

	Vertex v2;
	v2.x = -2.0f; //x axis goes across
	v2.y = 1.0f; //y axis is the depth
	v2.z = 3.0f; //z axis goes up and down

	Sphere *sphere2 = new Sphere(v2,1.0f); //put in another sphere

	//sphere->next = pm;

	//should that be commented? ^

	Ray ray;

	ray.position.x = 0.0001f;
	ray.position.y = 0.0f;
	ray.position.z = 0.0f;
/*
	DirectionalLight *dl = new DirectionalLight(Vector(0.51f, 1.0f, 0.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source
	DirectionalLight *dl1 = new DirectionalLight(Vector(-0.51f, 1.0f, 0.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source
	DirectionalLight *dl2 = new DirectionalLight(Vector(0.51f, 1.0f, 0.5f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source
	DirectionalLight *dl3 = new DirectionalLight(Vector(0.51f, 1.0f, -0.5f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source
	DirectionalLight *dl4 = new DirectionalLight(Vector(-0.51f, 1.0f, 0.5f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source
	DirectionalLight *dl5 = new DirectionalLight(Vector(-0.51f, 1.0f, -0.5f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source
	DirectionalLight *dl6 = new DirectionalLight(Vector(0.01f, 1.0f, 0.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source
	DirectionalLight *dl7 = new DirectionalLight(Vector(0.01f, 1.0f, -0.5f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source
	DirectionalLight *dl8 = new DirectionalLight(Vector(0.01f, 1.0f, 0.5f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source

	dl->next = dl1;
	dl1->next = dl2;
	dl2->next = dl3;
	dl3->next = dl4;
	dl4->next = dl5;
	dl5->next = dl6;
	dl6->next = dl7;
	dl6->next = dl8;
*/
	//PointLight *pl = new PointLight(Vertex (0,2,5), Colour(1.0f,1.0f,1.0f,0.0f));


	DirectionalLight *dl = new DirectionalLight(Vector(0.01f, 0.0f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source

	Phong bp1; //creates Phong surface illumination model for polymesh

	bp1.ambient.r = 0.1f;
	bp1.ambient.g = 0.0f;
	bp1.ambient.b = 0.0f;
	bp1.diffuse.r = 0.1f;
	bp1.diffuse.g = 0.0f;
	bp1.diffuse.b = 0.0f;
	bp1.specular.r = 0.4f;
	bp1.specular.g = 0.4f;
	bp1.specular.b = 0.4f;
	bp1.power = 40.0f;




	Phong bp3; //creates Phong surface illumination model for polymesh

	bp3.ambient.r = 0.2f;
	bp3.ambient.g = 0.2f;
	bp3.ambient.b = 0.2f;
	bp3.diffuse.r = 0.1f;
	bp3.diffuse.g = 0.1f;
	bp3.diffuse.b = 0.1f;
	bp3.specular.r = 0.4f;
	bp3.specular.g = 0.4f;
	bp3.specular.b = 0.4f;
	bp3.power = 40.0f;


	Phong bp2;

	bp2.ambient.r = 0.0f;
	bp2.ambient.g = 0.1f;
	bp2.ambient.b = 0.0f;
	bp2.diffuse.r = 0.0f;
	bp2.diffuse.g = 0.1f;
	bp2.diffuse.b = 0.0f;
	bp2.specular.r = 0.4f;
	bp2.specular.g = 0.4f;
	bp2.specular.b = 0.4f;
	bp2.power = 40.0f;


	Phong bp4;

	bp4.ambient.r = 0.1f;
	bp4.ambient.g = 0.0f;
	bp4.ambient.b = 0.1f;
	bp4.diffuse.r = 0.1f;
	bp4.diffuse.g = 0.0f;
	bp4.diffuse.b = 0.1f;
	bp4.specular.r = 0.4f;
	bp4.specular.g = 0.4f;
	bp4.specular.b = 0.4f;
	bp4.power = 40.0f;


	Phong bp5;

	bp5.ambient.r = 0.0f;
	bp5.ambient.g = 0.0f;
	bp5.ambient.b = 0.1f;
	bp5.diffuse.r = 0.0f;
	bp5.diffuse.g = 0.0f;
	bp5.diffuse.b = 0.1f;
	bp5.specular.r = 0.4f;
	bp5.specular.g = 0.4f;
	bp5.specular.b = 0.4f;
	bp5.power = 40.0f;




	pm->material = &bp5; //what does the & mean?

	pm->material->transparent = true;

	pm->material->reflective = true;

	pm->material->eta = 2.4; //glass refractive index

	pm->material->kr = 0;



	sphere->material = &bp3;

	sphere->material->transparent = false;

	sphere->material->reflective = false;

	sphere->material->eta = 0;

	sphere->material->kr = 0.4;

	//sphere->material->kt = 0;


	sphere2->material = &bp3;

	sphere2->material->transparent = false;

	sphere2->material->reflective = false;

	sphere2->material->eta = 0;

	sphere2->material->kr = 0.4;

	//sphere2->material->kt = 0;


	fl->material = &bp3;

	fl->material->transparent = true;

	fl->material->reflective = true;

	fl->material->eta = 1.52;

	fl->material->kr = 0;

	fl->material->kt = 0;


	ce->material = &bp2;

	ce->material->transparent = true;

	ce->material->reflective = true;

	ce->material->eta = 1.52;

	ce->material->kr = 0;

	ce->material->kt = 0;


	w1->material = &bp3;

	w1->material->transparent = false;

	w1->material->reflective = true;

	w1->material->eta = 1.52;

	w1->material->kr = 0.4;

	w1->material->kt = 0;


	w2->material = &bp1;

	w2->material->transparent = false;

	w2->material->reflective = false;

	w2->material->eta = 1.52;

	w2->material->kr = 0.1;

	w2->material->kt = 0;


	w3->material = &bp1;

	w3->material->transparent = false;

	w3->material->reflective = false;

	w3->material->eta = 1.52;

	w3->material->kr = 0.1;

	w3->material->kt = 0;


	w4->material = &bp4;

	w4->material->transparent = true;

	w4->material->reflective = true;

	w4->material->eta = 1.0003;

	w4->material->kr = 0;

	w4->material->kt = 0;


	//create list of objects
	sphere->next = pm;

	pm->next = sphere2;

	sphere2->next = fl;

	fl->next = ce;

	ce->next = w1;

	w1->next = w2;

	w2->next = w3;

	w3->next = w4;


	for (int y = 0; y < height; y += 1)
	{
		for (int x = 0; x < width; x += 1) //for each pixel
		{
			float fx = (float)x/(float)width;
			float fy = (float)y/(float)height;

			Vector direction;

			ray.direction.x = (fx-0.5f); //what is this part doing?
			ray.direction.y = (0.5f-fy);
			ray.direction.z = 0.5f;
			ray.direction.normalise();

			Colour colour;
			float depth;

			raytrace(ray, sphere, dl, colour, depth, 4);

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


// - may want to use scene.cpp later on, it deals with lights AND objects whereas main_lab4 deals with just objects

/* etas
 * Vacuum: 1.0
		• Air 1.0003
		• Water 1.33
		• Glass 1.52
		• Diamond 2.4
 */
