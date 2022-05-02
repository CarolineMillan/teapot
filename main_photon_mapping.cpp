/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

//    g++ -o pmexecutable main_photon_mapping.cpp framebuffer.cpp polymesh.cpp directional_light.cpp point_light.cpp sphere.cpp phong.cpp kdtree.cpp -lm

#include "framebuffer.h"
#include "ray.h"
#include "hit.h"
#include "polymesh.h"
#include "sphere.h"
#include "light.h"
#include "directional_light.h"
#include "point_light.h"
#include "material.h"
#include "phong.h"
#include "vector.h"
#include "vertex.h"
#include "photon.h"
#include "kdtree.hpp"
#include "BRDFread.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

void object_test(Ray ray, Object *objects, Hit &best_hit) {
	Object *obj = objects;

	best_hit.flag = false;

	while (obj != 0) {
		Hit obj_hit;
		obj_hit.flag = false;

		obj->intersection(ray, obj_hit);

		if (obj_hit.flag) {
			if (obj_hit.t > 0.0f) {
				if (best_hit.flag == false) {
					best_hit = obj_hit;
				} else if (obj_hit.t < best_hit.t) {
					best_hit = obj_hit;
				}
			}
		}

		obj = obj->next;
	}

	return;
}

void raytrace(Ray ray, Object *objects, Light *lights, Colour &colour,
		float &depth, Hit &best_hit) //gives depth and colour of hit point
		{
	// first step, find the closest primitive

	Hit shadow_hit;
	//Hit best_hit;
	object_test(ray, objects, best_hit);

	// if we found a primitive then compute the colour we should see
	if (best_hit.flag) {
		best_hit.what->material->compute_base_colour(colour); //this gives 0?
		depth = best_hit.t;
		Light *light = lights;
		//int no_of_reflections = 0; //depth counter to stop infinite recursion

		while (light != (Light*) 0) //there is a light source
		{
			Vector viewer;
			Vector ldir;

			viewer.x = -best_hit.position.x; //viewing ray position
			viewer.y = -best_hit.position.y;
			viewer.z = -best_hit.position.z;
			viewer.normalise();

			bool lit;
			lit = light->get_direction(best_hit.position, ldir); //need to write this function?

			if (ldir.dot(best_hit.normal) > 0) {
				lit = false; //light is facing wrong way.
			}

			if (lit) {

				Ray shadow_ray;

				shadow_ray.direction.x = -ldir.x;
				shadow_ray.direction.y = -ldir.y;
				shadow_ray.direction.z = -ldir.z;
				shadow_ray.position.x = best_hit.position.x
						+ (0.0001f * shadow_ray.direction.x);
				shadow_ray.position.y = best_hit.position.y
						+ (0.0001f * shadow_ray.direction.y);
				shadow_ray.position.z = best_hit.position.z
						+ (0.0001f * shadow_ray.direction.z);

				object_test(shadow_ray, objects, shadow_hit);

				if (shadow_hit.flag == true) {
					if (shadow_hit.t < 1000000000.0f) {
						lit = false; //there's a shadow so no lighting, if realistically close
					}
				}
			}

			if (lit) //do colour
			{
				Colour intensity;
				Colour scaling;

				light->get_intensity(best_hit.position, scaling);

				best_hit.what->material->compute_light_colour(viewer,
						best_hit.normal, ldir, intensity);

				intensity.scale(scaling); //scale by scaling to get intensity

				colour.add(intensity); //add intensity to the colour
			}

			light = light->next; //move on to next light. how to assign next light? Is that already done?
		}

		// TODO: compute refraction ray if material supports it.

		// TODO: compute reflection ray if material supports it.
	} else {
		depth = 7.0f;
		colour.r = 0.0f;
		colour.g = 0.0f;
		colour.b = 0.0f;
	}
}

void g_trace(Photon *photon, Object *objects, PointLight *pl, int ref_limit,
		Kdtree::KdNodeVector g_nodes) //traces a single photon through the scene, saves them as a vector of nodes
		{
	Ray photon_ray;
	photon->ray(photon_ray);

	Hit best_hit;

	object_test(photon_ray, objects, best_hit); //where does this photon hit first? maybe add in the shadow photons here

	photon->position = best_hit.position;

	Kdtree::KdNode p;

	std::vector<double> pos;

	pos[0] = photon->position.x;
	pos[1] = photon->position.y;
	pos[2] = photon->position.z;

	p.point = pos; //what point should i give it?
	p.data = photon;

	g_nodes.push_back(p);

	int new_w = 0;

	photon->g_russian_roulette(best_hit, photon->w, new_w);

	photon->w = new_w; //gets the new weight for russian_roulette

	while (!photon->absorbed) //russian roulette threshhold/weighting should mean the photon is eventually absorbed
	{
		if (photon->reflected) //calculate BDRF for reflection coefficient, use rendering equation later
		{

			Ray R;
			Hit rhit;

			//use BDRF to determine new intensity
			photon->BRDF_s = best_hit.what->material->BRDF_s;
			photon->BRDF_d = best_hit.what->material->BRDF_d;

			//add this photon to kd tree
			Kdtree::KdNode p;
			std::vector<double> pos;

			pos[0] = photon->position.x;
			pos[1] = photon->position.y;
			pos[2] = photon->position.z;

			p.point = pos;
			p.data = photon;

			g_nodes.push_back(p);
		}

		photon->g_russian_roulette(best_hit, photon->w, new_w); //need to do russian roulette again every iteration of the while loop

		photon->w = new_w;
	}

	//store last absorbed photon of the while loop in the tree - if statement?
	pos[0] = photon->position.x;
	pos[1] = photon->position.y;
	pos[2] = photon->position.z;

	p.point = pos;
	;
	p.data = photon;

	g_nodes.push_back(p);

	//continue path and store shadow photons
	while (objects != 0) {
		Hit obj_hit;
		obj_hit.flag = false;

		//does the photon intersect with this object?
		objects->intersection(photon_ray, obj_hit);

		if (obj_hit.flag) {
			if (obj_hit.position.x == best_hit.position.x
					&& obj_hit.position.y == best_hit.position.y
					&& obj_hit.position.y == best_hit.position.y) {
				photon->position = obj_hit.position;
				photon->shadow = true;

				Kdtree::KdNode p;
				std::vector<double> pos;

				pos[0] = photon->position.x;
				pos[1] = photon->position.y;
				pos[2] = photon->position.z;

				p.point = pos;
				p.data = photon;

				g_nodes.push_back(p);
			}
		}

	}
}


//doesn't emit towards specular surfaces, it emits in random directions then checks if the intersection is specular
void c_trace(Photon *photon, Object *objects, int ref_limit,
		Kdtree::KdNodeVector c_nodes) //traces a single photon through the scene, saves them as a vector of nodes
		{
	Ray photon_ray;
	photon->ray(photon_ray);

	Hit best_hit;

	object_test(photon_ray, objects, best_hit); //where does this photon hit first? maybe add in the shadow photons here

	if (best_hit.what->material->transparent) //only does it with photons that directly hit a transparent surface, what if they are reflected onto a transparent surface?
	{
		photon->position = best_hit.position;

		Kdtree::KdNode p;

		std::vector<double> pos;

		pos[0] = photon->position.x;
		pos[1] = photon->position.y;
		pos[2] = photon->position.z;

		p.point = pos;
		; //what point should i give it?
		p.data = photon;

		c_nodes.push_back(p);

		//continue path and store shadow photons - not in caustic

		int new_w;

		photon->c_russian_roulette(best_hit, photon->w, new_w);

		photon->w = new_w;

		while (!photon->absorbed) {
			if (photon->transmitted) {
				//direction is a function type of transmission - new direction is
				//if its transmitted we need a new direction, use raytracing transparency stuff
				Kdtree::KdNode p;

				std::vector<double> pos;

				pos[0] = photon->position.x;
				pos[1] = photon->position.y;
				pos[2] = photon->position.z;

				p.point = pos;
				p.data = photon;
				c_nodes.push_back(p);
				photon->c_russian_roulette(best_hit, photon->w, new_w);
				photon->w = new_w;
			}
		}

		//stores the last absorbed photon of the while loop

		pos[0] = photon->position.x;
		pos[1] = photon->position.y;
		pos[2] = photon->position.z;

		p.point = pos;
		p.data = photon;

		c_nodes.push_back(p);
	}
}


void g_photon_map(Object *objects, PointLight *lights, int ref_limit,
		Kdtree::KdNodeVector g_nodes) //creates a global photon map
		{

	Kdtree::KdNodeVector g_nodes1;

	for (int n = 0; n == 100000; n++) {
		while (lights != (Light*) 0) {
			Photon *photon = new Photon(*lights);

			//trace a photon through the scene and store all the hits in a kd tree
			g_trace(photon, objects, lights, 6, g_nodes1);

			for (int i = 0; i == g_nodes1.size(); i++) //add all these photons to the tree
					{
				Kdtree::KdNode kdnode;

				kdnode.data = g_nodes1[i].data;
				kdnode.point = g_nodes1[i].point;

				g_nodes.push_back(kdnode);
			}

			lights = lights->next; //move on to next light
		}
	}
		}


void c_photon_map(Object *objects, PointLight *lights, int ref_limit,
		Kdtree::KdNodeVector c_nodes) //creates a caustic photon map, only deals with transparency, don't need depth?
		{
	Kdtree::KdNodeVector c_nodes1;

	for (int n = 0; n == 100000; n++) //1 million photons?
			{
		while (lights != (Light*) 0) {
			Photon *photon = new Photon(*lights);

			//trace a photon through the scene and store all the hits in a kd tree - how to go through specular objects only?
			c_trace(photon, objects, 6, c_nodes1);

			for (int i = 0; i == c_nodes1.size(); i++) {
				c_nodes.push_back(c_nodes1[i]); //adds all photon hits to a tree
			}
			lights = lights->next;
		}
	}
}

int main(int argc, char *argv[]) {

	//Let's set the scene.
	int width = 1280;
	int height = 1280;
	// Create a framebuffer
	FrameBuffer *fb = new FrameBuffer(width, height);

	// The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
	Transform *transform = new Transform(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
			1.0f, -2.7f, 0.0f, 1.0f, 0.0f, 5.0f, 0.0f, 0.0f, 0.0f, 1.0f);

	Transform *t = new Transform();

	//  Read in the teapot model.
	PolyMesh *pm = new PolyMesh((char*) "teapot_smaller.ply", transform);

	// Read in the room model
	PolyMesh *fl = new PolyMesh((char*) "square.ply", t);
	PolyMesh *ce = new PolyMesh((char*) "ceiling.ply", t);
	PolyMesh *w1 = new PolyMesh((char*) "wall1.ply", t);
	PolyMesh *w2 = new PolyMesh((char*) "wall2.ply", t);
	PolyMesh *w3 = new PolyMesh((char*) "wall3.ply", t);

	Vertex v;
	v.x = 0.0f;
	v.y = 0.0f;
	v.z = 5.0f;

	Sphere *sphere = new Sphere(v, 1.0f);

	Vertex v2;
	v2.x = -1.0f;
	v2.y = 1.0f;
	v2.z = 3.0f;

	Sphere *sphere2 = new Sphere(v2, 0.5f); //put in another sphere

	Ray ray;

	ray.position.x = 0.0001f;
	ray.position.y = 0.0f;
	ray.position.z = 0.0f;

	//DirectionalLight *dl = new DirectionalLight(Vector(1.01f, -1.0f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f)); //creates a light source

	PointLight *pl = new PointLight(Vertex(0.0f, 3.0f, 5.0f), Colour(1.0f, 1.0f, 1.0f, 0.0f)); //a white light, above the teapot (hopefully)

	Phong bp1; //creates Phong surface illumination model for polymesh

	bp1.ambient.r = 0.0f;
	bp1.ambient.g = 0.0f;
	bp1.ambient.b = 0.2f;
	bp1.diffuse.r = 0.0f;
	bp1.diffuse.g = 0.0f;
	bp1.diffuse.b = 0.4f;
	bp1.specular.r = 0.4f;
	bp1.specular.g = 0.4f;
	bp1.specular.b = 0.4f;
	bp1.power = 40.0f;

	//defines what teapot looks like
	pm->material = &bp1;
	pm->material->transparent = true;
	pm->material->reflective = true;
	pm->material->eta = 1.33; //glass refractive index
	pm->material->kr = 0.7;

	pm->material->BRDF_s = bp1.specular;
	pm->material->BRDF_d = bp1.diffuse;

	Phong bp2;

	bp2.ambient.r = 0.0f;
	bp2.ambient.g = 0.2f;
	bp2.ambient.b = 0.0f;
	bp2.diffuse.r = 0.0f;
	bp2.diffuse.g = 0.4f;
	bp2.diffuse.b = 0.0f;
	bp2.specular.r = 0.4f;
	bp2.specular.g = 0.4f;
	bp2.specular.b = 0.4f;
	bp2.power = 40.0f;

	//defines what sphere looks like
	sphere->material = &bp2;
	sphere->material->transparent = false;
	sphere->material->reflective = true;
	sphere->material->kr = 0.4;

	sphere->material->BRDF_s = bp2.specular;
	sphere->material->BRDF_d = bp2.diffuse;

	//defines what sphere2 looks like
	sphere2->material = &bp2;
	sphere2->material->transparent = false;
	sphere2->material->reflective = true;
	sphere->material->kr = 0.8;
	sphere2->material->eta = 1.33; //refractive index for water

	sphere2->material->BRDF_s = bp2.specular;
	sphere2->material->BRDF_d = bp2.diffuse;

	Phong bp3; //pale green

	bp3.ambient.r = 0.0f;
	bp3.ambient.g = 0.01f;
	bp3.ambient.b = 0.0f;
	bp3.diffuse.r = 0.8f;
	bp3.diffuse.g = 1.0f;
	bp3.diffuse.b = 0.8f;
	bp3.specular.r = 0.2f;
	bp3.specular.g = 0.2f;
	bp3.specular.b = 0.2f;
	bp3.power = 40.0f;

	//defines what floor looks like
	fl->material = &bp2;
	fl->material->transparent = false;
	fl->material->reflective = true;
	fl->material->eta = 0;
	fl->material->kr = 0.41;
	fl->material->kt = 0;

	fl->material->BRDF_s = bp2.specular;
	fl->material->BRDF_d = bp2.diffuse;

	//defines what ceiling looks like
	ce->material = &bp2;
	ce->material->transparent = false;
	ce->material->reflective = true;
	ce->material->eta = 0;
	ce->material->kr = 0.42;
	ce->material->kt = 0;

	ce->material->BRDF_s = bp2.specular;
	ce->material->BRDF_d = bp2.diffuse;

	//defines what walls look like
	w1->material = &bp2;
	w1->material->transparent = false;
	w1->material->reflective = true;
	w1->material->eta = 0;
	w1->material->kr = 0.43;
	w1->material->kt = 0;

	w1->material->BRDF_s = bp2.specular;
	w1->material->BRDF_d = bp2.diffuse;

	w2->material = &bp2;
	w2->material->transparent = false;
	w2->material->reflective = true;
	w2->material->eta = 0;
	w2->material->kr = 0.44;
	w2->material->kt = 0;

	w2->material->BRDF_s = bp2.specular;
	w2->material->BRDF_d = bp2.diffuse;

	w3->material = &bp2;
	w3->material->transparent = false;
	w3->material->reflective = true;
	w3->material->eta = 0;
	w3->material->kr = 0.45;
	w3->material->kt = 0;

	w3->material->BRDF_s = bp2.specular;
	w3->material->BRDF_d = bp2.diffuse;


	//create a list of the objects
	sphere->next = pm;
	pm->next = sphere2;
	sphere2->next = fl;
	fl->next = ce;
	ce->next = w1;
	w1->next = w2;
	w2->next = w3;

	//create a list of bounding objects

	sphere->c_next = sphere2;

	//use a combination of raytracing and photon mapping to render
	for (int y = 0; y < height; y += 1) {
		for (int x = 0; x < width; x += 1) //for each pixel
				{
			float fx = (float) x / (float) width;
			float fy = (float) y / (float) height;

			Vector direction;

			ray.direction.x = (fx - 0.5f);
			ray.direction.y = (0.5f - fy);
			ray.direction.z = 0.5f;
			ray.direction.normalise();

			//Raytrace to get closest object

			Colour colour, L_l, L_s, L_c, L_d, L_0, L_e;
			float depth;
			Hit best_hit;

			object_test(ray, sphere, best_hit);

			Vector L_i, V;

			L_e.r = 0;
			L_e.g = 0;
			L_e.b = 0;

			//is the hit point a light source?
			while (pl != (Light*) 0) {

				if (pl->point.x == best_hit.position.x
						&& pl->point.y == best_hit.position.y
						&& pl->point.z == best_hit.position.z) //if the hit position is on a light, then it emits radiance
								{
					L_e.r += pl->intensity.r; //addition in case two lights are in the same position (unlikely)
					L_e.g += pl->intensity.g;
					L_e.b += pl->intensity.b;
				}
				pl = pl->next;
			}

			L_0 = L_e;

			//DIRECT - L_l
			raytrace(ray, pm, pl, L_l, depth, best_hit);

			L_l.r *= (best_hit.what->material->BRDF_s.r + best_hit.what->material->BRDF_d.r) * (-1* (pl->direction).dot(best_hit.normal));
			L_l.g *= (best_hit.what->material->BRDF_s.g + best_hit.what->material->BRDF_d.g) * (-1* (pl->direction).dot(best_hit.normal));
			L_l.b *= (best_hit.what->material->BRDF_s.b + best_hit.what->material->BRDF_d.b) * (-1* (pl->direction).dot(best_hit.normal));

			//SPECULAR - L_s
			if (best_hit.what->material->reflective) //material supports reflection
			{
				Ray R;
				float kr, ref_depth;
				Hit rhit;

				//calculate reflected eye ray
				R.direction = ray.direction
						- 2 * ((best_hit.normal).dot(ray.direction))
								* best_hit.normal;

				R.position.x = best_hit.position.x + 0.001 * R.direction.x;
				R.position.y = best_hit.position.y + 0.001 * R.direction.y;
				R.position.z = best_hit.position.z + 0.001 * R.direction.z;

				object_test(R, pm, rhit); //tests if there is an object to be reflected in the surface

				if (rhit.flag)
				{
					Hit ref_hit;
					raytrace(R, pm, pl, L_s, ref_depth, ref_hit);

					L_s.r *= best_hit.what->material->BRDF_s.r * (-1* (pl->direction).dot(best_hit.normal));
					L_s.g *= best_hit.what->material->BRDF_s.g * (-1* (pl->direction).dot(best_hit.normal));
					L_s.b *= best_hit.what->material->BRDF_s.b * (-1* (pl->direction).dot(best_hit.normal));
				}
			} else {
				L_s.r = 0;
				L_s.g = 0;
				L_s.b = 0;
			}

			//CAUSTICS - L_c

			//create a photon map
			Kdtree::KdNodeVector nodes;
			c_photon_map(sphere, pl, 6, nodes);

			//create kd tree
			Kdtree::KdTree c_tree(&nodes);

			Kdtree::KdNodeVector c_nodes_knn;
			std::vector<double> v;

			v[0] = best_hit.position.x;
			v[1] = best_hit.position.y;
			v[2] = best_hit.position.z;

			const std::vector<double> v0 = v;

			c_tree.k_nearest_neighbors(v0, 50, c_nodes_knn);

			float r; //get radius of the sphere from max distance between best_hit.position and rest of nodes
			float a, b, c, d, max_d = 0;

			for (int i = 0; i == c_nodes_knn.size(); i++) //finds max distance
					{
				a = best_hit.position.x - c_nodes_knn[i].point[0];
				b = best_hit.position.y - c_nodes_knn[i].point[1];
				c = best_hit.position.z - c_nodes_knn[i].point[2];
				d = sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2));

				if (d > max_d)
				{
					max_d = d;
				}
			}

			r = max_d;

			float dA = M_PI * pow(r, 2); //area of disc associated with sphere of radius r

			for (int i = 0; i == 49; i++) //each photon in sphere
					{
				Photon *p = new Photon();
				p = c_nodes_knn[i].data;

				//divide flux that photon represents by area of s and multiply by photon's BDRF to get radiance here
				L_c.r += p->BRDF_d.r * (p->intensity.r) * (-1* (p->direction).dot(best_hit.normal));
				L_c.g += p->BRDF_d.g * (p->intensity.g) * (-1* (p->direction).dot(best_hit.normal));
				L_c.b += p->BRDF_d.b * (p->intensity.b) * (-1* (p->direction).dot(best_hit.normal));
			}
			L_c.r *= 1/dA;
			L_c.g *= 1/dA;
			L_c.b *= 1/dA;

			//SOFT INDIRECT - L_d

			//create a photon map
			Kdtree::KdNodeVector nodes1;
			g_photon_map(sphere, pl, 6, nodes1);

			//create kd tree
			Kdtree::KdTree g_tree(&nodes1);

			Kdtree::KdNodeVector g_nodes_knn;

			g_tree.k_nearest_neighbors(v0, 50, g_nodes_knn);

			for (int i = 0; i == g_nodes_knn.size(); i++) //finds max distance
					{
				a = best_hit.position.x - g_nodes_knn[i].point[0];
				b = best_hit.position.y - g_nodes_knn[i].point[1];
				c = best_hit.position.z - g_nodes_knn[i].point[2];
				d = sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2));

				if (d > max_d) {
					max_d = d;
				}
			}

			r = max_d;

			dA = M_PI * pow(r, 2); //area of disc associated with sphere

			//find where the photons came from
			vector<Ray> nxt_rays;

			for (int i = 0; i == 49; i++) {
				Photon *p;

				p = g_nodes_knn[i].data;

				L_d.r += p->BRDF_d.r * (p->intensity.r) * (-1* (p->direction).dot(best_hit.normal));
				L_d.g += p->BRDF_d.g * (p->intensity.g) * (-1* (p->direction).dot(best_hit.normal));
				L_d.b += p->BRDF_d.b * (p->intensity.b) * (-1* (p->direction).dot(best_hit.normal));

				nxt_rays[i].position.x = g_nodes_knn[i].point[0];
				nxt_rays[i].position.y = g_nodes_knn[i].point[1];
				nxt_rays[i].position.z = g_nodes_knn[i].point[2];
				nxt_rays[i].direction = p->direction; //direction the photon came from
			}

			L_d.r *= 1/dA;
			L_d.g *= 1/dA;
			L_d.b *= 1/dA;

			vector<float> kr;

			for (int k = 0; k == nxt_rays.size(); k++) {
				kr[k] = best_hit.what->material->kr; //BRDF?
			}
			//store the intensity at this point, multiply by kr of the object

			//use importance sampling to decide which to follow and then resample at the next hit point (reflection)

			//bounce limit for photons
			int ref_limit = 3;

			Colour col, L_g_ref; //global reflections colour, will then be added onto L_d.
			float g_depth;
			Hit nxt_hit;

			while (ref_limit > 0) {
				for (int k = 0; k == nxt_rays.size(); k++) //ideally for imp_rays
				{
					object_test(nxt_rays[k], pm, nxt_hit); //get next hit point of each ray

					//do k nearest neighbours sampling again, but scale each reflection by kr (BRDF?) and add them all together
					Kdtree::KdNodeVector g_nodes_knn_1;

					std::vector<double> v1;

					v1[0] = nxt_hit.position.x;
					v1[1] = nxt_hit.position.y;
					v1[2] = nxt_hit.position.z;

					const std::vector<double> v2;

					g_tree.k_nearest_neighbors(v2, 50, g_nodes_knn_1);

					float r; //get radius of the sphere from max distance between best_hit.position and rest of nodes
					float a, b, c, d, max_d = 0;

					for (int i = 0; i == g_nodes_knn_1.size(); i++) //finds max distance
					{

						a = best_hit.position.x - g_nodes_knn_1[i].point[0];
						b = best_hit.position.y - g_nodes_knn_1[i].point[1];
						c = best_hit.position.z - g_nodes_knn_1[i].point[2];
						d = sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2));

						if (d > max_d)
						{
							max_d = d;
						}
					}

					r = max_d;

					float dA;
					dA = M_PI * pow(r, 2); //area of disc associated with sphere

					//now work out average colour
					for (int i = 0; i == 49; i++)
					{
						Photon *p;

						p = g_nodes_knn_1[i].data;

						col.r += kr[k] * p->BRDF * (p->intensity.r) * (-1* (p->direction).dot(best_hit.normal));
						col.g += kr[k] * p->BRDF * (p->intensity.g) * (-1* (p->direction).dot(best_hit.normal));
						col.b += kr[k] * p->BRDF * (p->intensity.b) * (-1* (p->direction).dot(best_hit.normal));

						nxt_rays[i].position.x = g_nodes_knn_1[i].point[0];
						nxt_rays[i].position.y = g_nodes_knn_1[i].point[1];
						nxt_rays[i].position.z = g_nodes_knn_1[i].point[2];
						nxt_rays[i].direction = p->direction; //direction the photon came from
					}

					col.r *= 1/dA;
					col.g *= 1/dA;
					col.b *= 1/dA;

					//get next kr value
					kr[k] = nxt_hit.what->material->kr;

					ref_limit -= 1;
				}
				L_d.r += col.r;
				L_d.g += col.g;
				L_d.b += col.b;
			}
			//TOTAL
			colour.r = L_e.r + L_l.r + L_s.r + L_c.r + L_d.r;
			colour.g = L_e.g + L_l.g + L_s.g + L_c.g + L_d.g;
			colour.b = L_e.b + L_l.b + L_s.b + L_c.b + L_d.b;
			//plot the point
			fb->plotPixel(x, y, colour.r, colour.g, colour.b);
			fb->plotDepth(x, y, depth);

		}

		cerr << "*" << flush;
	}

	// Output the framebuffer.
	fb->writeRGBFile((char*) "test.ppm");
	//  fb->writeDepthFile((char *)"depth.ppm");
	return 0;

}
