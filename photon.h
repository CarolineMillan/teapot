/*
 * photon_map.h
 *
 *  Created on: 24 Nov 2019
 *      Author: carolinemillan
 */

#ifndef PHOTON_H_
#define PHOTON_H_

#pragma once
#include "vertex.h"
#include "vector.h"
#include "colour.h"
#include "ray.h"
#include "hit.h"
#include "point_light.h"
//#include "photon.h"

class Photon {
public:

	Vertex start_pt;
	Vertex position; //position of intersection
	Vector direction;
	Colour intensity;

	float w; //for the weight for russian roulette
	float BRDF;

	Colour BRDF_s, BRDF_d;

	bool shadow;
	bool reflected;
	bool absorbed;
	bool transmitted;


	Photon();
	Photon(PointLight L); //generate a photon at the light source
	void ray(Ray &r); //generates a ray to model the photon
	void g_russian_roulette(Hit h, int w, int &new_w);
	void c_russian_roulette(Hit h, int w, int &new_w);
	void outcomes(); //unsure on this one - delete
	void main(); //delete
};



#endif /* PHOTON_H_ */

