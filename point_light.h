/*
 * point_light.h
 *
 *  Created on: 26 Nov 2019
 *      Author: carolinemillan
 */

#ifndef POINT_LIGHT_H_
#define POINT_LIGHT_H_

// PointLight is a child class of Light and implements a light
// that emmits from a single point.

#pragma once
#include "light.h"
#include "vertex.h"
#include "vector.h"

class PointLight : public Light {
public:

	Vertex point;
	Vector direction;
	Colour intensity;
	PointLight *next;

	PointLight();
	PointLight(Vertex pt, Colour col);
	PointLight(Vertex pt, Colour col, Vector dir);

	bool get_direction(Vertex &surface, Vector &dir); //don't need this?

	void get_intensity(Vector D, Colour &intensity);

};


#endif /* POINT_LIGHT_H_ */
