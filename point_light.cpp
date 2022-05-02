/*
 * point_light.cpp
 *
 *  Created on: 26 Nov 2019
 *      Author: carolinemillan
 */

#include "point_light.h"
#include "vector.h"
#include <random>
#include <ctime>

PointLight::PointLight()
{
	Light();
}

PointLight::PointLight(Vertex pt, Colour col, Vector dir)
{
	Light();

	direction = dir;
	direction.normalise();
	intensity = col;
	point = pt;
}

PointLight::PointLight(Vertex pt, Colour col)
{
	Light();

	this->intensity = col;
	this->point = pt;

	srand (time(NULL)); //initialises random seed
	this->direction.x = (rand() % 100000)/100000;
	this->direction.y = (rand() % 100000)/100000;
	this->direction.z = (rand() % 100000)/100000;


}

bool PointLight::get_direction(Vertex &surface, Vector &dir) //unsure about this
{
	//dir = direction;

	if (direction.x == dir.x && direction.y == dir.y && direction.y == dir.y)
	{
	return true;
	}
	else
	{
		return false;
	}
}

void PointLight::get_intensity(Vector D, Colour &level)
{
	if (direction.x == 0 && direction.y == 0 && direction.z == 0)
	{
		level = intensity;
	}
	else
	{
		//a function of the angle between light direction and light to object vector
		//I_l= I_i*(L.D)^n

		//get angle between two vectors
		float	cos_theta, d;
		cos_theta = direction.dot(D);
		//cos_theta.normalise();

		float n = 40;

		if(cos_theta >= 90) //is this the right test? radians?
		{
			level.r = 0;
			level.g = 0;
			level.b = 0;
		}
		else
		{
			level.r = intensity.r * pow(direction.dot(D), n);
			level.g = intensity.g * pow(direction.dot(D), n);
			level.b = intensity.b * pow(direction.dot(D), n);
		}
	}
}



