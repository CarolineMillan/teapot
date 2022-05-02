/*
 * photon_map.cpp
 *
 *  Created on: 24 Nov 2019
 *      Author: carolinemillan
 */


//SCALE THE COLOUR SO ITS BETWEEN 0 AND 1
/* TODO:
 * - barycentric coordinates to interpolate normals to get a smooth surface (calculate normals at shared vertices of triangles and interpolate these)
 * - scale the colour at some point t keep it between 0 and 1
 */

#include "photon.h"
#include "vertex.h"
#include "colour.h"
#include "vector.h"
#include "point_light.h"
#include "object.h"
#include <random>
#include <ctime>

Photon::Photon()
{
	position = (0,0,0);

	direction.x = 0;
	direction.y = 0;
	direction.z = 0;

	intensity.r = 0;
	intensity.g = 0;
	intensity.b = 0;

	shadow = false;
	reflected = false;
	absorbed = false;
	transmitted = false;

	w = 0;

	BRDF = BRDF_s+BRDF_d;

}

Photon::Photon(PointLight L)
{
	//generate a random photon at light source
	this->start_pt = L.point;
	this->intensity = L.intensity; //gets the colour of the photon
	this->direction = L.direction;

	//generate a random direction

	while (intensity == 0)
	{
		//generate random direction until you get one that's lit

		//generate random point and check its outside the sphere - x*x + y*y + z*z > 0 (?)
/*
		srand (time(NULL)); //initialises random seed
		direction.x = (rand() % 100000)/100000;
		direction.y = (rand() % 100000)/100000;
		direction.z = (rand() % 100000)/100000;
*/
		L.get_intensity(direction, intensity);
	}

	shadow = false;
	reflected = false;
	absorbed = false;
	transmitted = false;

	//from https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution
	srand (time(NULL)); //initialises random seed
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> dis(1, 100);

	w = dis(gen)/100;

	BRDF = BRDF_s + BRDF_d;
}

// the ray generated in the following is used to see where the photon intersects the scene

void Photon::ray(Ray &r)
{
	//our photon has a starting position and direction, so we can generate a ray
	r.position = start_pt;
	r.direction = direction;
}

// russian roulette decides whether the photon is absorbed, transmitted or reflected

//russian roulette on rgb photons seperately
void Photon::g_russian_roulette(Hit h, int w, int &new_w)
{
	if(h.flag)
	{
		float kr, kt, P;
		kr = h.what->material->kr;
		kt = h.what->material->kt;
		P = 1-kr;

		if(w <= 90)
		{
			//sample s uniformly
			//used https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution for this.
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> dis(1, 100);
			float s = dis(gen)/100;

			if(s<P)//terminate path
			{
				absorbed = true;
				reflected = false;
				transmitted = false;
				new_w = w;
			}
			else //its successful
			{
				new_w = w/(1-P);
				reflected = true;
				transmitted = false;
				absorbed = false;
			}
		}
		else //absorbed. i.e. no more photon to deal with
		{
			absorbed = true;
			reflected = false;
			transmitted = false;
		}
	}
}

void Photon::c_russian_roulette(Hit h, int w, int &new_w)
{
	if(h.flag)
	{
		float kt, P;
		kt = h.what->material->kt;
		P = 1-kt;

		if(w <= 90)
		{
			//sample s uniformly
			//used https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution for this.
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> dis(1, 100);
			float s = dis(gen)/100;

			if(s<P)//P any probability less than 1 - this is dependant on the surface.
			{
				//terminate path
				absorbed = true;
				reflected = false;
				transmitted = false;

			}
			else
			{
				new_w = w/(1-P);
				transmitted = true;
				reflected = false;
				absorbed = false;
			}
		}
		else //absorbed. i.e. no more photon to deal with
		{
			absorbed = true;
			reflected = false;
			transmitted = false;
		}
	}
}
