/*
 * csg.h
 *
 *  Created on: 6 Dec 2019
 *      Author: carolinemillan
 *
 *
 *      CSG for spheres
 */

#include "object.h"
#include "sphere.h"
#include "hit.h"

#ifndef CSG_H_
#define CSG_H_

//problem with defining object


class CSG : public Object {
public:
	Vertex v1, v2;
	float r1, r2;
	bool diff, uni, inter;

	CSG()
	{
		v1 = Vertex (0,0,0);
		r1 = 0;
		v2 = Vertex (0,0,0);
		r2 = 0;

		//give boolean values to implement more complete CSG
		diff = false;
		uni = false;
		inter = false;
	}

	CSG(Vertex v1, float r1, Vertex v2, float r2, bool diff, bool uni, bool inter) //try initialising it as this and see if it works
	{
		this->v1 = v1;
		this->r1 = r1;
		this->v2 = v2;
		this->r2 = r2;

		//decides what happens to these spheres. only one can be true at a time
		this->diff = diff;
		this->uni = uni;
		this->inter = inter;
	}

	//have another one that takes objects as an input, not vterices and floats, to do more layers of CSG.



public:
	virtual void intersection(Ray r, Hit h) //gives the intersection with the difference shape. maybe give a CSG three bools to determine whether to u/i/d.
	{

		std::cout << "a3" << std::endl;
		//s2-s1

		//try doing object_test, and take best_hit to be the closest s2 hit not contained in s1

		Sphere *s1 = new Sphere(v1, r1);
		Sphere *s2 = new Sphere(v2, r2);

		h.flag = false;

		Hit hit_2, hit_1, hit_3, hit_4;
		hit_2.flag=false;

		s2->intersection(r, hit_2);

		if (hit_2.flag) //ray intersects s2, did it intersect s1 first?
		{
			s1->intersection(r, hit_1);

			if(hit_1.flag) //ray also intersects s1
			{
				if (hit_2.t < hit_1.t) //if s2 hit is closer than s1 hit
				{
					h = hit_2; //s2 is closest hit, you can see s2
				}
				else //s1 is the closer hit
				{
					//travel all the way through the sphere and fire another ray to get other intersection with s1, is s2 intersection between these points?
					s1->intersection(r, hit_3);

					if(hit_3.flag) //reached the other side of the s1
					{
						if(hit_1.t < hit_2.t < hit_3.t)//if the s2 hit is inside s1
						{
							//check for the other end of s2
							r.position = hit_2.position;
							s2->intersection(r, hit_4);

							if(hit_4.flag)
							{
								if(hit_2.t < hit_4.t < hit_3.t) //if the whole s2 is within s1, we see nothing
								{
									h.flag = false;
								}
								else //csg_hit_4 is the other side of the sphere s2, so s2 becomes visible immediately after we leave s1
								{
									h.t = hit_3.t; //when you leave s1
									h.position.x = r.position.x + hit_3.t * r.direction.x;
									h.position.y = r.position.y + hit_3.t * r.direction.y;
									h.position.z = r.position.z + hit_3.t * r.direction.z;
									h.normal.x = h.position.x - s2->center.x;
									h.normal.y = h.position.y - s2->center.y;
									h.normal.z = h.position.z - s2->center.z;
									h.normal.normalise();
									h.flag = true;
									return;
								}
							}
							else
							{
								return;
							}
						}
						else //s2 hit is not within s1, and its not before s1, so it must be after s1 in which case we can see it
						{
							h = hit_2;
							return;
						}

					}

				}
			}
			else //only intersects s2, doesn't intersect s1 at all
			{
				h = hit_2;
				return;
			}
		}

		return;
	}
/*
	if(uni)
	{

	}

	if(inter)
	{

	}
*/

};


#endif /* CSG_H_ */
