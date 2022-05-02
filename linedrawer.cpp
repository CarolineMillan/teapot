/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the code you will need to replace for Lab 1.
 *
 * It contains two simple implementations that loop over the longest axis adding the gradient to the position on the other axis at each step.
 * The objective is for you to remove the use of floating point variables and work entirely with integers.
 * You should use Bresenhams Line Algorithm to achieve this.
 */




#include <iostream>
#include "linedrawer.h"
// defines set of x cordinates
int draw_x_line(FrameBuffer *fb, int sx, int sy, int ex, int ey) //takes four integers - sx = start x, ex = end x
{
  int dir = 1; // define dirivative = 1
  if (sx > ex) // if start x is more than end x then dir = -1
  {
    dir = -1;
  }

  int   x     = sx; // renames sx to just integer x
  float y     = (float)sy; //renames sy to just a float y
  float slope = ((float)ey-(float)sy)/((float)ex-(float)sx); // defines slope as a float dy/dx
        slope = slope * dir;  //gets direction of slope right, it should do this automatically by considering dx and dy seperately. I dont think you need dir
  
  while (x != ex) // x not equal to end x, plot the pixels and then add dy/dx to y and +/- 1 to x
  {
    fb->plotPixel(x, (int)y, 1.0f, 1.0f, 1.0f);

    y += slope;

    x += dir;
  }

}
// defines set of y coordinates
int draw_y_line(FrameBuffer *fb, int sx, int sy, int ex, int ey)
{
  int dir = 1;
  if (sy > ey)
  {
    dir = -1;
  }

  int   y     = sy;
  float x     = (float)sx;
  float slope = ((float)ex-(float)sx)/((float)ey-(float)sy);
        slope = slope * dir;
  
  while (y != ey)
  {
    fb->plotPixel((int)x, y, 1.0f, 1.0f, 1.0f);

    x += slope;

    y += dir;
  }

}


int draw_line(FrameBuffer *fb, int sx, int sy, int ex, int ey)
{
  if ((sx == ex) && (sy==ey))
  {
    return fb->plotPixel(sx, sy, 1.0f, 1.0f, 1.0f);
    
  } else if (((ex-sx)* (ex-sx)) >= ((ey-sy)* (ey-sy)))
  {
    return draw_x_line(fb, sx, sy, ex, ey);
    
  } else
  {
    return draw_y_line(fb, sx, sy, ex, ey);
  }
}


