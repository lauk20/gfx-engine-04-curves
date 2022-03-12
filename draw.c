#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"


/*======== void add_circle() ==========
  Inputs:   struct matrix * edges
            double cx
            double cy
            double r
            double step

  Adds the circle at (cx, cy) with radius r to edges
  ====================*/
void add_circle( struct matrix *edges,
                 double cx, double cy, double cz,
                 double r, double step ) {

  float i = 0;
  float x = r * cos(2 * M_PI * i) + cx;
  float y = r * sin(2 * M_PI * i) + cy;

  for (i = i + step; i <= 1; i = i + step){
    float x2 = r * cos(2 * M_PI * i) + cx;
    float y2 = r * sin(2 * M_PI * i) + cy;

    add_edge(edges, x, y, 0, x2, y2, 0);

    x = x2;
    y = y2;
  }
}

/*======== void add_curve() ==========
Inputs:   struct matrix *edges
         double x0
         double y0
         double x1
         double y1
         double x2
         double y2
         double x3
         double y3
         double step
         int type

Adds the curve bounded by the 4 points passsed as parameters
of type specified in type (see matrix.h for curve type constants)
to the matrix edges
====================*/
void add_curve( struct matrix *edges,
                double x0, double y0,
                double x1, double y1,
                double x2, double y2,
                double x3, double y3,
                double step, int type ) {

  if (type == HERMITE){
    struct matrix * hermite_inverse = make_hermite();
    struct matrix * givenX = generate_curve_coefs(x0, x1, x2, x3, HERMITE);
    struct matrix * givenY = generate_curve_coefs(y0, y1, y2, y3, HERMITE);
    /*
    //X
    given->m[0][0] = x0; //x0
    given->m[1][0] = x1; //x1
    given->m[2][0] = x2; //rx0
    given->m[3][0] = x3; //rx1
    //Y
    given->m[0][1] = y0; //y0
    given->m[1][1] = y1; //y1
    given->m[2][1] = y2; //ry0
    given->m[3][1] = y3; //ry1

    given->lastcol = 2;
    */

    matrix_mult(hermite_inverse, givenX);
    matrix_mult(hermite_inverse, givenY);

    float i = 0;
    float x = givenX->m[0][0] * (i * i * i) + givenX->m[1][0] * (i * i) + givenX->m[2][0] * (i) + givenX->m[3][0];
    float y = givenY->m[0][1] * (i * i * i) + givenY->m[1][1] * (i * i) + givenY->m[2][1] * (i) + givenY->m[3][1];

    for (i = i + step; i <= 1; i = i + step){;
      float x1 = givenX->m[0][0] * (i * i * i) + givenX->m[1][0] * (i * i) + givenX->m[2][0] * (i) + givenX->m[3][0];
      float y1 = givenY->m[0][1] * (i * i * i) + givenY->m[1][1] * (i * i) + givenY->m[2][1] * (i) + givenY->m[3][1];

      add_edge(edges, x, y, 0, x1, y1, 0);

      x = x1;
      y = y1;
    }
  } else if (type == BEZIER){
    struct matrix * bezier_inverse = make_bezier();
    struct matrix * givenX = generate_curve_coefs(x0, x1, x2, x3, BEZIER);
    struct matrix * givenY = generate_curve_coefs(y0, y1, y2, y3, BEZIER);

    matrix_mult(bezier_inverse, givenX);
    matrix_mult(bezier_inverse, givenY);

    float i = 0;
    float x = givenX->m[0][0] * (i * i * i) + givenX->m[1][0] * (i * i) + givenX->m[2][0] * (i) + givenX->m[3][0];
    float y = givenY->m[0][1] * (i * i * i) + givenY->m[1][1] * (i * i) + givenY->m[2][1] * (i) + givenY->m[3][1];

    for (i = i + step; i <= 1; i = i + step){;
      float x1 = givenX->m[0][0] * (i * i * i) + givenX->m[1][0] * (i * i) + givenX->m[2][0] * (i) + givenX->m[3][0];
      float y1 = givenY->m[0][1] * (i * i * i) + givenY->m[1][1] * (i * i) + givenY->m[2][1] * (i) + givenY->m[3][1];

      add_edge(edges, x, y, 0, x1, y1, 0);

      x = x1;
      y = y1;
    }
  }
}


/*======== void add_point() ==========
Inputs:   struct matrix * points
         int x
         int y
         int z
Returns:
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point( struct matrix * points, double x, double y, double z) {

  if ( points->lastcol == points->cols )
    grow_matrix( points, points->lastcol + 100 );

  points->m[0][ points->lastcol ] = x;
  points->m[1][ points->lastcol ] = y;
  points->m[2][ points->lastcol ] = z;
  points->m[3][ points->lastcol ] = 1;
  points->lastcol++;
} //end add_point

/*======== void add_edge() ==========
Inputs:   struct matrix * points
          int x0, int y0, int z0, int x1, int y1, int z1
Returns:
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/
void add_edge( struct matrix * points,
	       double x0, double y0, double z0,
	       double x1, double y1, double z1) {
  add_point( points, x0, y0, z0 );
  add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
         screen s
         color c
Returns:
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/
void draw_lines( struct matrix * points, screen s, color c) {

 if ( points->lastcol < 2 ) {
   printf("Need at least 2 points to draw a line!\n");
   return;
 }

 int point;
 for (point=0; point < points->lastcol-1; point+=2)
   draw_line( points->m[0][point],
	      points->m[1][point],
	      points->m[0][point+1],
	      points->m[1][point+1],
	      s, c);
}// end draw_lines









void draw_line(int x0, int y0, int x1, int y1, screen s, color c) {

  int x, y, d, A, B;
  //swap points if going right -> left
  int xt, yt;
  if (x0 > x1) {
    xt = x0;
    yt = y0;
    x0 = x1;
    y0 = y1;
    x1 = xt;
    y1 = yt;
  }

  x = x0;
  y = y0;
  A = 2 * (y1 - y0);
  B = -2 * (x1 - x0);

  //octants 1 and 8
  if ( abs(x1 - x0) >= abs(y1 - y0) ) {

    //octant 1
    if ( A > 0 ) {

      d = A + B/2;
      while ( x < x1 ) {
	plot( s, c, x, y );
	if ( d > 0 ) {
	  y+= 1;
	  d+= B;
	}
	x++;
	d+= A;
      } //end octant 1 while
      plot( s, c, x1, y1 );
    } //end octant 1

    //octant 8
    else {
      d = A - B/2;

      while ( x < x1 ) {
	//printf("(%d, %d)\n", x, y);
	plot( s, c, x, y );
	if ( d < 0 ) {
	  y-= 1;
	  d-= B;
	}
	x++;
	d+= A;
      } //end octant 8 while
      plot( s, c, x1, y1 );
    } //end octant 8
  }//end octants 1 and 8

  //octants 2 and 7
  else {

    //octant 2
    if ( A > 0 ) {
      d = A/2 + B;

      while ( y < y1 ) {
	plot( s, c, x, y );
	if ( d < 0 ) {
	  x+= 1;
	  d+= A;
	}
	y++;
	d+= B;
      } //end octant 2 while
      plot( s, c, x1, y1 );
    } //end octant 2

    //octant 7
    else {
      d = A/2 - B;

      while ( y > y1 ) {
	plot( s, c, x, y );
	if ( d > 0 ) {
	  x+= 1;
	  d+= A;
	}
	y--;
	d-= B;
      } //end octant 7 while
      plot( s, c, x1, y1 );
    } //end octant 7
  }//end octants 2 and 7
} //end draw_line
