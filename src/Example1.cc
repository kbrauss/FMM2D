/*
 * Example1.cc
 *
 *  Created on: Aug 10, 2016
 *      Author: dbpc
 */

#include <vector>
#include <complex>
#include <cmath>

#include "Point.h"
#include "Example1.h"

/*
 * Header Interface for Class Example1
 *
class Example1
{
  public:


    Example1(unsigned int upper_refinement_level_L);

    std::vector<Point>  getX() { return x; };
    std::vector<Point>  getY() { return y; };
    std::vector<double> getU() { return u; };

  private:

    std::vector<Point>  x;
    std::vector<Point>  y;
    std::vector<double> u;

    // highest refinement level (refinement level count starts on l = 1)
    // since C++ indexing starts with count zero, we will work with the
    // first level being level 0 ( l = 0 ) and therefore convert this
    // level L passed to the constructor to (C++ indexing) L - 1 for
    // highest refinement level to work with C++ vectors index that
    // correspond to those levels
    unsigned int L;

    // number of cells along each side of unit square domain
    double cells_per_side;

    // dividing up length of unit square domain into intervals
    // having length of cells of lowest refinement level L
    double cell_length;
    double quarter_length;
    double three_quarter_length;

    // 4 particles per cell times pow(4,L-1) cells
    unsigned int nTotalParticles;

    void createSourcePointsX();
    void createTargetPointsY();
    void createSourceParticleChargesU();
};
*/

Example1::Example1(unsigned int upper_refinement_level_L)
        :
    	L(upper_refinement_level_L)
{
  // number of cells along each side of unit square domain
  cells_per_side = pow(2.0,L-1);

  // dividing up length of unit square domain into intervals
  // having length of cells of lowest refinement level L
  cell_length = (1.0-0.0)/cells_per_side;
  quarter_length = 0.25*cell_length;
  three_quarter_length = 0.75*cell_length;

  // 4 particles per cell times pow(4,L-1) cells
  nTotalParticles = 4 * pow(4, L - 1);

  x.resize(nTotalParticles);
  y.resize(nTotalParticles);
  u.resize(nTotalParticles);

  createSourcePointsX();
  createTargetPointsY();
  createSourceParticleChargesU();

}


/**
 * Explanation createSourcePointsX and createTargetPointsY:
 *
 * In the code below, coordinates for source (x vector) and target (y vector) particles are
 * set at the quarter and 3/4 lengths of each cell in the x and y direction
 * There are four particles per cell with coordinates
 * -  ll_corner_point = P1(x,y) = P1(ll_corner_cell + quarter_length, ll_corner + quarter_length)
 *    lr_corner_point = P2(x,y) = P2(ll_corner_cell + three_quarter_length, ll_corner + quarter_length)
 *    ul_corner_point = P3(x,y) = P3(ll_corner_cell + quarter_length, ll_corner + three_quarter_length)
 *    ur_corner_point = P4(x,y) = P4(ll_corner_cell + three_quarter_length, ll_corner + three_quarter_length)
 * where ll, lr, ul, and ur means lower-left, lower-right, upper-left, and upper-right respectively
 *
 * If the lowest refinement level is L = 3
 * Then the number of cells at L = 3 is 4^L = 4^3 = 64
 * And the number of cells along the length of the unit square domain is 2^L = 2^3 = 8 (8 * 8 = 64)
 * The length of the side of a cell at refinement level l = 3 is therefore (1.0-0.0)/8 = 0.125
 * And if each cell has cell_length = 0.125, then dividing the side length (1.0-0.0) of the unit square domain
 * by the cell length results in (1.0-0.0)/0.125 = 8

 * We work through each cell, incrementing from lower-left cell corner to lower-left cell corner
 * starting at the lower left corner (0.0,0.0) of the unit square domain.
 * The inner for loop increments along the x-axis from lower-left cell corner x-coordinate to lower-left cell
 * corner x-coordinate until reaching the last lower-left cell corner of the 8th cell along the bottom of the
 * unit square domain.  The outer for loop then increments up along the y-axis to the y-coordinate of the lower-
 * left corner for the next set of 8 cells that the inner loop works across.
 *
 * Four target and four source points are collect at each stop of the inner loop.  The points are as described
 * above: ll_corner_point, lr_corner_point, ul_corner_point, and ur_corner_point.  For our example l = 3, there
 * are 4 points per cell * 64 cells = 256 source and target points.  Therefore, the each loop does 8 iterations
 * x_coord and y_coord increment from 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, to 0.875
 *
 */


void Example1::createSourcePointsX()
{
//  double cells_per_side = pow(2.0,level);         // number of cells along each side of unit square domain
//  double cell_length = (1.0-0.0)/cells_per_side;  // dividing up length of unit square domain into intervals
	                                                // having length of cells of lowest refinement level L
//  double quarter_length = 0.25*cell_length;
//  double three_quarter_length = 0.75*cell_length;
  unsigned int n_points = 0;

  for (double y_coord=0.0; y_coord<1.0; y_coord+=cell_length)
	for (double  x_coord=0.0; x_coord<1.0; x_coord+=cell_length)
	{
	  x[n_points].setX(x_coord + quarter_length);       // source particle (point)
	  x[n_points].setY(y_coord + quarter_length);       // lower left corner of cell
	  ++n_points;
	  x[n_points].setX(x_coord + three_quarter_length); // source particle (point)
	  x[n_points].setY(y_coord + quarter_length);       // lower right corner of cell
	  ++n_points;
	  x[n_points].setX(x_coord + quarter_length);       // source particle (point)
	  x[n_points].setY(y_coord + three_quarter_length); // upper left corner of cell
	  ++n_points;
	  x[n_points].setX(x_coord + three_quarter_length); // source particle (point)
	  x[n_points].setY(y_coord + three_quarter_length); // upper right corner of cell
	  ++n_points;
	}

}


void Example1::createTargetPointsY()
{
//  double cells_per_side = pow(2.0,level);         // number of cells along each side of unit square domain
//  double cell_length = (1.0-0.0)/cells_per_side;  // dividing up length of unit square domain into intervals
	                                                // having length of cells of lowest refinement level L
//  double quarter_length = 0.25*cell_length;
//  double three_quarter_length = 0.75*cell_length;
  unsigned int n_points = 0;

  for (double y_coord=0.0; y_coord<1.0; y_coord+=cell_length)
	for (double  x_coord=0.0; x_coord<1.0; x_coord+=cell_length)
	{
	    y[n_points].setX(x_coord + quarter_length);       // target particle (point)
	    y[n_points].setY(y_coord + quarter_length);       // lower left corner of cell
	    ++n_points;
	    y[n_points].setX(x_coord + three_quarter_length); // target particle (point)
	    y[n_points].setY(y_coord + quarter_length);       // lower right corner of cell
	    ++n_points;
	    y[n_points].setX(x_coord + quarter_length);       // target particle (point)
	    y[n_points].setY(y_coord + three_quarter_length); // upper left corner cell
		++n_points;
	    y[n_points].setX(x_coord + three_quarter_length); // target particle (point)
	    y[n_points].setY(y_coord + three_quarter_length); // upper right corner of cell
	    ++n_points;
	}

}

void Example1::createSourceParticleChargesU()
{
  for (unsigned int i=0; i<u.size(); ++i)
	  u[i] = 1.0;

}



