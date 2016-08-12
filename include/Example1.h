/*
 * Example1.h
 *
 *  Created on: Aug 10, 2016
 *      Author: dbpc
 */

#ifndef EXAMPLE1_H_
#define EXAMPLE1_H_

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




#endif /* EXAMPLE1_H_ */
