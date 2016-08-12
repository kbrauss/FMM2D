/*
 * FmmTree.cc
 *
 *  Created on: Jun 21, 2016
 *      Author: Keith D Brauss
 */

#include <complex>
#include <vector>
#include <iostream>
#include <limits>
#include <cmath>
#include <cassert>

#include "FmmTree.h"
#include "Box.h"
#include "Point.h"


using namespace std;

/**
 * Header Interface for Class Point
 *
class FmmTree
{
  public:
    int MAX_NUM_LEVEL=8;
    int DEFAULT_NUM_LEVEL=3;

    int dimension = 2;

    int numOfLevels;                         // number of levels is l - 1
    int currLevel;                           // current level is l

    std::vector<Point> x;
    std::vector<Point> y;

    Potential potential;

    std::vector<std::vector<Box> > tree_structure;              // an array of structs

    long numOpsIndirect;
    long numOpsDirect;

    FmmTree();                                // Constructor
    FmmTree(int level, Potential potential);  // Constructor
    FmmTree(int level, std::vector<Point> source, std::vector<Point> target, Potential potential);

    void initStruct();

    int getNumOfLevels() { return this->numOfLevels; };
    int getClusterThreshold();
    int getIndex(std::vector<Point> &z, Point &p);
    Box getBox(int level, int index) { return tree_structure[level][index]; };

    void printX ();
    void printY ();
    void printBoxInformation();
    void printTreeStructure();
    std::vector<double> solve(std::vector<double> &u);
    std::vector<double> solveDirect(std::vector<double> &u);

  private:
    void upwardPass(std::vector<double> &u);
    void downwardPass1();
    void downwardPass2();

};
*/

// construct the most basic tree (number of refinement levels is 4)
// that can still use FMM
FmmTree::FmmTree()
       :
       numOfLevels(4),
       currLevel(numOfLevels-1),
       tree_structure(numOfLevels),
       numOpsIndirect(0)
{}


// Explanation of Constructor FmmTree:
//
// the first argument 'level' of the constructor passes in the largest
// refinement level L.  However, this value is only true when counting of
// levels starts at the index l = 1
// Indexing in C++ starts with zero, therefore with respect to
// starting the count at l = 0, the largest refinement level
// (the first argument 'level' of the constructor below is)
// is L = level - 1
// We see this taken into account with currLevel and the
// numOfLevels-1 in the initializer list
FmmTree::FmmTree(int level, std::vector<Point> &sources, std::vector<Point> &targets, Potential &potential)
       :
       numOfLevels(level),
       currLevel(numOfLevels-1),
       potential(potential),
       tree_structure(numOfLevels),
       numOpsIndirect(0)
{
  // need to assert that levelOfBox is between 0 and 8
  // if refinement level is greater than 8, then bitwise operations will be incorrect (8 bits)
  // a refinement level less than 0 does not make sense (not defined)
  // (series approximation convergence not guaranteed)
  // Note: level gives refinement level based on count beginning with 1
  //       levels l = 1, 2, 3, 4, ....
  assert(level>0 && "FmmTree level < 1");
  assert(level<=8 && "FmmTree level > 8");

  unsigned int length = sources.size();
  x.resize(length);
  y.resize(length);
  for (unsigned int i=0; i<length; ++i)
  {
    x[i] = sources[i];
    y[i] = targets[i];
  }
  initStruct();
}

/**
 * Explanation of initStruct()
 *
 * First For Loop:
 *
 * Incrementing on i (i is row number of the struct matrix as well as refinement level)
 *
 * Each element of struct (let's say row number i) is initialized as a vector (column corresponding to the row)
 * of Box objects with its size being the number of cells for the refinement level i (4^i)
 *
 * Inside the first for loop is a second loop incrementing on j
 * (j is the column number of the struct matrix as well as the cell index n)
 *
 * Thinking of struct as a matrix.  Each element of the (matrix) struct is also
 * a Box object - Box(i,j,potential.getP()) where
 * - i is the row number
 *   - with respect to the Box constructor this the refinement level
 * - j is the column number
 *   - with respect to the Box constructor this is the cell index
 * - p is an integer from the potential object
 *   - with respect to the Box constructor this is the truncation index for the series approximation
 *
 * Second and Third For Loop:
 *
 * The second loop increments over i through each the source particle x[i].
 *  - For the refinement level numOfLevels-1, getBoxIndex determines the cell index n for each source particle x[i]
 *  - The addX member function of Box is called to add the source point x[i] to that cell (or box)
 *
 * The third loop increments over i through each the target particle y[i].
 *  - For the refinement level numOfLevels-1, getBoxIndex determines the cell index n for each source particle y[i]
 *  - The addY member function of Box is called to add the target point y[i] to that cell (or box)
 *
 */

void FmmTree::initStruct()
{
  int cells_in_tree_structure_i;

  for (unsigned int i=0; i<tree_structure.size(); ++i)
  {
    // creating all 4^i cells (boxes) for refinement level i
	cells_in_tree_structure_i = std::pow(4,i);
    tree_structure[i].resize(cells_in_tree_structure_i);
    // setting index, level and last term p in series sum for cell
    // in tree_structure at level i
    for (unsigned int j=0; j<tree_structure[i].size(); ++j)
    {
      tree_structure[i][j].setLevel(i);              // refinement level of box (cell) is i
      tree_structure[i][j].setIndex(j);              // index of box (cell) is j
      tree_structure[i][j].setP(potential.getP());

    }
  }
  // using getBoxIndex to perform sorting of source and target particles
  // into boxes (cells) for currLevel (numOfLevel-1)
  std::cout << "x.size() = " << x.size() << "\n";
  for (unsigned int i=0; i<x.size(); ++i)
  {
    tree_structure[numOfLevels-1][x[i].getBoxIndex(numOfLevels-1)].addX(x[i]);
  }
  for (unsigned int i=0; i<y.size(); ++i)
    tree_structure[numOfLevels-1][y[i].getBoxIndex(numOfLevels-1)].addY(y[i]);

}

// getting the largest number of points (source x or target y) in a cell
// for level numOfLevels-1 (highest refinement level when counting of l
// starts with l = 0
int FmmTree::getClusterThreshold()
{
  int ans = 0;
  for (unsigned int i=0; i<tree_structure[numOfLevels-1].size(); ++i)
  {
    int xlength = tree_structure[numOfLevels-1][i].getSizeX();
    int ylength = tree_structure[numOfLevels-1][i].getSizeY();
    if (xlength>ans)
      ans = xlength;
    if (ylength>ans)
      ans = ylength;
  }
  return ans;
}

// returns index of p in vector z
// point p must be a point in vector z
// else index returned is ans = -1
int FmmTree::getIndex(std::vector<Point> &z, Point &p)
{
  int ans = -1;
  for (unsigned int i=0; i<z.size(); ++i)
    if (z[i].equals(p))
    {
      ans = i;
      i = z.size();
    }
  return ans;
}

void FmmTree::printX()
{
  for (unsigned int i=0; i<x.size(); ++i)
  {
    std::string x_coords = this->x[i].coordToString();
	std::cout << x_coords << '\n';
  }

}

void FmmTree::printY()
{
  for (unsigned int i=0; i<x.size(); ++i)
  {
    std::string y_coords = this->y[i].coordToString();
	std::cout << y_coords << '\n';
  }

}

void FmmTree::printBoxInformation()
{
  for (unsigned int i=0; i<tree_structure[numOfLevels-1].size(); ++i)
  {
    std::cout << "For the Box at tree_structure[" << numOfLevels-1 << "][" << i << "]" << "\n";
    tree_structure[numOfLevels-1][i].printLevel();
    tree_structure[numOfLevels-1][i].printIndex();
    tree_structure[numOfLevels-1][i].printSizeX();
    tree_structure[numOfLevels-1][i].printSizeY();
  }
}


//void FmmTree::printTreeStructure()
//{
//  for (int i=0; i<numOfLevels; ++i)
//  {
//    std::cout << 'tree structure level' << i
//    		  << 'has ' << tree_structure[i].size() << 'cells' << '\n';
    //for (unsigned int j=0; j<tree_structure[i].size(); ++j)
    //  std::cout << 'cell[' << i << '][' << j << '] has '
    //            << tree_structure[i][j].getSizeX() << 'points.'
    //            << '\n';
//  }
//}

std::vector<double> FmmTree::solve(std::vector<double> &u)
{
  std::vector<double> answer;
  std::cout << "Starting updward pass..." << "\n";
  upwardPass(u);
  std::cout << "Completed updward pass..." << "\n";
  std::cout << "Starting downward pass 1..." << "\n";
  downwardPass1();
  std::cout << "Completed downward pass 1..." << "\n";
  std::cout << "Starting downward pass 2..." << "\n";
  downwardPass2();
  std::cout << "Completed downward pass 2..." << "\n";

  // v is the answer to the potential calculation using FMM
  // for each target y[i] (element of y), we will have calculated the potential
  // v[i] due to all the sources x using FMM
  std::vector<double> v;
  v.resize(y.size());

  // Explanation of Nested For Loops in Code Below - could be separate member function of FmmTree
  //
  // [0] - for each box at the highest refinement level numOfLevels
  //       (index starts on zero, so numOfLevels-1)
  //   [1] - getting a reference thisBox for the box to be worked on
  //   [2] - getting the target points yPoints of this box
  //   [3] - if there are target points in this box
  //     [4] - for each target point yPoints[j]
  //       [5] - getting a reference thisY for the target point
  //       [6] - declaring the regular part of the potential calculation
  //             where the source points x[i] are far enough away from thisBox
  //             that the potential calculation can be approximated by a series
  //       [7] - retrieving series coefficients D (could be done outside this loop)
  //       [8] - retrieving powers of R-expansion for thisY
  //       [9] - for each of p terms of R-expansion series
  //        [10] - calculating and adding the first p terms of the series - a
  //               truncated approximation to the infinite series - and only
  //               taking real part of series?
  //      [11] - initializing the singular part of the potential calculation
  //             where the source points x[i] are too close to approximate
  //             the potential calculation with a series and the calculation
  //             must be done directly
  //      [12] - declaring and collecting the neighbors and this box's indices
  //             in this line and the two lines above
  //      [13] - for each near neighbor (including this box)
  //        [14-15] - obtaining a reference thisNeighborsBox to neighor's box
  //                  and x-values (sources) thisNeighborsX in neighbor's box
  //        [16] - if there are source terms thisNeighborsX in neighbor's box
  //          [17] - for each of the source terms thisNeighborsX[q]
  //            [18-19] - getting reference thisX and thisU for thisNeighborsX[q]
  //                      and the charge u for that source particle thisX
  //            [20-22] - declaring and initializing coordinates for thisX and thisY
  //            [23-25] - if thisX and thisY are not the same
  //                      using relative and absolute comparison for the cases
  //                      where thisXCoord and thisYCoord may be large or small
  //                      if thisXCoord and thisYCoord are small then maxXYOne is 1
  //                      and we are comparing absolutely and not relatively
  //                      (relative comparison - percentages - could be smaller than
  //                      machine epsilon
  //              [26-28] - calculating the potential directly (singular part)
  //                        for thisX on thisY (and taking only real part?)
  //      [29] - once completing direct calculations on thisY for all sources thisX
  //             in the near neighbors, then adding the result (singular part) to
  //             the result from the series approximations to the potential calculation
  //             for sources far enough away (regular part)
  //             Making sure to put this final result in the same location (have same index value)
  //             as the corresponding location (index value) of yPoints[j] = thisY in the vector
  //             of target points y
  //
  for (unsigned int i=0; i<tree_structure[numOfLevels-1].size(); ++i)                       // 0
  {
    Box& thisBox = tree_structure[numOfLevels-1][i];                                        // 1
    std::vector<Point> yPoints = thisBox.getY();                                            // 2
    if (yPoints.size() > 0)                                                                 // 3
    {
      for (unsigned int j=0; j<yPoints.size(); ++j)                                         // 4
      {
        Point& thisY = yPoints[j];                                                          // 5
        double regPart = 0.0;                                                               // 6
        std::vector<std::complex<double> > d = thisBox.getD();                              // 7

        std::vector<std::complex<double> > rVec                                             // 8
              = potential.getRVector(thisY.getCoord(), thisBox.getCenter().getCoord());
        numOpsIndirect+=potential.getP();
        for (unsigned int k=0; k<d.size(); ++k)                                             // 9
        {
          regPart += (d[k] * rVec[k]).real();                                               // 10
          numOpsIndirect++;
        }

        double sinPart = 0.0;                                                               // 11
        std::vector<int> neighbors_indexes_and_me;
        thisBox.getNeighborsIndex(neighbors_indexes_and_me);
        neighbors_indexes_and_me.push_back(thisBox.getIndex());                             // 12

        for (unsigned int m=0; m<neighbors_indexes_and_me.size(); ++m)                      // 13
        {
          Box& thisNeighborsBox
              = tree_structure[numOfLevels-1][neighbors_indexes_and_me[m]];                 // 14
          std::vector<Point> thisNeighborsX = thisNeighborsBox.getX();                      // 15
          if (thisNeighborsX.size() > 0)                                                    // 16
          {
            for (unsigned int q=0; q<thisNeighborsX.size(); ++q)                            // 17
            {
              Point& thisX = thisNeighborsX[q];                                             // 18
              double thisU = u[getIndex(x, thisX)];                                         // 19
              std::complex<double> inc, thisYCoord, thisXCoord;                             // 20
              thisYCoord = thisY.getCoord();                                                // 21
              thisXCoord = thisX.getCoord();                                                // 22
              double maxXY = std::max(std::abs(thisXCoord),
            	                         std::abs(thisYCoord));
              double maxXYOne = std::max(1.0,maxXY);                                        // 23
              if (std::abs(thisYCoord-thisXCoord)
                      <= std::numeric_limits<double>::epsilon()*maxXYOne)                   // 24
              {
                // Do nothing - thisYCoord and thisXCoord are the same point
            	// (up to machine epsilon)
            	// This can happen when the target and source points are the same
            	// Specifically, this happens when thisNeighborsBox is thisBox
            	// (last q-value), and thisX and thisY are a target point and
            	// a source point in the same box (and the target and source
            	// points are the same).
              }
              else // target and source points thisY and thisX are not the same             // 25
              {
                  inc = potential.direct(thisY.getCoord(), thisX.getCoord());               // 26
                  inc = inc * thisU;                                                        // 27
                  sinPart += inc.real();                                                    // 28
                  numOpsIndirect++;
              }
            }
          }
        }

        v[getIndex(y, thisY)] = sinPart + regPart;                                          // 29
      }
    }
  }

  return v;

}



// Explanation of upwardPass:
//
// It seems like the original version of this function
// should be separate from FmmTree and in another class,
// since FmmTree has to pass itself as argument in thisBox.getParent(this)
// However, a work around this pass of tree_structure
// as an argument to thisBox.getParent(this) was made by only
// getting indices from class Box member function getParent of parent.
// This avoids the need to pass tree_structure since no
// boxes are now being returned from thisBox.getParent(this).
// The ability to obtain the actual box of the parent was the
// purpose of passing tree_structure as an argument.
// The same work arounds for getting children, parent's neighbors
// and parent's neighbor's children were also done for the related
// member functions of class Box

void FmmTree::upwardPass(std::vector<double> &u)
{
  for (unsigned int i=0; i<tree_structure[numOfLevels-1].size(); ++i)
  {
    Box& thisBox = tree_structure[numOfLevels-1][i];
    std::vector<Point> xPoints = thisBox.getX();
    if (xPoints.size() > 0)
      for (unsigned int j=0; j<xPoints.size(); ++j)
      {
    	Point& thisX = xPoints[j];
        std::vector<std::complex<double> >
           B = potential.getSCoeff(thisX.getCoord(), thisBox.getCenter().getCoord());
        numOpsIndirect++;
        numOpsIndirect+=potential.getP();        // O(p) flops in getSCoeff
                                                 // p*(1 subtraction, 1 pow, 1 division, 1 mult by -1)
        double thisU = u[getIndex(x,thisX)];

        for (unsigned int k=0; k<B.size(); ++k)
        {
          B[k] = B[k]*thisU;
          std::cout << "B[" << k << "] = " << B[k] << "  ";
          numOpsIndirect++;
        }
        thisBox.addToC(B);
        numOpsIndirect+=potential.getP();
      }

  }

  // Here el stands for refinement level and
  //      k  stands for box (cell) index
  for (int el = numOfLevels-1; el>=2; --el)
  {
    std::cout << "Upward pass level " << +el << "\n";
    for (unsigned int k=0; k<tree_structure[el].size(); ++k)
    {
      Box& thisBox = tree_structure[el][k];
      int parentBoxLevel = el-1;
      Box& parentBox = tree_structure[parentBoxLevel][thisBox.getParentIndex()];
      numOpsIndirect++;
      std::complex<double> from = thisBox.getCenter().getCoord();
      numOpsIndirect++;
      std::complex<double> to = parentBox.getCenter().getCoord();
      numOpsIndirect++;

      // translating the thisBox's series that has coeffs thisBoxC
      // from its center at location 'from' = thisBox.getCenter().getCoord()
      // to its parent's center at location 'to' = parentBox.getCenter().getCoord()
      // The new series with parent center can be added to the parent's
      // C series since the powers for each term of the two series are now the same
      // see Math.cc file notes for the details
      std::vector<std::complex<double> > newCoeffs
        = potential.getSS(from, to, thisBox.getC());
      parentBox.addToC(newCoeffs);
      numOpsIndirect+=pow(potential.getP(),2);
      numOpsIndirect+=potential.getP();
    }
  }

}

// this pass is for the iteraction list E_4
// (not the neighbors)
void FmmTree::downwardPass1()
{
  std::vector<int> thisBoxNeighborsE4Indexes;
  std::complex<double> from;
  std::complex<double> to;

  for (int el=2; el<numOfLevels; ++el)
  {

	for (unsigned int k=0; k<tree_structure[el].size(); ++k)
    {
      // getting the neighbor indices and interaction list index
      thisBoxNeighborsE4Indexes.resize(0);
      Box& thisBox = tree_structure[el][k];
      thisBox.getNeighborsE4Index(thisBoxNeighborsE4Indexes);

      // translating the far field series with old coefficients C to
      // a near field series with new coefficients Dtilde (see Main.cc notes)
      for (unsigned int j=0; j<thisBoxNeighborsE4Indexes.size(); ++j)
      {
        Box& thisBoxE4Neighbor = tree_structure[el][thisBoxNeighborsE4Indexes[j]];
        from = thisBoxE4Neighbor.getCenter().getCoord();
        ++numOpsIndirect;
        to = thisBox.getCenter().getCoord();
        ++numOpsIndirect;

        std::vector<std::complex<double> > newCoeffs
                                           = potential.getSR(from, to, thisBoxE4Neighbor.getC());

        thisBox.addToDtilde(newCoeffs);
      }

    }

  }
}


void FmmTree::downwardPass2()
{
  std::complex<double> from, to;
  std::vector<int> children_indexes;

  for (unsigned int i=0; i<tree_structure[2].size(); ++i)
  {
	std::vector<std::complex<double> > DtildeCoeffs
	  = tree_structure[2][i].getDtilde();
	tree_structure[2][i].addToD(DtildeCoeffs);
    numOpsIndirect+=potential.getP();
  }

  for (int el=2; el<numOfLevels-1; ++el)
  {
    for (unsigned int m=0; m<tree_structure[el].size(); ++m)
    {
      Box& thisBox = tree_structure[el][m];
      from = thisBox.getCenter().getCoord();
      numOpsIndirect++;
      children_indexes.resize(0);
      thisBox.getChildrenIndex(children_indexes);
      for (unsigned int k=0; k<children_indexes.size(); ++k)
      {
    	Box& thisBoxChild = tree_structure[el+1][children_indexes[k]];
        to = thisBoxChild.getCenter().getCoord();
        numOpsIndirect++;
    	std::vector<std::complex<double> > newCoeffs
    	  = potential.getRR(from, to, thisBox.getD());

        thisBoxChild.addToD(newCoeffs);
        numOpsIndirect+=std::pow(potential.getP(),2);

    	std::vector<std::complex<double> > DtildeCoeffs
    	  = thisBoxChild.getDtilde();
        thisBoxChild.addToD(DtildeCoeffs);
        numOpsIndirect+=potential.getP()*2;
       }
    }
  }

}


std::vector<double> FmmTree::solveDirect(std::vector<double> &u)
{
  std::vector<double> v(y.size());
  std::complex<double> potential_direct_calculation;

  for (unsigned int j=0; j<v.size(); ++j)
  {
    for (unsigned int i=0; i<x.size(); ++i)
    {
      // taking care of relative and absolute difference
      // issues for when x[i] and y[j] are both small
      // or both large (see explanation in FmmTree member function solve
      // above)
      double maxXY = std::max(std::abs(y[j].getCoord()),
      	                         std::abs(x[i].getCoord()));
      double maxXYOne = std::max(1.0,maxXY);
      if (std::abs(y[j].getCoord()-x[i].getCoord())
                <= std::numeric_limits<double>::epsilon()*maxXYOne)
      {
        // Do nothing - y[j] and x[i] are the same point
    	// (up to machine epsilon)
    	// This happens when the target and source points are the same.
    	// Specifically, this happens when y[j] and x[i] are a target
    	// point and a source point in the same box (and the target
    	// and source points are the same).
    	// Also, physically it does not make sense for a particle y[j] to act
    	// on itself
      }
      else // target and source points y[j] and x[i] are not the same
      {
        potential_direct_calculation
          = u[i] * potential.direct(y[j].getCoord(),
      		                        x[i].getCoord());
        v[j] += potential_direct_calculation.real();
        numOpsDirect++;
      }
    }
  }

  return v;
}
