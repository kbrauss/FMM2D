/** Box.cc
 *
 *  Created on: Jun 21, 2016
 *      Author: Keith D Brauss
 */

#include <complex>
#include <vector>
#include <iostream>
#include <string>

#include "Box.h"
#include "Util.h"

using namespace std;

/**
 * The Problem Type
 *
 * We want to apply the FMM to a problem that has the following properties
 * - The domain is the unit square [0,1]x[0,1]
 * - The refinement is quadtree where one cell is a square and is split into
 *   4 squares of equal size
 * - The source points x and the target points y will be the same
 *   (the source points x will also be the target points y).
 *   The cells at each refinement level are squares
 * - The number of source (or target) particles in each cell at the lowest refinement level (l = L) are known
 *   - We will set this number of particles per cell at l = L to nParticlesPerCell = 4
 * - The particles' locations are determined by formulas that use the coordinates of the corners of their cells
 *   - Example: - Suppose that the corners are ll=(0.0,0.0), lr=(1.0,0.0), ul=(0.0,1.0), and ur=(1.0,1.0),
 *                where ll stands for lower-left, lr stands for lower-right, ul stands for upper left,
 *                ur stands for upper right
 *              - We wish the four particles to be inside the cell and 1/4 away from each corner w.r.t. the x and y axis
 *              - Therefore, assuming that we are working with squares for cells, we can measure the length
 *                of one side of the cell \f$s = \sqrt{ (ll[0] - lr[0])^2 + (ll[1] - lr[1])^2 }\f$.
 *              - We want 1/4 of this distance.  Let d = 0.25s
 *              - Lets refer to the particles as points in the source vector x, say x1, x2, x3, x4
 *              - Then let x1[0] = ll[0] + d;  x1[1] = ll[1] + d;
 *                         x2[0] = lr[0] - d;  x2[1] = lr[1] + d;
 *                         x3[0] = ul[0] + d;  x3[1] = ul[1] - d;
 *                         x4[0] = ur[0] - d;  xr[1] = ur[1] - d;
 *
 *                ul(0,1)  *------------------------* ur(1,1)
 *                         |                        |
 *                         |                        |
 *                         |      *          *      |
 *                         |      x3         x4     |
 *                         |                        |
 *                         |                        |
 *                         |                        |
 *                         |      *          *      |
 *                         |      x1         x2     |
 *                         |                        |
 *                ll(0,0)  *------------------------* lr(1,0)
 *
 *                         |------------s-----------|
 *                         |--d---|          |---d--|
 *
 *
 * The series representing the mother function is truncated at the index value p.
 * - Therefore, the coefficient arrays: c, dtilde, and d will have a size of p+1
 *   (since the series index starts counting at zero).
 * For each refinement level l, there are 4^l cells
 * - for l = 1, there are 4^1 = 4 cells
 * - for l = 2, there are 4^2 = 16 cells
 * - for l = 3, there are 4^3 = 64 cells
 * - ...
 * - for l = L, there are 4^L cells
 * Therefore, since the number of particles per cell at the lowest refinement level is set,
 * the total number of points in the domain is known and equals nParticlesPerCell*pow(4,L)
 * - Example: If the lowest refinement level is L = 3 and nParticlesPerCell = 4,
 *            then the total number of particles is 4*4^3 = 4*64 = 256 particles in the domain
 *
 */

/**
 * Header Interface for Class Point
 *
class Box
{
  public:
    int DEFAULT_P = 12;
    int DEFAULT_LEVEL=3;
    int DEFAULT_INDEX=0;

    int level;                             // refinement level of box
    int index;                             // cell index of box (or cell)
    int p;                                 // p is the index at which the series are truncated

//    unsigned int nParticlesPerCell = 4;    // points in each cell at lowest refinement level (l = L)
//    int totalParticles;                     // total particles in domain (unit square)


    bool empty;                            // is box empty


    std::vector<std::complex<double> > c;
    std::vector<std::complex<double> > dtilde;
    std::vector<std::complex<double> > d;

    //
    // The number of source points (and target points) in a box will depend
    // on the refinement level.  The lowest level l = L will have 4 particles
    // per box.  As we work up through the levels to level l = 0, the number
    // of particles will increase by a multiple of 4 when going from one level
    // up to the next (four children for each parent)

    std::vector<Point> x;  // source points
    std::vector<Point> y;  // target points

    // Creates a new instance of Node
    Box();
    Box(int level, int index, int p);

    int       getLevel() { return level; };
    void      printLevel() { std::cout << "Box level is " << level << '\n'; };
    void      setLevel(int i) { this->level = i; };
    int       getIndex() { return index; };
    void      printIndex() { std::cout << "Box index is " << index << '\n'; };
    void      setIndex(int i) { this->index = i; };
    Point     getCenter();
    double    getSize() { return std::pow(2.0, -level); };

    void      setP(int p) { this->p = p; };
    bool      isEmpty() { return empty; };

    std::vector<std::complex<double> > getC() {return c; };
    void                               addToC(std::vector<std::complex<double> > &increment);
    std::string                        cToString();
    void                               printC();

    std::vector<std::complex<double> > getD() {return d; };
    void                               addToD(const std::vector<std::complex<double> > &increment);
    std::string                        dToString();
    void                               printD();

    std::vector<std::complex<double> > getDtilde() {return dtilde; };
    void                               addToDtilde(const std::vector<std::complex<double> > &increment);
    std::string                        dtildeToString();
    void                               printDtilde();


    void               addX(Point &p);
    int                getSizeX() { return this->x.size(); };
    std::vector<Point> getX() { return this->x; };
    void               printSizeX() { std::cout << "Box sizeX is " << x.size() << "\n"; };

    void               addY(Point &p);
    int                getSizeY() { return this->y.size(); };
    std::vector<Point> getY() { return this->y; };
    void               printSizeY() { std::cout << "Box sizeY is " << y.size() << "\n"; };

    std::string        toString();
    int                getParentIndex();
    void               getNeighborsIndex        (std::vector<int> &neighbor_indexes);
    void               getParentsNeighborsIndex (std::vector<int> &parents_neighbor_indexes);

    void               getNeighborsE4Index      (std::vector<int> &neighborE4_indexes);

    void               getChildrenIndex(std::vector<int> &children_indexes);

  private:

    void               getChildrenIndexOfBox(int levelOfBox, int indexOfBox, std::vector<int> &children_indexes_of_box);

};
*/


Box::Box()
   :
   level(DEFAULT_LEVEL),
   index(DEFAULT_INDEX),
   p(DEFAULT_P),
   empty(true),
   c(p+1),
   dtilde(p+1),
   d(p+1)
{
  for (int i=0; i<p; ++i)
  {
	// intiializing coefficients
    c[i] = 0.0;
    dtilde[i] = 0.0;
    d[i] = 0.0;
  }
}


Box::Box(int level, int index, int p)
   :
   level(level),
   index(index),
   p(p),
   empty(true),
   c(p+1),
   dtilde(p+1),
   d(p+1)
{
  for (int i=0; i<p; ++i)
  {
	// intiializing coefficients
    c[i] = 0.0;
    dtilde[i] = 0.0;
    d[i] = 0.0;
  }
}

Point Box::getCenter()
{
  Util util;
  std::complex<double> ll_corner = util.uninterleave(this->index, this->level);
  std::complex<double> middle(0.5,0.5);
  ll_corner += middle;
  ll_corner *= getSize();
  Point center(ll_corner);
  return center;
}


void Box::addToC(std::vector<std::complex<double> > &increment)
{
  for (unsigned int i=0; i<c.size(); ++i)
    c[i] = c[i] + increment[i];
}

std::string Box::cToString()
{
  std::string ansc = "        C: ";
  for (unsigned int i=0; i<c.size(); ++i)
  {
    ansc+= "(" + std::to_string(c[i].real()) + " " + std::to_string(c[i].imag()) + ") ";
  }
  return ansc + "\n";
}

void Box::printC()
{
  std::cout << cToString() << "\n";
}

void Box::addToD(const std::vector<std::complex<double> > &increment)
{
  for (unsigned int i=0; i<d.size(); ++i)
    d[i] = d[i] + increment[i];
}

std::string Box::dToString()
{
  std::string ansd = "        D: ";
  for (unsigned int i=0; i<d.size(); ++i)
  {
    ansd+= "(" + std::to_string(d[i].real()) + " " + std::to_string(d[i].imag()) + ") ";
  }
  return ansd + "\n";
}

void Box::printD()
{
  std::cout << dToString() << "\n";
}


void Box::addToDtilde(const std::vector<std::complex<double> > &increment)
{
  for (unsigned int i=0; i<dtilde.size(); ++i)
    dtilde[i] = dtilde[i] + increment[i];
}

std::string Box::dtildeToString()
{
  std::string ansdt = "        Dtilde: ";
  for (unsigned int i=0; i<dtilde.size(); ++i)
  {
    ansdt+= "(" + std::to_string(dtilde[i].real()) + " " + std::to_string(dtilde[i].imag()) + ") ";
  }
  return ansdt + "\n";
}

void Box::printDtilde()
{
  std::cout << dtildeToString() << "\n";
}


void Box::addX(Point &p)
{
  this->x.push_back(p);
}

void Box::addY(Point &p)
{
  this->y.push_back(p);
}

std::string Box::toString()
{
  std::string ans = "box (l = " + std::to_string(level) + ", n = " + std::to_string(index) + ") \n";
  std::string ansc = "        C: ";
  std::string ansdt = "   Dtilde: ";
  std::string ansd = "        D: ";
  for (unsigned int i=0; i<c.size(); ++i)
  {
    ansc+= "(" + std::to_string(c[i].real()) + " " + std::to_string(c[i].imag()) + ") ";
    ansdt+="(" + std::to_string(dtilde[i].real()) + " " + std::to_string(dtilde[i].imag()) + ") ";
    ansd+= "(" + std::to_string(d[i].real()) + " " + std::to_string(d[i].imag()) + ") ";
  }
  ans += ansc + "\n" + ansdt + "\n" + ansd + "\n";
  return ans;
}

/**
 * Explanation of getParentIndex
 *
 *  *       ^ y-axis
 *       |
 *       |
 *       |
 *       |
 *       |
 *
 * 1.0   ------------------------------------------------------------------
 *       |       |       |       |       ||       |       |       |       |
 *       |  21   |   23  |   29  |   31  ||   53  |   55  |   61  |   63  |
 *       |       |       |       |       ||       |       |       |       |
 * 0.875 --------5---------------7----------------13--------------15-------
 *       |       |       |       |       ||       |       |       |       |
 *       |  20   |   22  |   28  |   30  ||   52  |   54  |   60  |   62  |
 *       |       |       |       |       ||       |       |       |       |
 * 0.75  ----------------1--------------------------------3----------------
 *       |       |       |       |       ||       |       |       |       |
 *       |  17   |   19  |   25  |   27  ||   49  |   51  |   57  |   59  |
 *       |       |       |       |       ||       |       |       |       |
 * 0.625 --------4---------------6----------------12--------------14-------
 *       |       |       |       |       ||       |       |       |       |
 *       |  16   |   18  |   24  |   26  ||   48  |   50  |   56  |   58  |
 *       |       |       |       |       ||       |       |       |       |
 * 0.5   __________________________________________________________________
 *
 *       |       |       |       |       ||       |       |       |       |
 *       |   5   |   7   |   13  |   15  ||   37  |   39  |   45  |   47  |
 *       |       |       |       |       ||       |       |       |       |
 * 0.375 --------1---------------3----------------9---------------11-------
 *       |       |       |       |       ||       |       |       |     o |
 *       |   4   |   6   |   12  |   14  ||   36  |   38  |   44  |   46  |
 *       |       |       |       |       ||       |       |       |       |
 * 0.25  ----------------0--------------------------------2----------------
 *       |       |       |       |       ||       |       |       |       |
 *       |   1   |   3   |   9   |   11  ||   33  |   35  |   41  |   43  |
 *       |       |       |       |       ||       |       |       |       |
 * 0.125 --------0---------------2----------------8---------------10-------
 *       |       |       |       |       ||       |       |       |       |
 *       |   0   |   2   |   8   |   10  ||   32  |   34  |   40  |   42  |
 *       |       | x     |       |       ||       |       |       |       |
 * 0.0   ------------------------------------------------------------------   --------------->
 *      0.0     0.125   0.25    0.375    0.5     0.625   0.75    0.875   1.0                      x-axis
 *
 *     This is level l=3 refinement
 *
 *
 * going to parent level (level - 1)
 * shifting the index of this box to right 2 bits (syntax is index >> 2)
 *
 * Example: index n = 8 and level l = 2
 *
 *          At level l = 2 have 4^l = 4^2 = 16 cells
 *          At parent level l = 2-1 = 1 have 4^l = 4^1 = 4 cells
 *
 *          In binary n = 8 is  00001000
 *          Shifting 00001000 to right 2 bits gives 00000010 = 2^1 = 2
 *          We can see above that cell n=8 at l=2 has parent cell n=2 at l=1
 *
 * Example: index n = 12 and level l = 2
 *
 *          At level l = 3 have 4^l = 4^3 = 64 cells
 *          At parent level l = 3-1 = 2 have 4^2 = 4^2 = 16 cells
 *
 *          In binary n = 12 is  00001100
 *          Shifting 00001100 to right 2 bits gives 00000011 = 2^1 + 2^0 = 2 + 1 = 3
 *          We can see above that cell n=8 at l=2 has parent cell n=2 at l=1
 *
 */
int Box::getParentIndex()
{
  return index >> 2;
}

/**
 * Explanation of getNeighborsIndex
 *
 * The nested for loop
 *  for (int i=-1; i<=1; i++)
      for (int j=-1; j<=1; j++)
        if (  (i!=0||j!=0) && x+i>=0 && x+i<Math.pow(2, level) && y+j>=0 && y+j<Math.pow(2, level))
          ans.addElement( t.getBox(level, Util.interleave(x+i, y+j, level)) );
 *
 * Looking at the nested for loop, we see that there are 3*3 = 9 cases and that the case
 * (i,j) = (0,0) is removed as a possibility by the if statement.  Therefore, there are 8 cases and
 * these match up with the 8 possible neighbors.  The case (0,0) matches up with the center cell
 * for which we are trying to find the neighbors.
 * We can see that -1, 0, and 1 for i corresponds to left, center (zero) and right one cell length (one unit)
 * in the x direction
 * and -1, 0, and 1 for j corresponds to up, center (zero) and down one unit in the y direction from the cell of interest.
 *
 * Example: Let tmp =  0.15625 + 0.03125i
 *          Then uninterleave function of class Util returns the location of the lower left corner
 *          of the cell where this point is located.  The location is given in terms of cell lengths from the
 *          lower left corner of the domain.  For this example the point (0.15625,0.03125) is located in the
 *          lower left corner of cell n = 2. The uninterleave function would return 1 + 0i or x = 1 and y = 0
 *          implying left one cell length from the lower left corner of the domain and up zero cell lengths
 *          to reach the lower left corner of cell n = 2.
 *
 *          Let level = 3
 *          For the cell n = 2 containing the point (0.15625,0.03125) and having x = 1 and y = 0
 *          The neighbors of n = 0,1,2,3,8, and 9
 *
 * i = -1
 *   j = -1
 *     ((i=)-1!=0 || (j=)-1!=0) && ((x+i=1+-1=)0>=0) && ((x+i=1+-1=)0<8(=2^3)) && (y+j=0+-1=-1>=0) && ((y+j=0+-1=)-1<8(=2^3))
 *                T                             T                    T                        F                     T         = False
 * (False) There is no neighbor diagonally off the lower left corner of the cell n = 2
 *
 * i = -1
 *   j = 0
 *     ((i=)-1!=0 || (j=)0!=0) && ((x+i=1+-1=)0>=0) && ((x+i=1-1=)0<8(=2^3)) && (y+j=0+0=0>=0) && (y+j=0+0=0<8(=2^3))
 *                T                            T                   T                      T                 T                 = True
 * (True) There is a neighbor directly to the left of the cell n = 2
 * Recall that ans is a vector of boxes
 *    ans.addElement( t.getBox(level, Util.interleave(x+i, y+j, level)) );
 *    ans.addElement( t.getBox(3, Util.interleave(1+-1, 0+0, 3)) );
 *    ans.addElement( t.getBox(3, Util.interleave(0, 0, 3)) );
 *
 * The interleave function returns the index n = 0.
 * Since Util.interleave(0,0,3) returns the box that is zero increments horizontally and zero increments
 * vertically upward from the lower left hand corner of the domain.  The index for this cell is n = 0.
 *
 * i = -1
 *   j = 1
 *     ((i=)-1!=0 || (j=)1!=0) && ((x+i=1+-1=)0>=0) && ((x+i=1+-1=)-1<8(=2^3)) && ((y+j=0+1=)1>=0) && ((y+j=0+1=)1<8(=2^3))
 *                T                            T                     T                       T                    T           = True
 * (True) There is a neighbor diagonally off the upper left corner of the cell n = 2
 * Recall that ans is a vector of boxes
 *    ans.addElement( t.getBox(level, Util.interleave(x+i, y+j, level)) );
 *    ans.addElement( t.getBox(3, Util.interleave(1+-1, 0+1, 3)) );
 *    ans.addElement( t.getBox(3, Util.interleave(0, 1, 3)) );
 *
 * The interleave function returns the index n = 1.
 * Since Util.interleave(0,1,3) returns the box that is zero increments horizontally and one increment
 * vertically upward from the lower left hand corner of the domain.  The index for this cell is n = 1.
 *
 * i = 0
 *   j = -1
 *     ((i=)0!=0 || (j=)-1!=0) && ((x+i=1+0=)1>=0) && ((x+i=1+0=)1<8(=2^3)) && ((y+j=0+-1=)-1>=0) && ((y+j=0+-1=)-1<8(=2^3))
 *               T                           T                    T                          F                     T          = False
 * (False) There is no neighbor directly below the cell n = 2
 *
 * i = 0
 *   j = 0
 *     ((i=)0!=0 || (j=)0!=0) && ((x+i=1+0=)1>=0) && ((x+i=1+0=)1<8(=2^3)) && ((y+j=0+0=)0>=0) && ((y+j=0+0=)0<8(=2^3))
 *               F                           T                   T                        T                   T               = False
 * (False) The cell located zero units left and zero units up is the cell itself (it is its own neighbor)
 * and is automatically accounted for
 *
 * i = 0
 *   j = 1
 *     ((i=)0!=0 || (j=)1!=0) && ((x+i=1+0=)1>=0) && ((x+i=1+0=)1<8(=2^3)) && ((y+j=0+1=)1>=0) && ((y+j=0+1=)1<8(=2^3))
 *               T                           T                   T                        T                   T               = True
 * (True) There is a neighbor directly above the cell n = 2
 * Recall that ans is a vector of boxes
 *    ans.addElement( t.getBox(level, Util.interleave(x+i, y+j, level)) );
 *    ans.addElement( t.getBox(3, Util.interleave(1+0, 0+1, 3)) );
 *    ans.addElement( t.getBox(3, Util.interleave(1, 1, 3)) );
 *
 * The interleave function returns the index n = 3.
 * Since Util.interleave(1,1,3) returns the box that is one horizontal cell length and one vertical
 * cell length upward from the lower left hand corner of the domain.  The index for this cell is n = 3.
 *
 * i = 1
 *   j = -1
 *     ((i=)1!=0 || (j=)-1!=0) && ((x+i=1+1=)2>=0) && ((x+i=1+1)=2<8(=2^3)) && ((y+j=0+-1=)-1>=0) && ((y+j=0+-1=)-1<8(=2^3))
 *               T                            T                   T                          F                     T          = False
 * (False) There is no neighbor diagonally off the lower right corner of the cell n = 2
 *
 * i = 1
 *   j = 0
 *     ((i=)1!=0 || (j=)0!=0) && ((x+i=1+1=)2>=0) && ((x+i=1+1=)2<8(=2^3)) && ((y+j=0+0=)0>=0) && ((y+j=0+0=)0<8(=2^3))
 *               T                           T                   T                        T                   T               = True
 * (True) There is a neighbor directly to the right of the cell n = 2
 * Recall that ans is a vector of boxes
 *    ans.addElement( t.getBox(level, Util.interleave(x+i, y+j, level)) );
 *    ans.addElement( t.getBox(3, Util.interleave(1+1, 0+0, 3)) );
 *    ans.addElement( t.getBox(3, Util.interleave(2, 0, 3)) );
 *
 * The interleave function returns the index n = 8.
 * Since Util.interleave(1,1,3) returns the box that is two horizontal cell lengths and zero vertical
 * cell lengths upward from the lower left hand corner of the domain.  The index for this cell is n = 8.
 *
 * i = 1
 *   j = 1
 *     ((i=)1!=0 || (j=)1!=0) && ((x+i=1+1=)2>=0) && ((x+i=1+1=)2<8(=2^3)) && ((y+j=0+1=)1>=0) && ((y+j=0+1)=1<8(=2^3))
 *                T                          T                   T                        T                   T               = True
 * (True) There is a neighbor diagonally off the upper right corner of the cell n = 2
 * Recall that ans is a vector of boxes
 *    ans.addElement( t.getBox(level, Util.interleave(x+i, y+j, level)) );
 *    ans.addElement( t.getBox(3, Util.interleave(1+1, 0+1, 3)) );
 *    ans.addElement( t.getBox(3, Util.interleave(2, 1, 3)) );
 *
 * The interleave function returns the index n = 9.
 * Since Util.interleave(2,1,3) returns the box that is two increments horizontally and one increment
 * vertically upward from the lower left hand corner of the domain.  The index for this cell is n = 9.
 *
 */

void Box::getNeighborsIndex(std::vector<int> &neighbor_indexes)
{
  // Finding x,y increments (cell lengths) from the lower left corner of
  // the domain to the lower left corner of this Box (cell) using the uninterleave
  // function from class Util and placing these increments in tmp in terms of
  // real term = horizontal units from lower left corner of domain and
  // imag term = vertical units from lower left corner of the domain
  neighbor_indexes.resize(0);
  Util util;

  std::complex<double> tmp = util.uninterleave(index,level);
  int x = (int)tmp.real();
  int y = (int)tmp.imag();
  for (int i=-1; i<=1; i++)
    for (int j=-1; j<=1; j++)
      if (  (i!=0||j!=0) && x+i>=0 && x+i<std::pow(2, level) && y+j>=0 && y+j<std::pow(2, level))
        neighbor_indexes.push_back(util.interleave(x+i, y+j, level));
}

 // see getNeighborsIndex above for explanation
 void Box::getParentsNeighborsIndex(std::vector<int> &parents_neighbor_indexes)
 {
   // Finding x,y increments (cell lengths) from the lower left corner of
   // the domain to the lower left corner of this Box (cell) using the uninterleave
   // function from class Util and placing these increments in tmp in terms of
   // real term = horizontal units from lower left corner of domain and
   // imag term = vertical units from lower left corner of the domain
   parents_neighbor_indexes.resize(0);
   Util util;

   int parent_index;
   parent_index = this->getParentIndex();

   std::complex<double> tmp = util.uninterleave(parent_index, level-1);
   int x = (int)tmp.real();
   int y = (int)tmp.imag();
   for (int i=-1; i<=1; i++)
     for (int j=-1; j<=1; j++)
       if (  (i!=0||j!=0) && x+i>=0 && x+i<std::pow(2, level-1) && y+j>=0 && y+j<std::pow(2, level-1))
         parents_neighbor_indexes.push_back(util.interleave(x+i, y+j, level-1));
 }

/**
 * Explanation of getNeighborsE4Index
 *
 * The member function gets the interaction list for this box
 * The interaction list or set E_4 is made up of this box's parent's
 * nearest neighbor's children (peers of this box).
 * However, the list does not include the nearest neighbors of this box.
 * We obtain the interaction list E_4 by
 * (1) getting the indices of the near neighbors of this box
 *     - we will check to make sure that these cell indexes are
 *     - not in the interaction list set E_4 that will be returned
 * (2) getting the indexes of this box's parent's nearest neighbors
 *     - we will then get the children of each of these parent boxes
 *     - the children then form the interaction list - minus the near neighbors of this box
 * (3) getting the set of indices for the children of this box's parent's nearest neighbors
 * (4) comparing and removing the indices of this box's near neighbors from the list (set)
 *     of step (3)
 */
void Box::getNeighborsE4Index(std::vector<int> &neighborE4_indexes)
{
  neighborE4_indexes.resize(0);

  // getting the indexes of the near neighbors of this box
  std::vector<int> neighbor_indexes;
  neighbor_indexes.resize(0);
  this->getNeighborsIndex(neighbor_indexes);

  // getting the (parent) indexes of the near neighbors of the parent of this box
  std::vector<int> parents_neighbor_indexes;
  parents_neighbor_indexes.resize(0);
  this->getParentsNeighborsIndex(parents_neighbor_indexes);

//for (unsigned int i=0; i<neighbor_indexes.size(); ++i)
//    std::cout << "neighbor_indexes[" << i << "] = " << neighbor_indexes[i] << "\n";
//for (unsigned int j=0; j<parents_neighbor_indexes.size(); ++j)
//    std::cout << "parents_neighbor_indexes[" << j << "] = " << parents_neighbor_indexes[j] << "\n";


  // getting the (peer) indexes of the children of the neighbors of the parent
  std::vector<int> parent_neighbor_children_indexes;
  std::vector<int> box_children_indexes;
  parent_neighbor_children_indexes.resize(0);
  for (unsigned int i=0; i<parents_neighbor_indexes.size(); ++i)
  {
	box_children_indexes.resize(0);
    getChildrenIndexOfBox(level-1, parents_neighbor_indexes[i],
    		              box_children_indexes);
    for (unsigned int m=0; m<box_children_indexes.size(); ++m)
    parent_neighbor_children_indexes.push_back(box_children_indexes[m]);
  }

  // comparing and removing the indexes from the parent_neighbor_children that
  // are also near neighbors of this box
  bool flag_not_in_list = false;  // true that parent_neighbor_children_indexes[i] is not
                                  // in the interaction list (is a near neighbor) of this box
  for (unsigned int i=0; i<parent_neighbor_children_indexes.size(); ++i)
  {
	flag_not_in_list =  false;

	for (unsigned int j=0; j<neighbor_indexes.size(); ++j)
      if (parent_neighbor_children_indexes[i] == neighbor_indexes[j])
        flag_not_in_list = true;

    if (flag_not_in_list == true)
    {
      // not adding the index 'parent_neighbor_children_indexes[i]'
      // to the interaction list neighborE4_indexes
    }
    else
    {
      neighborE4_indexes.push_back(parent_neighbor_children_indexes[i]);
    }

  }

}

/** Explanation of getChildrenIndex
 *
 *  For 2D, quadtree structure we have four children indices to get
 *
 *  In the for loop
 *    for (int i=0; i<4; ++i)
 *      children_indexes.push_back((index<<2)+i);
 *
 *  the command
 *      (index<<2)+i;
 *
 *  shifts the index to the left two bits
 *
 *  Example: l = 2 and index n = 9
 *
 *           n = 9 = 2^3 + 2^0 = 00001001
 *
 *           Shifting n = 9 to the left two bits (index << 2)
 *           we have 00100100 and this is
 *           00100100 = 2^2 + 2^5 = 4 + 32 = 36
 *
 *           We can see that this is the index of the lower
 *           left (ll) child (l=3 and n=36) of cell l=2 and n = 9
 *           (See the grid diagram above.)
 *
 *           The other three children are obtained by adding
 *           1, 2, and 3 to the index of this child
 *           (36 + 1 = 37 (ul), 36 + 2 = 38 (lr), and 36 + 3 = 39 (ur))
 *
 */
void Box::getChildrenIndex(std::vector<int> &children_indexes)
{
  children_indexes.resize(0);
  for (int i=0; i<4; ++i)
    children_indexes.push_back((this->index<<2)+i);
}

void Box::getChildrenIndexOfBox(int levelOfBox, int indexOfBox, std::vector<int> &children_indexes_of_box)
{
  // need to assert that levelOfBox is between 2 and 8
  // if refinement level is greater than 8, then bitwise operations will be incorrect (8 bits)
  // if refinement level is less than 2, then FMM does not apply (series approximation convergence not guaranteed)
  children_indexes_of_box.resize(0);
  for (int i=0; i<4; ++i)
    children_indexes_of_box.push_back((indexOfBox<<2)+i);
}
