/** Point.cc
 *  Created on: June 20, 2016
 *      Author: Keith D Brauss
 * */

#include <complex>
#include <cmath>
#include <sstream>
#include <string>
#include <limits>  // DBL_EPSILON
#include <vector>

#include "Point.h"
#include "Util.h"

/**
 * Header Interface for Class Point
 *
class Point
{
  public:
	std::complex<double> coord; // coordinates for point in plane as complex number

	// constructors
	Point();
    //Point(double x, double y);
    Point(std::complex<double> coord);
    //Point(const Point &p);

    std::complex<double> getCoord() { return coord; };
    void setCoord(std::complex<double> coord) { this->coord = coord; };
    void setY(double y_coord) { this->coord.imag(y_coord); };
    void setX(double x_coord) { this->coord.real(x_coord); };

    std::string coordToString();
    bool equals(Point &p);
    int getBoxIndex(unsigned int level);
};
 */

Point::Point()
    :
    coord(0.0,0.0)
{}

Point::Point(std::complex<double> coord)
    :
    coord(coord)
{}

//inline Point::Point(double x, double y)
//    : coord(x,y)
//{}

//inline Point::Point(std::complex<double> c)
//    : coord(c.real(),c.imag())
//{}

//inline Point::Point(const Point &p)
//    : coord(p.coord.real(),p.coord.imag())
//{}

std::string Point::coordToString()
{
  std::string result;     // string which will contain the result
  std::ostringstream convert;  // stream used for the conversion
  convert  << "Point("<< coord.real() << "," << coord.imag() << ")";
  result = convert.str(); // set 'result' to the contents of the stream
  return result;
}

bool Point::equals(Point &p)
{
  double diff = std::abs(this->coord - p.getCoord());
  if (diff  < std::numeric_limits<double>::epsilon() )
    return true;
  else
    return false;
}

/**
 * Explanation of getBoxIndex member function
 *
 * The function uses the Util class member function interleave to
 * get the index number (n) of the cell (box) that this Point point is located in
 * for the refinement level that has been specified.
 *
 * getBoxIndex consists of just one call to the interleave member function of
 * the class Util.  The member function interleave of class Util is therefore
 * explained here using the example(s) below.  Our examples are based on
 * a level l = 3 refinement (enough refinements to understand the method).
 * The indexing for level 1, 2 and 3 is shown below as well as the location
 * of the two example points x and o.  We show how getBoxIndex locates
 * the level l=3 cells where the points are located.
 *
 * Before going to the examples below, it may be helpful to look at the
 * explanation for the setbit and getbit member functions of the class Util
 * in the Util.cc file.  There you can find a review of the conversion from
 * binary to decimal.  The setting of the 8 bits that define the decimal value of an
 * integer (corresponding to the index number here) is the main tool used in
 * locating the cell where the point lies.
 *
 *       ^ y-axis
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
 *     Example 1:
 *
 *     Let
 *        x - has coordinates (0.15625,0.03125) or 0.03125 + 0.03125i
 *        x is in cell (n,l) = (2,3)
 *        where cell number n = 2 and refinement level l = 3
 *
 *     Note that in terms of bits of an integer
 *     00000010 = 2^1 = 2
 *
 *
 *     Then
 *        floor(coord.real()*std::pow(2,level)) = floor(0.15625*2^3)
 *                                              = floor(0.15625*8)
 *                                              = floor(1.25)
 *                                              = 1
 *        floor(coord.imag()*std::pow(2,level)) = floor(0.03125*2^3)
 *                                              = floor(0.03125*8)
 *                                              = floor(0.25)
 *                                              = 0
 *     Therefore applying getBoxIndex to this point results in
 *     the call the Util class's member function interleave
 *        Util::interleave(int x, int y, int level)
 *        Util.interleave(1,0,3)
 *     The interleave function passes in x=1, y=0, level=3
 *     Then calls setbit of the Util class in the for loop below
 *
 *          int ans = 0;
 *          for (unsigned int i=0; i<level=3; ++i)
 *          {
 *            ans = setbit(ans, (level-i)*2-1, getbit(x, level-i-1));
 *            ans = setbit(ans, (level-i)*2-2, getbit(y, level-i-1));
 *          }
 *          return ans;
 *
 *     On the first pass (i=0) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-1, getbit(x, level-i-1))
 *            = setbit(0, (3-0)*2-1, getbit(1, 3-0-1))
 *            = setbit(0, 3*2-1, getbit(1,2))     // 1 = 00000001
 *            = setbit(0, 5, getbit(00000*01, 2)) // getting bit 2 of int 1 (starred)
 *            = setbit(0, 5, 0)  // set int 0 bit 3 to 0 (but already 0)
 *            = 0
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-2, getbit(y, level-i-1))
 *            = setbit(0, (3-0)*2-2, getbit(0,3-0-1))
 *            = setbit(0, 4, getbit(0,2))         // 0 = 00000000
 *            = setbit(0, 4, getbit(00000*00, 2)) // getting bit zero of int 0 (starred)
 *            = setbit(0, 4, 0)  // set int 0 bit 4 to 0 (but already 0)
 *            = setbit(0, 4, 0)
 *            = 0
 *
 *     On the second pass (i=1) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-1, getbit(x, level-i-1))
 *            = setbit(0, (3-1)*2-1, getbit(1, 3-1-1))
 *            = setbit(0, 2*2-1, getbit(1,1))     // 1 = 00000001
 *            = setbit(0, 3, getbit(000000*1, 1)) // getting bit 1 of int 1 (starred)
 *            = setbit(0, 3, 0)  // set int 0 bit 3 to 0 (but already 0)
 *            = 0
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-2, getbit(y, level-i-1))
 *            = setbit(0, (3-1)*2-2, getbit(0,3-1-1))
 *            = setbit(0, 4, getbit(0,1))         // 0 = 00000000
 *            = setbit(0, 2, getbit(000000*0,1))  // getting bit 1 of int 0 (starred)
 *            = setbit(0, 2, 0)
 *            = 0
 *
 *     On the third pass (i=2) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-1, getbit(x, level-i-1))
 *            = setbit(0, (3-2)*2-1, getbit(1, 3-2-1))
 *            = setbit(0, 1*2-1, getbit(1,0))     // 1 = 00000001
 *            = setbit(0, 1, getbit(0000000*, 0)) // getting bit zero of int 1 (starred)
 *            = setbit(0, 1, 1)  // set int 0 bit 1 to 1 and get 00000010 = 2^1 = 2
 *            = 2
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-1, getbit(y, level-i-1))
 *            = setbit(2, (3-2)*2-2, getbit(0,3-2-1))
 *            = setbit(2, 0, getbit(0,0))
 *            = setbit(2, 0, getbit(00000000,0))  // getting bit zero of int 0
 *            = setbit(2, 0, 0)  // set int 2 bit 0 to 0 (but bit 0 of 00000010 already 0)
 *            = 2
 *
 *     And the position of point x at level three is in cell (n,l)
 *     with box index n = 2 and level l = 3 (see figure above for location of x)
 *
 *
 *     Example 2:
 *
 *     Let
 *        o - has coordinates (0.96875,0.34375) or 0.96875 + 0.34375i
 *        x is in cell (n,l) = (29,3)
 *        where cell number n = 29 and refinement level l = 3
 *     Note that
 *     00011011 = 2^4 + 2^3 + 2^2 + 2^0 = 16 + 8 + 4 + 1 = 29
 *
 *     Then
 *        floor(coord.real()*std::pow(2,level)) = floor(0.96875*2^3)
 *                                              = floor(0.96875*8)
 *                                              = floor(7.75)
 *                                              = 7
 *        floor(coord.imag()*std::pow(2,level)) = floor(0.34375*2^3)
 *                                              = floor(0.34375*8)
 *                                              = floor(2.75)
 *                                              = 2
 *     Looking at the illustration and the last two examples, it looks
 *     like the product and floor functions that return x and y
 *     (such as (x,y) = (7,2)) are giving the increments necessary
 *     to reach the cell from the lower left corner of the domain (0,0).
 *
 *     For example, reading the values (x,y) = (7,2) as 7 units (cell lengths)
 *     horizontally to the right from the lower left corner and 2 units up
 *     we arrive at the lower left hand corner of the cell containing the point
 *     of interest (here (0.96875,0.34375))
 *
 *     Therefore applying getBoxIndex to this point results in
 *     the call to the Util class's member function interleave
 *        Util::interleave(int x, int y, int level)
 *        Util.interleave(7,2,3)
 *     The interleave function passes in x=7, y=2, level=3
 *     Then calls setbit of the Util class in the for loop below
 *
 *          int ans = 0;
 *          for (unsigned int i=0; i<level=3; ++i)
 *          {
 *            ans = setbit(ans, (level-i)*2-1, getbit(x, level-i-1));
 *            ans = setbit(ans, (level-i)*2-2, getbit(y, level-i-1));
 *          }
 *          return ans;
 *
 *     On the first pass (i=0) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-1, getbit(x, level-i-1))
 *            = setbit(0, (3-0)*2-1, getbit(7, 3-0-1))
 *            = setbit(0, 3*2-1, getbit(7,2))     // 7 = 00000111
 *            = setbit(0, 5, getbit(00000111, 2)) // getbit 2 of int 7 (which is third 1 from right)
 *            = setbit(0, 5, 1)  // set int 0 bit 5 to 1 -> setting 00000000 to 00100000 = 2^5 = 32)
 *            = 32
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-2, getbit(y, level-i-1))
 *            = setbit(32, (3-0)*2-2, getbit(2,3-0-1))
 *            = setbit(32, 4, getbit(2,2))
 *            = setbit(32, 4, getbit(000000010, 2)) // getbit 2 of int 2 (which is second 0 from right)
 *            = setbit(32, 4, 0)  // set int 32 bit 4 to 0 -> setting 00100000 to 00100000 (bit 4 already 0)
 *            = 32
 *
 *     On the second pass (i=1) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-1, getbit(x, level-i-1))
 *            = setbit(32, (3-1)*2-1, getbit(7, 3-1-1))
 *            = setbit(32, 2*2-1, getbit(7,1))     // 7 = 00000111
 *            = setbit(32, 3, getbit(00000111, 1)) // getbit 1 of int 7 (which is second 1 from right)
 *            = setbit(32, 3, 1)  // set int 0 bit 3 to 1 (00100000 to 00101000 = 2^5 + 2^3 = 32 + 8 = 40)
 *            = 40
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-1, getbit(y, level-i-1))
 *            = setbit(40, (3-1)*2-2, getbit(2,3-1-1))
 *            = setbit(40, 2, getbit(2,1))         // 2 = 00000010
 *            = setbit(40, 2, getbit(00000010, 1)) // getbit 1 of int 2 (which is first 1 from right)
 *            = setbit(40, 2, 1)  // set int 40 bit 3 to 1 (00101000 to 00101100 = 2^5 + 2^3 + 2^2 = 32 + 8 + 4 = 44)
 *            = 44
 *
 *     On the third pass (i=2) through the for loop we have
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-1, getbit(x, level-i-1))
 *            = setbit(44, (3-2)*2-1, getbit(7, 3-2-1))
 *            = setbit(44, 1*2-1, getbit(7,0))     // 7 = 00000111
 *            = setbit(44, 1, getbit(00000111, 0)) // getting bit zero of int 7
 *            = setbit(44, 1, 1)  // set int 44 bit 1 to 1 (00101100 to 00101110 = 44 + 2^1 = 46)
 *            = 46
 *
 *        ans = setbit(int n, int pos, int setto)
 *            = setbit(ans, (level-i)*2-2, getbit(y, level-i-1))
 *            = setbit(46, (3-2)*2-2, getbit(2,3-2-1))
 *            = setbit(46, 0, getbit(2,0))
 *            = setbit(46, 0, getbit(00000010,0))  // getting bit zero of int 2
 *            = setbit(46, 0, 0)  // set int 46 bit 0 to 0 (00101110 to 00101110 = 46)
 *            = 46
 *
 *     And the position of o is 46 (see figure above)
 *
 *     Remarks:
 *
 *     From the two examples above can see that for level l = 3
 *     the for loop goes from i = 0 up to i = 2 < 3
 *
 *     Within the two statements in the for loop
 *            ans = setbit(ans, (level-i)*2-1, getbit(x, level-i-1));
 *            ans = setbit(ans, (level-i)*2-2, getbit(y, level-i-1));
 *     are the arguments
 *       (level-i)*2-1
 *       (level-i)*2-2
 *
 *     We can see from the examples that these two arguments locate the bit that is
 *     to be set to zero or one
 *     For the level l = 3
 *       For i = 0
 *         we have (level - i)*2 - 1 = (3 - 0)*2 - 1 = 6 - 1 = 5
 *         and     (level - i)*2 - 2 = (3 - 0)*2 - 2 = 6 - 2 = 4
 *       For i = 1
 *         we have (level - i)*2 - 1 = (3 - 1)*2 - 1 = 4 - 1 = 3
 *         and     (level - i)*2 - 2 = (3 - 1)*2 - 2 = 4 - 2 = 2
 *       For i = 2
 *         we have (level - i)*2 - 1 = (3 - 2)*2 - 1 = 2 - 1 = 1
 *         and     (level - i)*2 - 2 = (3 - 2)*2 - 2 = 2 - 2 = 0
 *
 *     Further, they count down the bit locations as the loop progresses (as shown above).
 *     For level l=3 they count down from bit location 5 to bit location 0
 *     In terms of bit positions of an integer, this is
 *     00*00000 (bit position 5), 000*0000 (bit position 4), 0000*000 (bit position 3)
 *     00000*00 (bit position 2), 000000*0 (bit position 1), 0000000* (bit position 0)
 *
 *     Depending on the coordinates of the point of interest, each of these bits will
 *     be set to 0 or 1.  We can see that the largest number possible for level l=3
 *     is when all the bits are set to 1 (from position 5 to position 0) resulting in
 *     the number
 *                00111111 = 0*2^7 + 0*2^6 + 1*2^5 + 1*2^4 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0
 *                         = 0     + 0     + 32    + 16    + 8     + 4     + 2     + 1
 *                         = 40 + 20 + 3 = 63
 *     And this is exactly the number of cells (64 cells counted from 0 to 63) in the refinement level l=3
 *
 *     We also see that the largest amount of cells possible occurs at level l=7 and this
 *     is
 *                11111111 = 1*2^7 + 1*2^6 + 1*2^5 + 1*2^4 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0
 *                         = 128   + 64    + 32    + 16    + 8     + 4     + 2     + 1
 *                         = 192 + 40 + 20 + 3 = 192 + 63 = 255
 *     So at level 7 (l = 7) we have 256 cells counted from 0 to 255
 *     It looks like for higher refinements we will need a different data type than int
 *     so that we can take into account more cells (maybe a long int).
 *
 *     Notes
 *
 *     The first pass through the for loop
 *      - sets the 5th bit or does not
 *        - if the bit is set to 0 (0*2^5 = 0)
 *          - the location of the point will lie in the left half
 *            (set bit to 0 implies cell number n = 0*2^5 = 0 or higher (up to 31))
 *        - if the bit is set to 1 (1*2^5 = 32)
 *          - the location of the point will be in right half of the domain
 *            (set bit to 1 implies cell number n = 1*2^5=32 or higher (up to 63))
 *      - sets the 4th bit or does not
 *        - if the bit is set to 1 ( 1*2^4 = 16) the location of the point can be in two different regions
 *          - (i) the location of the point will now be in the upper left quadrant
 *            (cell number n = 0*2^5 + 1*2^4 = 0 + 16 = 16 or higher (up to n = 16 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 31) )
 *          - (ii) the location of the point will be in upper right quadrant
 *            (cell number n = 1*2^5 + 1*2^4 = 32 + 16 = 48 or higher (up to n = 48 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 63))
 *        - if the bit is set to 0 (0*2^4 = 0) the location of the point can be in two different regions
 *          - (i) the location of the point will now be in the lower left quadrant
 *            (cell number n = 0*2^5 + 0*2^4 = 0 + 0 = 0 or higher (up to n = 0 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 15) )
 *          - (ii) the location of the point will be in lower right quadrant
 *            (cell number n = 1*2^5 + 0*2^4 = 32 + 0 = 32 or higher (up to 32 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 47))
 *
 *     The second pass through the loop
 *      - sets the 3rd bit or does not (4 cases for each situation)
 *        - if the bit is set to 1 ( 1*2^3 = 8) the location of the point can be in four different regions
 *          - (i) the location of the point will be on the right side of the lower left quadrant
 *            (cell number n =0*2^5 + 0*2^4 + 1*2^3 = 0 + 0 + 8 = 8 or higher
 *            (up to n = 8 + 1*2^2 + 1*2^1 + 1*2^0 = 15))
 *          - (ii) the location of the point will now be in the right side of the upper left quadrant
 *            (cell number n = 0*2^5 + 1*2^4 + 1*2^3 = 0 + 16 + 8 = 24 or higher
 *            (up to n = 16 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 31) )
 *          - (iii) the location of the point will be on the right side of the lower right quadrant
 *            (cell number n = 1*2^5 + 0*2^4 + 1*2^3 = 32 + 0 + 8 = 40 or higher
 *            (up to n = 40 + 1*2^2 + 1*2^1 + 1*2^0 = 47))
 *          - (iv) the location of the point will be on the right side of the upper right quadrant
 *            (cell number n = 1*2^5 + 1*2^4 + 1*2^3 = 32 + 16 + 8 = 56 or higher
 *            (up to n = 48 + 1*2^3 + 1*2^2 + 1*2^1 + 1*2^0 = 63))
 *        - if the bit is set to 0 (0*2^3 = 0)
 *          - (i) the location of the point will be on the left side of the lower left quadrant
 *            (cell number n =0*2^5 + 0*2^4 + 0*2^3 = 0 + 0 + 0 = 0 or higher
 *            (up to n = 0 + 1*2^2 + 1*2^1 + 1*2^0 = 7))
 *          - (ii) the location of the point will now be in the left side of the upper left quadrant
 *            (cell number n = 0*2^5 + 1*2^4 + 0*2^3 = 0 + 16 + 0 = 16 or higher
 *            (up to n = 16 + 1*2^2 + 1*2^1 + 1*2^0 = 23) )
 *          - (iii) the location of the point will be on the left side of the lower right quadrant
 *            (cell number n = 1*2^5 + 0*2^4 + 0*2^3 = 32 + 0 + 0 = 32 or higher
 *            (up to n = 32 + 1*2^2 + 1*2^1 + 1*2^0 = 39))
 *          - (iv) the location of the point will be on the left side of the upper right quadrant
 *            (cell number n = 1*2^5 + 1*2^4 + 0*2^3 = 32 + 16 + 0 = 48 or higher
 *            (up to n = 48 + 1*2^2 + 1*2^1 + 1*2^0 = 55))
 *      - sets the 2nd bit or does not (8 cases for each situation)
 *        - if the bit is set to 1 ( 1*2^2 = 4) the location of the point can be in 8 different regions
 *          - (i) the location of the point will be on upper left side of the lower left quadrant
 *            (cell number n =0*2^5 + 0*2^4 + 0*2^3 + 1*2^2 = 0 + 0 + 0 + 4 = 4 or higher
 *            (up to n = 4 + 1*2^1 + 1*2^0 = 7))
 *          - (ii) the location of the point will be on upper right side of the lower left quadrant
 *            (cell number n = 0*2^5 + 0*2^4 + 1*2^3 + 1*2^2 = 0 + 0 + 8 + 4 = 12 or higher
 *            (up to n = 12 + 1*2^1 + 1*2^0 = 15))
 *          - (iii) the location of the point will now be in the upper left side of the upper left quadrant
 *            (cell number n = 0*2^5 + 1*2^4 + 0*2^3 + 1*2^2 = 0 + 16 + 0 + 4 = 20 or higher
 *            (up to n = 20 + 1*2^1 + 1*2^0 = 23) )
 *          - (iv) the location of the point will now be in the upper right side of the upper left quadrant
 *            (cell number n = 0*2^5 + 1*2^4 + 1*2^3 + 1*2^2 = 0 + 16 + 8 + 4 = 28 or higher
 *            (up to n = 28 + 1*2^1 + 1*2^0 = 31) )
 *          - (v) the location of the point will be on the upper left side of the lower right quadrant
 *            (cell number n = 1*2^5 + 0*2^4 + 0*2^3 + 1*2^2 = 32 + 0 + 0 + 4 = 36 or higher
 *            (up to n = 36 + 1*2^1 + 1*2^0 = 39))
 *          - (vi) the location of the point will be on the upper right side of the lower right quadrant
 *            (cell number n = 1*2^5 + 0*2^4 + 1*2^3 + 1*2^2 = 32 + 0 + 8 + 4 = 44 or higher
 *            (up to n = 44 + 1*2^1 + 1*2^0 = 47))
 *          - (vii) the location of the point will be on the upper left side of the upper right quadrant
 *            (cell number n = 1*2^5 + 1*2^4 + 0*2^3 +1*2^2 = 32 + 16 + 0 + 4 = 52 or higher
 *            (up to n = 52 + 1*2^1 + 1*2^0 = 55))
 *          - (viii) the location of the point will be on the upper right side of the upper right quadrant
 *            (cell number n = 1*2^5 + 1*2^4 + 1*2^3 + 1*2^2 = 32 + 16 + 8 + 4 = 60 or higher
 *            (up to n = 60 + 1*2^1 + 1*2^0 = 63))
 *        - if the bit is set to 0 ( 0*2^2 = 0) the location of the point can be in 8 different regions
 *          - (i) the location of the point will be on lower left side of the lower left quadrant
 *            (cell number n =0*2^5 + 0*2^4 + 0*2^3 + 0*2^2 = 0 + 0 + 0 + 0 = 0 or higher
 *            (up to n = 0 + 1*2^1 + 1*2^0 = 3))
 *          - (ii) the location of the point will be on lower right side of the lower left quadrant
 *            (cell number n = 0*2^5 + 0*2^4 + 1*2^3 + 0*2^2 = 0 + 0 + 8 + 0 = 8 or higher
 *            (up to n = 8 + 1*2^1 + 1*2^0 = 11))
 *          - (iii) the location of the point will now be in the lower left side of the upper left quadrant
 *            (cell number n = 0*2^5 + 1*2^4 + 0*2^3 + 0*2^2 = 0 + 16 + 0 + 0 = 16 or higher
 *            (up to n = 16 + 1*2^1 + 1*2^0 = 19) )
 *          - (iv) the location of the point will now be in the lower right side of the upper left quadrant
 *            (cell number n = 0*2^5 + 1*2^4 + 1*2^3 + 0*2^2 = 0 + 16 + 8 + 0 = 24 or higher
 *            (up to n = 24 + 1*2^1 + 1*2^0 = 27) )
 *          - (v) the location of the point will be on the lower left side of the lower right quadrant
 *            (cell number n = 1*2^5 + 0*2^4 + 0*2^3 + 0*2^2 = 32 + 0 + 0 + 0 = 32 or higher
 *            (up to n = 32 + 1*2^1 + 1*2^0 = 35))
 *          - (vi) the location of the point will be on the lower right side of the lower right quadrant
 *            (cell number n = 1*2^5 + 0*2^4 + 1*2^3 + 0*2^2 = 32 + 0 + 8 + 0 = 40 or higher
 *            (up to n = 40 + 1*2^1 + 1*2^0 = 43))
 *          - (vii) the location of the point will be on the lower left side of the upper right quadrant
 *            (cell number n = 1*2^5 + 1*2^4 + 0*2^3 + 0*2^2 = 32 + 16 + 0 + 0 = 48 or higher
 *            (up to n = 48 + 1*2^1 + 1*2^0 = 51))
 *          - (viii) the location of the point will be on the lower right side of the upper right quadrant
 *            (cell number n = 1*2^5 + 1*2^4 + 1*2^3 + 0*2^2 = 32 + 16 + 8 + 0 = 56 or higher
 *            (up to n = 56 + 1*2^1 + 1*2^0 = 59))
 *
 *     The third pass through the for loop will end the loop for refinement level l = 3
 *
 *     From the first two passes through the for loop we can see that
 *      - each pass specifies a cell that the point is in.
 *      - the cell determined is with respect to a refinement level specified by the counter i of the for loop
 *        - i = 0 corresponds to the first refinement level l = 1
 *        - setbit 5
 *          - if bit 5 is set to 1, the point is in the right side of the domain
 *            setbit 4
 *            - if bit 4 is set to 1, the point is in the upper cell of the right side of the domain
 *            - if bit 4 is set to 0, the point is in the lower cell of the right side of the domain
 *        - if bit 5 is set to 0, the point is on the left side of the domain
 *          - setbit 4
 *            - if bit 4 is set to 1, the point is in the upper cell of the left side of the domain
 *            - if bit 4 is set to 0, the point is in the lower cell of the left side of the domain
 *
 *      - This takes care of refinement level l = 1
 *      - We know what cell (out of the 4 of level l = 1) the point is in with respect to refinement level l = 1
 *      - We then move on to the second pass through the for loop (i = 1)
 *        to determine in which of the four cells (w.r.t the cell determined for l = 1) the point is in
 *        for refinement level l = 2
 *        - We do this in two steps as well with setbit 3 and setbit 2
 *          - if setbit 3 is on (1) then the point is on the left side of the cell of l = 1
 *          - if setbit 3 is off (0) then the point is on the right side of the cell of l = 1
 *          - if setbit 2 is on (1) then the point is in the upper cell of the side determined by setbit 3
 *          - if setbit 2 is off (0) then the point is in the lower cell of the side determined by setbit 3
 *      - Finally we pass through the for loop the last time with i = 2 and l = 3
 *        and determine the cell that the point is in w.r.t the last refinement level l = 3
 */


int Point::getBoxIndex(unsigned int level)
{
  Util util;
  return util.interleave(std::floor(coord.real()*std::pow(2,level)),
		                 std::floor(coord.imag()*std::pow(2,level)),
		                 level);
}
