/**  Util.cc
 *
 *   Created on: Jun 20, 2016
 *       Author: Keith D Brauss
 */

#include <complex>
#include <cmath>

#include "Util.h"

/**
 * Header Interface for Class Util
 *
 *
class Util
{
  public:
    Util() {};
    int interleave(int x, int y, int level);
    std::complex<double> uninterleave(int n, int L);
    int setbit(int n, int pos, int setto);
    int getbit(int n, int pos);

};
*
*/

/** interleave function operation is explained in Point.cc */
// Interleave is used by class Point to return the cell
// (for a given refinement level) where the point is located
// In the getBoxIndex function of class Point, int x and int y
// are integers that are the floor of the product of the
// x and y coordinates of the point and the base 2^{level}
// These values appear to indicate the increments (cell lengths)
// necessary to reach the lower left hand corner of the cell
// of interest (containing the point with desired coordinates)
// Specifically, int x and int y appear to correspond to the
// horizontal increments (x) to the right of the lower left hand corner
// of the domain and vertical increments (y) upward to reach the
// lower left hand corner of the cell of interest
int Util::interleave(int x, int y, int level)
{
  if (x == 0 && y == 0)
    return 0;
  else
  {
    int ans = 0;
    for (int i=0; i<level; ++i)
    {
      ans = setbit(ans, (level-i)*2-1, getbit(x, level-i-1));
      ans = setbit(ans, (level-i)*2-2, getbit(y, level-i-1));
    }
    return ans;
  }
}

/**
 * Explanation of Uninterleave
 *
 * This function returns the location of the lower left corner of the cell whose index n
 * and level are given L.  The location is returned with respect to the lower left corner
 * of the domain and in units of cell length.
 * For example, at l = 3 and cell index n = 2,
 * the complex number 1 + 0i is returned.  This refers to one horizontal cell length from
 * the lower left corner of the domain and zero vertical cells lengths upward.
 * For example, at l = 3 and cell index n = 46,
 * the complex number returned is 7 + 2i.  This refers to the 7 horizontal cell lengths
 * and 2 vertical cell lengths to reach the lower left corner of cell n = 46 when starting
 * at the lower left corner of the domain at (0,0).
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
 *
 * Example:
 *
 * For the call uninterleave(index n, level L) = uninterleave(3, 3)
 * the refinement level  L = 3 and n = 3 are passed in
 * (see the location of cell three (index n = 3) for level L=3) in the figure above.
 *
 * xInt and yInt are initialized to 0 at the beginning of the code
 * The complex number returned by the function is xInt + i yInt
 * We look at the for loop for this example to understand what is being
 * returned in terms of the refinement and its geometry
 *
 * In the for loop below, setbit turns on/off (sets to 1 or 0)
 * the L-i-1 bit in the integer
 * - xInt if getbit(n, (L-i-1)*2+1) is 1/0 (on or off)
 * - yInt if getbit(n, (L-i-1)*2) is 1/0 (on or off)
 *
 * The largest bit positions that can be set is L-i-1 = 3 - 0 - 1 = 2
 * The largest that xInt and yInt can be after passing through the
 * for loop occurs if all the bits in positions 2, 1, and 0 are turned on (set to 1)
 * This results in the value 2^2 + 2^1 + 2^0 = 4 + 2 + 1 = 7
 * This returned complex number would then be (xInt,yInt) = (7,7)
 * This is the position of the lower left corner (ll_corner) of the cell in the upper-right
 * corner of the level L = 3 refinement.
 * From this we see that this method is able to locate all ll_corners of the cells
 * in the refinement level
 * To obtain a particular point in the cell we can multiply the values of xInt and yInt
 * by the length of the cells in the refinement.
 * For the example with L = 3 and (xInt,yInt) = (7,7)
 * - The cell lengths are 0.125 (see figure above).
 * - The lower-left corner of the cell in the upper-right corner of the refinement
 *   has position (x,y) = (0.875,0.875)
 * - This can be obtained by multiplying (xInt,yInt) = (7,7) by 0.125
 *   i.e., (7*0.125,7*0.125) = (0.875,0.875)
 *
 * Passing through the for loop with index i we have
 * i = 0:
 *   xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*2+1));
 *   yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*2));
 *   becomes
 *   xInt = setbit(0, 3-0-1, getbit(3, (3-0-1)*2+1));
 *   yInt = setbit(0, 3-0-1, getbit(3, (3-0-1)*2));
 *   implies
 *   xInt = setbit(0, 2, getbit(3, (2)*2+1)); //
 *   yInt = setbit(0, 2, getbit(3, (2)*2));   // 3 = 2^1 + 2^0 = 00000011
 *   implies
 *   xInt = setbit(0, 2, getbit(3, 5));     // getbit(3,5) = 0
 *   yInt = setbit(0, 2, getbit(3, 4));     // getbit(3,4) = 0
 *   implies
 *   xInt = setbit(0, 2, 0);                // setbit(0,2,0) = 00000000 = 0
 *   yInt = setbit(0, 2, 0);                // setbit(0,2,0) = 00000000 = 0
 *
 * i = 1:
 *   xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*2+1));
 *   yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*2));
 *   becomes
 *   xInt = setbit(0, 3-1-1, getbit(3, (3-1-1)*2+1));
 *   yInt = setbit(0, 3-1-1, getbit(3, (3-1-1)*2));
 *   implies
 *   xInt = setbit(0, 1, getbit(3, (1)*2+1));
 *   yInt = setbit(0, 1, getbit(3, (1)*2)); // 3 = 2^1 + 2^0 = 00000011
 *   implies
 *   xInt = setbit(0, 1, getbit(3, 3));     // getbit(3,3) = 0
 *   yInt = setbit(0, 1, getbit(3, 2));     // getbit(3,2) = 0
 *   implies
 *   xInt = setbit(0, 1, 0);                // setbit(0,2,0) = 00000000 = 0
 *   yInt = setbit(0, 1, 0);                // setbit(0,2,0) = 00000000 = 0
 *
 * i = 2:
 *   xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*2+1));
 *   yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*2));
 *   becomes
 *   xInt = setbit(0, 3-2-1, getbit(3, (3-2-1)*2+1));
 *   yInt = setbit(0, 3-2-1, getbit(3, (3-2-1)*2));
 *   implies
 *   xInt = setbit(0, 0, getbit(3, (0)*2+1));
 *   yInt = setbit(0, 0, getbit(3, (0)*2)); // 3 = 2^1 + 2^0 = 00000011
 *   implies
 *   xInt = setbit(0, 0, getbit(3, 1));     // getbit(3,1) = 1
 *   yInt = setbit(0, 0, getbit(3, 0));     // getbit(3,0) = 1
 *   implies
 *   xInt = setbit(0, 0, 1);                // setbit(0,0,1) = 00000001 = 1
 *   yInt = setbit(0, 0, 1);                // setbit(0,0,1) = 00000001 = 1
 *
 *   the for loop ends
 *   and therefore
 *   xInt = 1
 *   yInt = 1
 *   Therefore, the complex number xInt + i yInt = 1 + 1i is returned
 *   by the function
 */

std::complex<double> Util::uninterleave(int n, int L)
{
  std::complex<double> ll_box_corner;
  if (n==0)
  {
    ll_box_corner.real(0.0); ll_box_corner.imag(0.0);
  }
  else
  {
    int xInt = 0;
    int yInt = 0;
    for (int i=0; i<L; ++i)
    {
      xInt = setbit(xInt, L-i-1, getbit(n, (L-i-1)*2+1));
      yInt = setbit(yInt, L-i-1, getbit(n, (L-i-1)*2));
    }
    ll_box_corner.real(xInt); ll_box_corner.imag(yInt);
  }
  return ll_box_corner;
}




/** Explanation of setbit and getbit member functions
 *
 *  Recall the formula for going from binary to decimal
 *  Ex: byte 10010110   < - - - - - - - bit zero  : 0 x 2^0
 *           |   |  |                   bit one   : 1 x 2^1
 *           |   |  bit 0               bit two   : 1 x 2^2
 *           |   bit 3                  bit three : 0 x 2^3
 *           bit 7                      bit four  : 1 x 2^4
 *                                      bit five  : 0 x 2^5
 *                                      bit six   : 0 x 2^6
 *                                  +   bit seven : 1 x 2^7
 *                                -----------------------------
 *                                      2^1 + 2^2 + 2^4 + 2^7 = 2 + 4 + 16 + 128
 *                                                            = 22 + 128
 *                                                            = 150
 *
 *  1u is an unsigned value (int) with the single bit 0 set (00000001 - 8 bits w/ bit 0 set = 1)
 *  here we are using 1 instead of 1u and ints instead of unsigned ints
 *  now 1u << 0 means shift 1=00000001 to the left zero units
 *      1u << 1 means shift 1=00000001 to the left one unit = 00000010 = 1 x 2^1 = 2
 *      1u << 2 means shift 1=00000001 to the left two units = 00000100 = 1 x 2^2 = 4
 *
 *  Therefore, for the setbit function below
 *  In the if statement
 *    n | (1u<<pos) means "n union (00000001 shifted to the left pos units)"
 *    and results in a union of n with 1u<<pos (a union of all the bits that are on (ones))
 *    ex:         10100001 | 1u<<3 = 10100001 | 00001000 = 10101000
 *    in decimal:      161 | 8     =                     = 128 + 32 + 8 = 168
 *  And the code below in the else statement
 *    n & ~(1u<<pos) = n and (00000001 shifted to the left pos units)
 *    results in an intersection of n with 1u<<pos (an intersection of all bits that are on (ones))
 *    ex:         10101001 & ~(1u<<3) = 10101001 & ~(00001000) = 10101001 & 11110111 = 10100001
 */

/** turns on (set to 1) bit pos of int n if setto == 1 else turns off (set to 0) bit pos of int n */
int Util::setbit(int n, int pos, int setto)
{
  if (setto == 1)
    return n | (1<<pos);  // turning on bit 'pos' in number n (setting bit pos equal to 1)
  else
    return n & ~(1<<pos); // turning off bit 'pos' in number n (setting bit pos equal to 0)
}

/** checks whether a bit at position 'pos' of an int 'n' is on or off and returns 1 if on and 0 if off */
int Util::getbit(int n, int pos)
{
  if ((n & (1 << pos)) != 0)  // checking whether bit pos of int n is on or off (1 or 0)
    return 1;                 // returning 1 if bit pos of int n is 1 (on)
  else
    return 0;                 // returning 0 if bit pos of int n is 0 (off)
}
