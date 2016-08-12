/*
 * Point.h
 *
 *  Created on: Jul 8, 2016
 *      Author: dbpc
 */

#ifndef POINT_H_
#define POINT_H_

class Point
{
  public:
	std::complex<double> coord; // coordinates for point in plane as complex number

	/** constructors */
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




#endif /* POINT_H_ */
