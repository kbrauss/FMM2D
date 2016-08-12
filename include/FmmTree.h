/*
 * FMMTree.h
 *
 *  Created on: Jul 13, 2016
 *      Author: dbpc
 */

#ifndef FMMTREE_H_
#define FMMTREE_H_

#include <complex>
#include <vector>

#include "Box.h"
#include "Potential.h"


class FmmTree
{
  public:
    int MAX_NUM_LEVEL=8;
    int DEFAULT_NUM_LEVEL=3;

    int dimension = 2;

    int numOfLevels;
    int currLevel;

    std::vector<Point> x;
    std::vector<Point> y;

    Potential potential;

    std::vector<std::vector<Box> > tree_structure;              // an array of structs

    long numOpsIndirect;
    long numOpsDirect;

    FmmTree();                                // Constructor
    FmmTree(int level, std::vector<Point> &source, std::vector<Point> &target, Potential &potential);

    void initStruct();

    int getClusterThreshold();
    int getNumOfLevels() { return this->numOfLevels; };
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




#endif /* FMMTREE_H_ */
