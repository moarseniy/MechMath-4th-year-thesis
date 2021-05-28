#ifndef FEMFUNC_H
#define FEMFUNC_H

#include "Linal2.h"
#include <iostream>
#include <vector>

using namespace std;


struct Element
{
    void CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ);
    void FindSparseSize(std::vector<couple> &Sparse);
    Matrix B = Matrix(6, 12);
    int nodesIds[4];
};

struct Constraint {
    enum Type {
        UX = 1,
        UY = 2,
        UZ = 3,
        UXY = 4,
        UXZ = 5,
        UYZ = 6,
        UXYZ = 7
    };
    int node;
    Type type;
};

void FindNonZeroNumbers(int *elementId0color, int *elementId1color, int *elementId2color, int *elementId3color, int elementCount,
                        float *nodesXcuda, float *nodesYcuda, float *nodesZcuda, int nodesCount);

void FindSetsColor(std::vector<Element> &elements, std::vector<Element> &elementsColor, std::vector<int> &elementsOrderColor);
//void FindSetsColor(int *elementId0, int *elementId1, int *elementId2, int *elementId3, int *elementId0color, int *elementId1color, int *elementId2color, int *elementId3color, int elementsCount);








#endif
