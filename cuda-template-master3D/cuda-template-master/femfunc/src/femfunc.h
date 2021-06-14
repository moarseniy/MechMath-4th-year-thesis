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

struct ElementLight
{
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

struct Couple {
    int index;
    int degree;
};

struct Coords {
    Coords(int x, int y){
        this->x = x;
        this->y = y;
    }
    int x;
    int y;
};

void FindNonZeroNumbers(int *elementId0color, int *elementId1color, int *elementId2color, int *elementId3color, int elementCount,
                        float *nodesXcuda, float *nodesYcuda, float *nodesZcuda, int nodesCount);

void CreateStiffStructure(int *ConnectMatrix, int *x, int *y_s, int size, int nodesCount);
int FindConnectivityMatrix(int *ConnectMatrix, int *rowSizes, int *elements, int *element0, int *element1, int *element2, int *element3, int nodesCount, int elementsCount);
void FindSetsColor(std::vector<Element> &elements, std::vector<Element> &elementsColor, std::vector<int> &elementsOrderColor);
void FindSetsColorSuperNew(std::vector<ElementLight> &elements, std::vector<ElementLight> &elementsColor, std::vector<int> &elementsOrderColor);
//void FindSetsColor(int *elementId0, int *elementId1, int *elementId2, int *elementId3, int *elementId0color, int *elementId1color, int *elementId2color, int *elementId3color, int elementsCount);
void FindSetsColorNew(std::vector<Element> &elements, int *allcolors, int *nodesDegrees, int nodesCount);
void FindColors(std::vector<ElementLight> elements, int *colors, int elementsCount);






#endif
