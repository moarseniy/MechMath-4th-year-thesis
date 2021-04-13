
#include <iostream>
#include <fstream>
#include <vector>
//#include "Linal2.h"
#include "VTKfile.h"

using namespace std;

void MakeVTKfile(){

}

//void MakeVTKfile(std::string output_vtk, MyArray nodesX, MyArray nodesY,  std::vector<Element> elements, MyArray displacements) {
//    fstream outvtk;
//    outvtk.open(output_vtk, fstream::out);
//    outvtk << "# vtk DataFile Version 1.0\nresults.vtk  2D Unstructured Grid of Linear Triangles\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS "
//           << nodesX.get_size() << " double\n";
//    for (int i = 0; i < nodesX.get_size(); i++) {
//        outvtk << nodesX[i] << " " << nodesY[i] << " " << 0.0 << "\n";
//    }

//    outvtk << "CELLS " << elements.size() << " " << elements.size() * 4 << "\n";
//    for (int i = 0; i < elements.size(); i++) {
//        outvtk << 3 << " " << elements[i].nodesIds[0] << " " << elements[i].nodesIds[1] << " " << elements[i].nodesIds[2] << "\n";
//    }

//    outvtk << "CELL_TYPES " << elements.size() << "\n";
//    for (int i = 0; i < elements.size(); i++) {
//        outvtk << 5 << "\n";
//    }

//    outvtk << "\nPOINT_DATA " << nodesX.get_size() << "\n";

//    outvtk << "VECTORS displacements double\n";
//    for (int i = 0; i < displacements.get_size() - 1; i += 2) {
//        outvtk << displacements[i] << " " << displacements[i + 1] << " 0.0\n";
//    }

//    outvtk << "\nSCALARS summary float\nLOOKUP_TABLE default\n";
//    for (int i = 0; i < displacements.get_size() - 1; i += 2) {
//        outvtk << std::sqrt(displacements[i] * displacements[i] + displacements[i + 1] * displacements[i + 1]) << "\n";
//    }
//}
