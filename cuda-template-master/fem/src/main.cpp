#include "init.h"
#include "Linal2.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
//#include <cmath>


using namespace std;

struct Element
{
    void CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY);
    Matrix B = Matrix(3, 6);
    int nodesIds[3];
};

struct Constraint {
    enum Type {
        UX = 1,
        UY = 2,
        UXY = 3
    };
    int node;
    Type type;
};


void Element::CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY) {
    MyArray x(3), y(3);
    x[0] = nodesX[nodesIds[0]]; x[1] = nodesX[nodesIds[1]]; x[2] = nodesX[nodesIds[2]];
    y[0] = nodesY[nodesIds[0]]; y[1] = nodesY[nodesIds[1]]; y[2] = nodesY[nodesIds[2]];

//    x.Show();
//    y.Show();
//    cout << endl;

    //Matrix B(3, 6);

    Matrix C(3, 3);
    C(0, 0) = C(1, 0) = C(2, 0) = 1.0;
    C(0, 1) = x[0]; C(1, 1) = x[1]; C(2, 1) = x[2];
    C(0, 2) = y[0]; C(1, 2) = y[1]; C(2, 2) = y[2];

    Matrix IC(3, 3);
    //C.Show();
    C.inverse(IC, 3, 0);

    for (int i = 0; i < 3; i++) {
        B(0, 2 * i + 0) = IC(1, i);
        B(0, 2 * i + 1) = 0.0;
        B(1, 2 * i + 0) = 0.0;
        B(1, 2 * i + 1) = IC(2, i);
        B(2, 2 * i + 0) = IC(2, i);
        B(2, 2 * i + 1) = IC(1, i);
    }
    //B.Show();

    Matrix K(6, 6);
    Matrix temp1(6, 3);
    Matrix temp2(6, 3);
    Matrix temp_B(3, 6);
    double determinant = C.det(3);

    //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

    temp_B = B;
    //temp_B.Show();

    temp_B.transpose2();
    //temp_B.Show();

    temp1 = temp_B.Product(D);
    //temp2.Show();

    K = temp1.Product(B);

    //K.Show();

    K.scale(std::abs(determinant) * 0.5);
    //K.Show();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Triplet trplt11(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 0, K(2 * i + 0, 2 * j + 0));
            Triplet trplt12(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 1, K(2 * i + 0, 2 * j + 1));
            Triplet trplt21(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 0, K(2 * i + 1, 2 * j + 0));
            Triplet trplt22(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 1, K(2 * i + 1, 2 * j + 1));

            if (trplt11.get_value() != 0) {
                triplets.push_back(trplt11);
            }
            if (trplt12.get_value() != 0) {
                triplets.push_back(trplt12);
            }
            if (trplt21.get_value() != 0) {
                triplets.push_back(trplt21);
            }
            if (trplt22.get_value() != 0) {
                triplets.push_back(trplt22);
            }
        }
    }
}


double SetConstraints(int i, int j, double v, int index) {
    if (i == index || j == index) {
        return i == j ? 1.0 : 0.0;
    } else {
        return v;
    }
}

void ApplyConstraints(SparseMatrix& K, const std::vector<Constraint>& constraints) {
    std::vector<int> indicesToConstraint;

    for (std::vector<Constraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it) {
        if (it->type & Constraint::UX) {
            indicesToConstraint.push_back(2 * it->node + 0);
        }
        if (it->type & Constraint::UY) {
            indicesToConstraint.push_back(2 * it->node + 1);
        }
    }


//    for (int i = 0; i < K.get_row(); i++) {
//        for (int j = 0; j < K.get_col(); j++) {
//            for (std::vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit) {
//                K(i, j) = SetConstraints(i, j, K(i, j), *idit);
//            }
//        }
//    }
    for (int i = 0; i < K.get_size(); i++) {
        for (std::vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit) {
            K.set_value(K.get_x(i), K.get_y(i), SetConstraints(K.get_x(i), K.get_y(i), K.get_value(i), *idit));
        }
    }
}

void CalculateStressAndDeformation(std::vector<MyArray> &Deformation,
                                   std::vector<MyArray> &Stress,
                                   std::vector<double> &epsilon_mises,
                                   std::vector<double> &sigma_mises,
                                   Matrix D,
                                   std::vector<Element> elements,
                                   MyArray displacements) {
    MyArray StressVector(3);
    MyArray DeformationVector(3);
    MyArray delta(6);

    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it) {
        delta[0] = displacements[2 * it->nodesIds[0] + 0];
        delta[1] = displacements[2 * it->nodesIds[0] + 1];
        delta[2] = displacements[2 * it->nodesIds[1] + 0];
        delta[3] = displacements[2 * it->nodesIds[1] + 1];
        delta[4] = displacements[2 * it->nodesIds[2] + 0];
        delta[5] = displacements[2 * it->nodesIds[2] + 1];

        DeformationVector = it->B.Product(delta);
        StressVector = D.Product(DeformationVector);

        double sigma = sqrt(StressVector[0] * StressVector[0] - StressVector[0]
                    * StressVector[1] + StressVector[1] * StressVector[1] + 3.0 * StressVector[2] * StressVector[2]);
        sigma_mises.push_back(sigma);

        double epsilon = sqrt(DeformationVector[0] * DeformationVector[0] - DeformationVector[0]
                * DeformationVector[1] + DeformationVector[1] * DeformationVector[1] + 3.0 * DeformationVector[2] * DeformationVector[2]);
        epsilon_mises.push_back(epsilon);

        Deformation.push_back(DeformationVector);
        Stress.push_back(StressVector);

        //DeformationVector.Show();
        //StressVector.Show();
        //cout << endl;
    }
}

void MakeVTKfile(std::string output_vtk,
                 MyArray nodesX,
                 MyArray nodesY,
                 std::vector<Element> elements,
                 MyArray displacements,
                 std::vector<MyArray> Stress,
                 std::vector<double> sigma_mises,
                 std::vector<MyArray> Deformation,
                 std::vector<double> epsilon_mises) {
    fstream outvtk;
    outvtk.open(output_vtk, fstream::out);
    outvtk << "# vtk DataFile Version 1.0\nresults.vtk  2D Unstructured Grid of Linear Triangles\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS "
           << nodesX.get_size() << " double\n";
    for (int i = 0; i < nodesX.get_size(); i++) {
        outvtk << nodesX[i] << " " << nodesY[i] << " " << 0.0 << "\n";
    }

    outvtk << "CELLS " << elements.size() << " " << elements.size() * 4 << "\n";
    for (int i = 0; i < elements.size(); i++) {
        outvtk << 3 << " " << elements[i].nodesIds[0] << " " << elements[i].nodesIds[1] << " " << elements[i].nodesIds[2] << "\n";
    }

    outvtk << "CELL_TYPES " << elements.size() << "\n";
    for (int i = 0; i < elements.size(); i++) {
        outvtk << 5 << "\n";
    }

    outvtk << "\nPOINT_DATA " << nodesX.get_size() << "\n";
    outvtk << "VECTORS displacements double\n";
    for (int i = 0; i < displacements.get_size() - 1; i += 2) {
        outvtk << displacements[i] << " " << displacements[i + 1] << " 0.0\n";
    }

    outvtk << "\nSCALARS summary double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < displacements.get_size() - 1; i += 2) {
        outvtk << std::sqrt(displacements[i] * displacements[i] + displacements[i + 1] * displacements[i + 1]) << "\n";
    }

    outvtk << "\nCELL_DATA " << elements.size() << "\n";
    outvtk << "VECTORS stress double\n";
    for (int i = 0; i < Stress.size(); i++) {
        outvtk << Stress[i][0] << " " << Stress[i][1] << " " << Stress[i][2] << "\n";
    }

    outvtk << "\nSCALARS mises_stress double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < sigma_mises.size(); i++) {
        outvtk << sigma_mises[i] << "\n";
    }

    outvtk << "\nVECTORS deformation double\n";
    for (int i = 0; i < Deformation.size(); i++) {
        outvtk << Deformation[i][0] << " " << Deformation[i][1] << " " << Deformation[i][2] << "\n";
    }

    outvtk << "\nSCALARS mises_deformation double\nLOOKUP_TABLE default\n";
    for (int i = 0; i < epsilon_mises.size(); i++) {
        outvtk << epsilon_mises[i] << "\n";
    }
}




int main(void) {


    std::string name = "rect";


    std::string directory = "D:/FiniteElementMethod/prepared_meshes/";
    std::string output_vtk = "D:/FiniteElementMethod/final_results/results.vtk";
    std::string output_results = "D:/FiniteElementMethod/final_results/output.txt";

    fstream nodes_file, elements_file, loads_file, constraints_file, outfile;
    nodes_file.open(directory + name + "/nodes.txt", fstream::in);
    elements_file.open(directory + name + "/elements.txt", fstream::in);
    loads_file.open(directory + name + "/loads.txt", fstream::in);
    constraints_file.open(directory + name + "/constraints.txt", fstream::in);
    outfile.open(output_results, fstream::out);

    //callCudaKernel();

    double poissonRatio = 0.29, youngModulus = 2068.0;
    Matrix D(3, 3);

    D(0,0) = 1.0;			D(0, 1) = poissonRatio;	D(0, 2) = 0.0;
    D(1, 0) = poissonRatio;	D(1, 1) = 1.0; 			D(1, 2) = 0.0;
    D(2, 0) = 0.0;        	D(2, 1) = 0.0;        	D(2, 2) = (1.0 - poissonRatio) / 2.0;
    D.scale(youngModulus / (1.0 - pow(poissonRatio, 2.0)));

//    D(0,0)=D(1,1)=1.0;
//    D(0, 1) = D(1, 0) = poissonRatio / (1.0 - poissonRatio);
//    D(2, 2) = (1.0 - 2.0 * poissonRatio) / (2.0 * (1.0 - poissonRatio));
//    D.scale(youngModulus * (1.0 - poissonRatio) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)));


    int nodesCount;
    nodes_file >> nodesCount;

    MyArray nodesX(nodesCount);
    MyArray nodesY(nodesCount);

    for (int i = 0; i < nodesCount; ++i) {
        nodes_file >> nodesX[i] >> nodesY[i];
    }

    //nodesX.Show();
    //cout<<"\n";
    //nodesY.Show();

    std::vector<Element>   	elements;
    std::vector<Constraint>	constraints;

    int elementCount;
    elements_file >> elementCount;
    for (int i = 0; i < elementCount; ++i) {
        Element element;
        elements_file >> element.nodesIds[0] >> element.nodesIds[1] >> element.nodesIds[2];
        elements.push_back(element);
    }

    int constraintCount;
    constraints_file >> constraintCount;

    for (int i = 0; i < constraintCount; ++i) {
        Constraint constraint;
        int type;
        constraints_file >> constraint.node >> type;
        constraint.type = static_cast<Constraint::Type>(type);
        constraints.push_back(constraint);
    }

    MyArray loads(2 * nodesCount);

    int loadsCount;
    loads_file >> loadsCount;

    for (int i = 0; i < loadsCount; ++i) {
        int node;
        float x, y;
        loads_file >> node >> x >> y;
        loads[2 * node + 0] = x;
        loads[2 * node + 1] = y;
    }


    cout << "Data read success\n";

    std::vector<Triplet> triplets;
    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
    {
        it->CalculateStiffnessMatrix(D, triplets, nodesX, nodesY);
    }
    cout << "CalculateStiffnessMatrix success\n";

//     cout <<"TRIPLETS!\n";
//     for (std::vector<Triplet>::iterator it = triplets.begin(); it != triplets.end(); ++it)
//     {
//        it->Show();
//     }
//    cout<<triplets.size()<<endl;


    Matrix K_matrix(2 * nodesCount);
    MyArray result(2 * nodesCount);
    MyArray displacements(loads.get_size());


    SparseMatrix globalK(triplets.size());

    globalK.ConvertTripletToSparse(triplets);
    globalK.resize();
    cout << "new size= " << globalK.get_size() << "\n";
    globalK.Show();
    cout << "\n\n";

    globalK.SortIt();
    globalK.Show();
    //globalK.ShowAsMatrix(2 * nodesCount);
    ApplyConstraints(globalK, constraints);
    globalK.Show();
    cout << "\n\n";
    cout << "ApplyConstraints success\n";

    int nonzero = globalK.CountNonZero();
    SparseMatrix globalK2(nonzero);
    globalK2 = globalK.DeleteZeros();
    globalK2.Show();

    globalK2.CGM_solve(loads, displacements, loads.get_size());

//    callCGM_GPU(globalK2.get_x(),
//                globalK2.get_y(),
//                globalK2.get_data(),
//                loads.get_data(),
//                displacements.get_data(),
//                loads.get_size(),
//                nonzero);

    //globalK2.Show();
    //globalK2.ShowAsMatrix(2 * nodesCount);
    //displacements.Show();

    //globalK2.ConvertToMatrix(K_matrix);
    //K_matrix.LU_solve(K_matrix, loads, displacements, loads.get_size());
    //K_matrix.Show();

    //ApplyConstraints(K_matrix, constraints);
    //K_matrix.Show();



    cout << "LU_solve success\n";

    //displacements.WriteToFile();
    //outfile << displacements;


    ///////CHECK ANSWER
//    MyArray check(loads.get_size());
//    check = K_matrix.Product(displacements);
//    check.Show();
//    loads.Show();
    ///////


    for (int i = 0; i < displacements.get_size(); i++) {
        outfile << displacements[i] << " ";
        if (i % 2 != 0 && i != 0) {
            outfile << "\n";
        }
    }
    outfile << endl;

//    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it) {
//        Matrix delta(6, 1);
//        //std::cout<<2*it->nodesIds[0]<<" "<<2*it->nodesIds[1]<<" "<<2*it->nodesIds[2]<<endl;
//        delta(0, 0) = displacements[2 * it->nodesIds[0] + 0];
//        delta(1, 0) = displacements[2 * it->nodesIds[0] + 1];
//        delta(2, 0) = displacements[2 * it->nodesIds[1] + 0];
//        delta(3, 0) = displacements[2 * it->nodesIds[1] + 1];
//        delta(4, 0) = displacements[2 * it->nodesIds[2] + 0];
//        delta(5, 0) = displacements[2 * it->nodesIds[2] + 1];
//        //D * it->B * delta;
//        Matrix sigma(3, 1);
//        sigma = D.Product(it->B);
//        sigma = sigma.Product(delta);
//        double sigma_mises = sqrt(sigma(0, 0) * sigma(0, 0) - sigma(0, 0)
//            * sigma(1, 0) + sigma(1, 0) * sigma(1, 0) + 3.0 * sigma(2, 0) * sigma(2, 0));

//        outfile << sigma_mises << std::endl;
//    }


    cout << "SUCCESS\n";

    std::vector<MyArray> Deformation;
    std::vector<MyArray> Stress;
    std::vector<double> sigma_mises;
    std::vector<double> epsilon_mises;
    CalculateStressAndDeformation(Deformation, Stress, epsilon_mises, sigma_mises, D, elements, displacements);

    MakeVTKfile(output_vtk, nodesX, nodesY, elements, displacements, Stress, sigma_mises, Deformation, epsilon_mises);




    return 0;
}
