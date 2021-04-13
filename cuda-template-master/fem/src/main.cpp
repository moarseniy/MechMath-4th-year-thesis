#include "init.h"
#include "Linal2.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>
//#include <cmath>


using namespace std;

struct Element
{
    void CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ);
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


void Element::CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY, MyArray& nodesZ) {
    MyArray x(4), y(4), z(4);
    x[0] = nodesX[nodesIds[0]]; x[1] = nodesX[nodesIds[1]]; x[2] = nodesX[nodesIds[2]], x[3] = nodesX[nodesIds[3]];
    y[0] = nodesY[nodesIds[0]]; y[1] = nodesY[nodesIds[1]]; y[2] = nodesY[nodesIds[2]], y[3] = nodesY[nodesIds[3]];
    z[0] = nodesZ[nodesIds[0]]; z[1] = nodesZ[nodesIds[1]]; z[2] = nodesZ[nodesIds[2]], z[3] = nodesZ[nodesIds[3]];
    //x.Show();
    //y.Show();

    Matrix C(4, 4);
    C(0, 0) = C(1, 0) = C(2, 0) = C(3, 0) = 1.0;
    C(0, 1) = x[0]; C(1, 1) = x[1]; C(2, 1) = x[2]; C(3, 1) = x[3];
    C(0, 2) = y[0]; C(1, 2) = y[1]; C(2, 2) = y[2]; C(3, 2) = y[3];
    C(0, 3) = z[0]; C(1, 3) = z[1]; C(2, 3) = z[2]; C(3, 3) = z[3];


    Matrix IC(4, 4);
    //C.Show();
    C.inverse(IC, 4, 0);
    //IC.Show();

    //Matrix B(6, 12);

    for (int i = 0; i < 4; i++) {
        B(0, 3 * i + 0) = IC(1, i);
        B(0, 3 * i + 1) = 0.0;
        B(0, 3 * i + 2) = 0.0;

        B(1, 3 * i + 0) = 0.0;
        B(1, 3 * i + 1) = IC(2, i);
        B(1, 3 * i + 2) = 0.0;

        B(2, 3 * i + 0) = 0.0;
        B(2, 3 * i + 1) = 0.0;
        B(2, 3 * i + 2) = IC(3, i);

        B(3, 3 * i + 0) = IC(2, i);
        B(3, 3 * i + 1) = IC(1, i);
        B(3, 3 * i + 2) = 0.0;

        B(4, 3 * i + 0) = 0.0;
        B(4, 3 * i + 1) = IC(3, i);
        B(4, 3 * i + 2) = IC(2, i);

        B(5, 3 * i + 0) = IC(3, i);
        B(5, 3 * i + 1) = 0.0;
        B(5, 3 * i + 2) = IC(1, i);
    }

    //B.Show();

    Matrix K(12, 12);
    Matrix temp1(12, 6);
    Matrix temp_B(6, 12);
    double determinant = C.det(4);
    //cout << determinant << endl;

    //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

    temp_B = B;
    //temp_B.Show();

    temp_B.transpose2();
    //temp_B.Show();

    temp1 = temp_B.Product(D);
    //temp2.Show();

    K = temp1.Product(B);

    //K.Show();

    K.scale(std::abs(determinant) / 6.0);
    //K.Show();

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Triplet trplt11(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 0, K(3 * i + 0, 3 * j + 0));
            Triplet trplt12(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 1, K(3 * i + 0, 3 * j + 1));
            Triplet trplt13(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 2, K(3 * i + 0, 3 * j + 2));

            Triplet trplt21(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 0, K(3 * i + 1, 3 * j + 0));
            Triplet trplt22(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 1, K(3 * i + 1, 3 * j + 1));
            Triplet trplt23(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 2, K(3 * i + 1, 3 * j + 2));

            Triplet trplt31(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 0, K(3 * i + 2, 3 * j + 0));
            Triplet trplt32(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 1, K(3 * i + 2, 3 * j + 1));
            Triplet trplt33(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 2, K(3 * i + 2, 3 * j + 2));

            if (trplt11.get_value() != 0.0) {
                triplets.push_back(trplt11);
            }
            if (trplt12.get_value() != 0.0) {
                triplets.push_back(trplt12);
            }
            if (trplt13.get_value() != 0.0) {
                triplets.push_back(trplt13);
            }
            if (trplt21.get_value() != 0.0) {
                triplets.push_back(trplt21);
            }
            if (trplt22.get_value() != 0.0) {
                triplets.push_back(trplt22);
            }
            if (trplt23.get_value() != 0.0) {
                triplets.push_back(trplt23);
            }
            if (trplt31.get_value() != 0.0) {
                triplets.push_back(trplt31);
            }
            if (trplt32.get_value() != 0.0) {
                triplets.push_back(trplt32);
            }
            if (trplt33.get_value() != 0.0) {
                triplets.push_back(trplt33);
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

void ApplyConstraints(SparseMatrix& K, const std::vector<Constraint>& constraints)
{
    std::vector<int> indicesToConstraint;

    for (std::vector<Constraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it) {
        if (it->type & Constraint::UX) {
            indicesToConstraint.push_back(3 * it->node + 0);
        }
        if (it->type & Constraint::UY) {
            indicesToConstraint.push_back(3 * it->node + 1);
        }
        if (it->type & Constraint::UZ) {
            indicesToConstraint.push_back(3 * it->node + 2);
        }
//        if (it->type & Constraint::UXY) {
//            indicesToConstraint.push_back(3 * it->node + 0);
//            indicesToConstraint.push_back(3 * it->node + 1);
//        }
//        if (it->type & Constraint::UXZ) {
//            indicesToConstraint.push_back(3 * it->node + 0);
//            indicesToConstraint.push_back(3 * it->node + 2);
//        }
//        if (it->type & Constraint::UYZ) {
//            indicesToConstraint.push_back(3 * it->node + 1);
//            indicesToConstraint.push_back(3 * it->node + 2);
//        }
//        if (it->type & Constraint::UXYZ) {
//            indicesToConstraint.push_back(3 * it->node + 0);
//            indicesToConstraint.push_back(3 * it->node + 1);
//            indicesToConstraint.push_back(3 * it->node + 2);
//        }
    }


//    for (int i = 0; i < indicesToConstraint.size(); i++) {
//        cout << indicesToConstraint[i] << " ";
//    }


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
    MyArray StressVector(6);
    MyArray DeformationVector(6);
    MyArray delta(12);

    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it) {
        delta[0] = displacements[3 * it->nodesIds[0] + 0];
        delta[1] = displacements[3 * it->nodesIds[0] + 1];
        delta[2] = displacements[3 * it->nodesIds[0] + 2];
        delta[3] = displacements[3 * it->nodesIds[1] + 0];
        delta[4] = displacements[3 * it->nodesIds[1] + 1];
        delta[5] = displacements[3 * it->nodesIds[1] + 2];
        delta[6] = displacements[3 * it->nodesIds[2] + 0];
        delta[7] = displacements[3 * it->nodesIds[2] + 1];
        delta[8] = displacements[3 * it->nodesIds[2] + 2];
        delta[9] = displacements[3 * it->nodesIds[3] + 0];
        delta[10] = displacements[3 * it->nodesIds[3] + 1];
        delta[11] = displacements[3 * it->nodesIds[3] + 2];

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
                 MyArray nodesZ,
                 std::vector<Element> elements,
                 MyArray displacements,
                 std::vector<MyArray> Stress,
                 std::vector<double> sigma_mises,
                 std::vector<MyArray> Deformation,
                 std::vector<double> epsilon_mises) {
    fstream outvtk;
    outvtk.open(output_vtk, fstream::out);
    outvtk << "# vtk DataFile Version 1.0\nresults.vtk  3D Unstructured Grid of Triangles\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS "
           << nodesX.get_size() << " double\n";
    for (int i = 0; i < nodesX.get_size(); i++) {
        outvtk << nodesX[i] << " " << nodesY[i] << " " << nodesZ[i] << "\n";
    }

    outvtk << "CELLS " << elements.size() << " " << elements.size() * 5 << "\n";
    for (int i = 0; i < elements.size(); i++) {
        outvtk << 4 << " " << elements[i].nodesIds[0] << " " << elements[i].nodesIds[1] << " " << elements[i].nodesIds[2] << " " << elements[i].nodesIds[3] << "\n";
    }

    outvtk << "CELL_TYPES " << elements.size() << "\n";
    for (int i = 0; i < elements.size(); i++) {
        outvtk << 10 << "\n";
    }

    outvtk << "\nPOINT_DATA " << nodesX.get_size() << "\n";

    outvtk << "\nVECTORS displacements double\n";
    for (int i = 0; i < displacements.get_size() - 2; i += 3) {
        outvtk << displacements[i] << " " << displacements[i + 1] << " " << displacements[i + 2] << "\n";
    }

    outvtk << "\nSCALARS summary float\nLOOKUP_TABLE default\n";
    for (int i = 0; i < displacements.get_size() - 2; i += 3) {
        outvtk << std::sqrt(displacements[i] * displacements[i] + displacements[i + 1] * displacements[i + 1] + displacements[i + 2] * displacements[i + 2]) << "\n";
    }

//    outvtk << "\nCELL_DATA " << elements.size() << "\n";
//    outvtk << "VECTORS stress double\n";
//    for (int i = 0; i < Stress.size(); i++) {
//        outvtk << Stress[i][0] << " " << Stress[i][1] << " " << Stress[i][2] << " " << Stress[i][3] << " " << Stress[i][4] << " " << Stress[i][5] << "\n";
//    }

//    outvtk << "\nSCALARS mises_stress double\nLOOKUP_TABLE default\n";
//    for (int i = 0; i < sigma_mises.size(); i++) {
//        outvtk << sigma_mises[i] << "\n";
//    }

//    outvtk << "\nVECTORS deformation double\n";
//    for (int i = 0; i < Deformation.size(); i++) {
//        outvtk << Deformation[i][0] << " " << Deformation[i][1] << " " << Deformation[i][2] << " " << Deformation[i][3] << " " << Deformation[i][4] << " " << Deformation[i][5] << "\n";
//    }

//    outvtk << "\nSCALARS mises_deformation double\nLOOKUP_TABLE default\n";
//    for (int i = 0; i < epsilon_mises.size(); i++) {
//        outvtk << epsilon_mises[i] << "\n";
//    }
}




int main(void) {


    std::string name = "bulk_task";


    std::string directory = "D:/FiniteElementMethod3D/prepared_meshes/";
    std::string output_vtk = "D:/FiniteElementMethod3D/final_results/results.vtk";
    std::string output_results = "D:/FiniteElementMethod3D/final_results/output.txt";

    fstream nodes_file, elements_file, loads_file, constraints_file, outfile;
    nodes_file.open(directory + name + "/nodes.txt", fstream::in);
    elements_file.open(directory + name + "/elements.txt", fstream::in);
    loads_file.open(directory + name + "/loads.txt", fstream::in);
    constraints_file.open(directory + name + "/constraints.txt", fstream::in);
    outfile.open(output_results, fstream::out);

    clock_t t1, t2;
    int start_time = clock();
    double time1 = 0.0;
    t1 = clock();

    double poissonRatio = 0.3, youngModulus = 2000.0;//2e+11;
    //infile >> poissonRatio >> youngModulus;
    Matrix D(6, 6);

    D(0, 0) = D(1, 1) = D(2, 2) = 1.0;
    D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(2, 1) = D(1, 2) = poissonRatio / (1.0 - poissonRatio);
    D(3, 3) = D(4, 4) = D(5, 5) = (1.0 - 2.0 * poissonRatio) / (2.0 * (1.0 - poissonRatio));

    D.scale((youngModulus * (1.0 - poissonRatio)) / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)));

    int nodesCount;
    nodes_file >> nodesCount;

    MyArray nodesX(nodesCount);
    MyArray nodesY(nodesCount);
    MyArray nodesZ(nodesCount);

    for (int i = 0; i < nodesCount; ++i) {
        nodes_file >> nodesX[i] >> nodesY[i] >> nodesZ[i];
    }

//    nodesX.Show();
//    cout<<"\n";
//    nodesY.Show();
//    cout<<"\n";
//    nodesZ.Show();

    std::vector<Element>   	elements;
    std::vector<Constraint>	constraints;

    int elementCount;
    elements_file >> elementCount;
    for (int i = 0; i < elementCount; ++i) {
        Element element;
        elements_file >> element.nodesIds[0] >> element.nodesIds[1] >> element.nodesIds[2] >> element.nodesIds[3];
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

    MyArray loads(3 * nodesCount);

    int loadsCount;
    loads_file >> loadsCount;

    for (int i = 0; i < loadsCount; ++i) {
        int node;
        float x, y, z;
        loads_file >> node >> x >> y >> z;
        loads[3 * node + 0] = x;
        loads[3 * node + 1] = y;
        loads[3 * node + 2] = z;
    }
    //loads.Show();

    cout << "Data read success\n";

    std::vector<Triplet> triplets;
    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
    {
        it->CalculateStiffnessMatrix(D, triplets, nodesX, nodesY, nodesZ);
    }
    cout << "CalculateStiffnessMatrix success\n";

//     cout <<"TRIPLETS!\n";
//     for (std::vector<Triplet>::iterator it = triplets.begin(); it != triplets.end(); ++it)
//     {
//        it->Show();
//     }
    cout<<triplets.size()<<endl;
    //int sparse_size = triplets.size();


    //Matrix K_matrix(3 * nodesCount);
    MyArray result(3 * nodesCount);
    MyArray displacements(loads.get_size());
    //int nonzero;
    //nonzero = CountNonZero(triplets);
    //cout << "nonzero = " << nonzero;
    SparseMatrix globalK(triplets.size());

    globalK.ConvertTripletToSparse(triplets);
    cout << "new size= "<<globalK.get_size()<<"\n";
    globalK.resize();

    //globalK.Show();


    cout << "\n\n";
    globalK.SortIt();
    //globalK.Show();
    //globalK.ShowAsMatrix(2 * nodesCount);


    ApplyConstraints(globalK, constraints);
    cout << "ApplyConstraints success\n";


    unsigned long long int end_time = clock();
    cout<< "Time: "<< end_time - start_time<< " ms\n";
    t2 = clock();
    time1 = difftime(t2, t1);
    printf("%lf\n", time1);


    int nonzero = globalK.CountNonZero();
    cout << "nonzero = " << nonzero << endl;
    SparseMatrix globalK2(nonzero);
    globalK2 = globalK.DeleteZeros();

    end_time = clock();
        cout<< "Time: "<< end_time - start_time<< " ms\n";
        t2 = clock();
        time1 = difftime(t2, t1);
        printf("%lf\n", time1);


    //globalK2.Show();
    //globalK2.ShowAsMatrix(globalK2.get_size());

    //globalK2.CGM_solve(loads, displacements, loads.get_size());

    //globalK2.Show();
    //globalK2.ShowAsMatrix(2 * nodesCount);

    //displacements.Show();
    //cout<<"\n\n\n";


    //cout << "LU_solve success\n";


    //displacements.WriteToFile();
    //outfile << displacements;
//    for (int i = 0; i < displacements.get_size(); i++) {
//        outfile << displacements[i] << " ";
//        if ((i + 1) % 3 == 0 && i != 0) {
//            outfile << "\n";
//        }
//    }
//    outfile << endl;


    cout << "SUCCESS\n";




    int *ptr = new int[loads.get_size() + 1];
    int *ind = new int[globalK2.get_size()];
    double *data_csr = new double[globalK2.get_size()];
    double *result1 = new double[loads.get_size()];
    for (int i = 0; i < loads.get_size() + 1; i++) {
        //result1[i] = 0.0;
        ptr[i] = 0;
    }
    for (int i = 0; i < globalK2.get_size(); i++) {
        ind[i] = 0;
        data_csr[i] = 0.0;
    }


    globalK2.ConvertToCSR(ptr, ind, data_csr, loads.get_size());

    Prepare_CSR(globalK2.get_x(), ptr, globalK2.get_size(), loads.get_size());

//    for (int i = 0; i < loads.get_size() + 1; i++) {
//        //ptr[i]++;
//        //cout << ptr[i] << " ";
//    }

//    cout << "\n\n";
//    for (int i = 0; i < globalK2.get_size(); i++) {
//        //ind[i]++;
//        //cout << ind[i] << " ";
//    }

//    LU_GPU_SOLVE(ptr, ind, globalK2.get_data(), loads.get_size(), globalK2.get_size(), loads.get_data(), result1);



    TestCudaSolve(ptr, ind, globalK2.get_data(), loads.get_size(), globalK2.get_size(), loads.get_data(), displacements.get_data());

//    for (int i = 0; i < loads.get_size(); i++) {
//        cout << result1
//    }


    displacements.Show();

    std::vector<MyArray> Deformation;
    std::vector<MyArray> Stress;
    std::vector<double> sigma_mises;
    std::vector<double> epsilon_mises;
    CalculateStressAndDeformation(Deformation, Stress, epsilon_mises, sigma_mises, D, elements, displacements);

    MakeVTKfile(output_vtk, nodesX, nodesY, nodesZ, elements, displacements, Stress, sigma_mises, Deformation, epsilon_mises);

    return 0;
}
