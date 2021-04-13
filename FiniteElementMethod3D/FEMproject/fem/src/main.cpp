#include "init.h"
#include "Linal2.h"

#include <iostream>
#include <vector>
#include <fstream>
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
        UX = 1 << 0,
        UY = 1 << 1,
        UXY = UX | UY
    };
    int node;
    Type type;
};


void Element::CalculateStiffnessMatrix(Matrix& D, std::vector<Triplet>& triplets, MyArray& nodesX, MyArray& nodesY) {
    MyArray x(3), y(3);
    x[0] = nodesX[nodesIds[0]]; x[1] = nodesX[nodesIds[1]]; x[2] = nodesX[nodesIds[2]];
    y[0] = nodesY[nodesIds[0]]; y[1] = nodesY[nodesIds[1]]; y[2] = nodesY[nodesIds[2]];
    //x.Show();
    //y.Show();

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

void ApplyConstraints(Matrix& K, const std::vector<Constraint>& constraints)
{
    std::vector<int> indicesToConstraint;

    for (std::vector<Constraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it) {
        if (it->type & Constraint::UX) {
            indicesToConstraint.push_back(2 * it->node + 0);
        }
        if (it->type & Constraint::UY) {
            indicesToConstraint.push_back(2 * it->node + 1);
        }
    }


    for (int i = 0; i < K.get_row(); i++) {
        for (int j = 0; j < K.get_col(); j++) {
            for (std::vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit) {
                K(i, j) = SetConstraints(i, j, K(i, j), *idit);
            }
        }
    }
}


int main(void) {
    //callCudaKernel();
    fstream infile, outfile;
    infile.open("D:/FEMproject/fem/input2.txt", fstream::in);
    outfile.open("D:/results/output.txt", fstream::out);

    double poissonRatio, youngModulus;
    infile >> poissonRatio >> youngModulus;
    Matrix D(3, 3);

    D(0,0) = 1.0;			D(0, 1) = poissonRatio;	D(0, 2) = 0.0;
    D(1, 0) = poissonRatio;	D(1, 1) = 1.0; 			D(1, 2) = 0.0;
    D(2, 0) = 0.0;        	D(2, 1) = 0.0;        	D(2, 2) = (1.0 - poissonRatio) / 2.0;

    D.scale(youngModulus / (1.0 - pow(poissonRatio, 2.0)));

    int nodesCount;
    infile >> nodesCount;

    MyArray nodesX(nodesCount);
    MyArray nodesY(nodesCount);

    for (int i = 0; i < nodesCount; ++i) {
        infile >> nodesX[i] >> nodesY[i];
    }
    //nodesX.Show();
    //cout<<"\n";
    //nodesY.Show();
    std::vector<Element>   	elements;
    std::vector<Constraint>	constraints;

    int elementCount;
    infile >> elementCount;
    for (int i = 0; i < elementCount; ++i) {
        Element element;
        infile >> element.nodesIds[0] >> element.nodesIds[1] >> element.nodesIds[2];
        elements.push_back(element);
    }

    int constraintCount;
    infile >> constraintCount;

    for (int i = 0; i < constraintCount; ++i) {
        Constraint constraint;
        int type;
        infile >> constraint.node >> type;
        constraint.type = static_cast<Constraint::Type>(type);
        constraints.push_back(constraint);
    }

    MyArray loads(2 * nodesCount);

    int loadsCount;
    infile >> loadsCount;

    for (int i = 0; i < loadsCount; ++i) {
        int node;
        float x, y;
        infile >> node >> x >> y;
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
    cout<<triplets.size()<<endl;
    int sparse_size = triplets.size();


    Matrix K_matrix(2 * nodesCount);

    //int nonzero;
    //nonzero = CountNonZero(triplets);
    //cout << "nonzero = " << nonzero;
    SparseMatrix globalK(sparse_size);

    globalK.ConvertTripletToSparse(triplets);
    cout << "new size= "<<globalK.get_size()<<"\n";
    globalK.Show();

    cout << "\n\n";
    globalK.SortIt();
    globalK.Show();

    globalK.ConvertToMatrix(K_matrix);
    //K_matrix.Show();

    ApplyConstraints(K_matrix, constraints);
    K_matrix.Show();
    cout << "ApplyConstraints success\n";

    loads.Show();
    MyArray displacements(loads.get_size());
    K_matrix.LU_solve(K_matrix, loads, displacements, loads.get_size());

    cout << "LU_solve success\n";

    //displacements.WriteToFile();
    //outfile << displacements;
    for (int i = 0; i < displacements.get_size(); i++) {
        outfile << displacements[i] << " ";
        if (i % 2 != 0 && i != 0) {
            outfile << "\n";
        }
    }
    outfile << endl;

    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it) {
        Matrix delta(6, 1);
        //std::cout<<2*it->nodesIds[0]<<" "<<2*it->nodesIds[1]<<" "<<2*it->nodesIds[2]<<endl;
        delta(0, 0) = displacements[2 * it->nodesIds[0]];
        delta(1, 0) = displacements[2 * it->nodesIds[0] + 1];
        delta(2, 0) = displacements[2 * it->nodesIds[1]];
        delta(3, 0) = displacements[2 * it->nodesIds[1] + 1];
        delta(4, 0) = displacements[2 * it->nodesIds[2]];
        delta(5, 0) = displacements[2 * it->nodesIds[2] + 1];
        //D * it->B * delta;
        Matrix sigma(3, 1);
        sigma = D.Product(it->B);
        sigma = sigma.Product(delta);
        double sigma_mises = sqrt(sigma(0, 0) * sigma(0, 0) - sigma(0, 0)
            * sigma(1, 0) + sigma(1, 0) * sigma(1, 0) + 3.0 * sigma(2, 0) * sigma(2, 0));

        outfile << sigma_mises << std::endl;
    }

    cout << "SUCCESS\n";

    return 0;
}
