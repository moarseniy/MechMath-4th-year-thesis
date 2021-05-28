
#include <iostream>
#include "femfunc.h"
#include <vector>

using namespace std;

void FindNonZeroNumbers(int *element0, int *element1, int *element2, int *element3, int elementCount,
                        float *nodesXcuda, float *nodesYcuda, float *nodesZcuda, int nodesCount) {
    int k = 0;
    for (int i = 0; i < nodesCount; i++) {
        for (int j = 0; j < elementCount; j++) {
            if (i == element0[j] || i == element1[j] || i == element2[j] || i == element3[j]) {
                k += 4 * 3;
            }
        }
    }

    cout << "kolvo nonzero = " << k << endl;
}


void FindSetsColor(std::vector<Element> &elements, std::vector<Element> &elementsColor, std::vector<int> &elementsOrderColor) {

    int k = 0, k_prev = 0, i = 0;
    const int size = elements.size();
    bool add;

    while (elementsColor.size() != size) {
        while (i < elements.size()) {
            add = true;
            for (int j = k_prev; j < k + k_prev; j++) {
                for (int t = 0; t < 4; t++) {
                    for (int p = 0; p < 4; p++) {
                        if (elements[i].nodesIds[t] == elementsColor[j].nodesIds[p] ||
                                elements[i].nodesIds[t] == elementsColor[j].nodesIds[p] ||
                                elements[i].nodesIds[t] == elementsColor[j].nodesIds[p] ||
                                elements[i].nodesIds[t] == elementsColor[j].nodesIds[p]) {
                            add = false;
                        }
                    }
                }
            }
            if (add == true) {
                elementsColor.push_back(elements[i]);
//                cout << elements[i].nodesIds[0] << " "
//                        << elements[i].nodesIds[1] << " "
//                        << elements[i].nodesIds[2] << " "
//                        << elements[i].nodesIds[3] << endl;
                elements.erase(elements.begin() + i);
                k++;
            } else {
                i++;
            }
        }
        cout << k << " ";
        elementsOrderColor.push_back(k);
        k_prev += k;
        k = 0;
        i = 0;
    }
}


void Element::FindSparseSize(std::vector<couple> &Sparse) {
    int SIZE = Sparse.size();

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            couple c11(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 0);
            couple c12(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 1);
            couple c13(3 * nodesIds[i] + 0, 3 * nodesIds[j] + 2);

            couple c21(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 0);
            couple c22(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 1);
            couple c23(3 * nodesIds[i] + 1, 3 * nodesIds[j] + 2);

            couple c31(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 0);
            couple c32(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 1);
            couple c33(3 * nodesIds[i] + 2, 3 * nodesIds[j] + 2);

            bool a11 = false, a12 = false, a13 = false, a21 = false, a22 = false, a23 = false,
                    a31 = false, a32 = false, a33 = false;
            if (SIZE == 0) {
                Sparse.push_back(c11);
                Sparse.push_back(c12);
                Sparse.push_back(c13);
                Sparse.push_back(c21);
                Sparse.push_back(c22);
                Sparse.push_back(c23);
                Sparse.push_back(c31);
                Sparse.push_back(c32);
                Sparse.push_back(c33);
            } else {
            for (int k = 0; k < SIZE; k++) {
                if (Sparse[k].x == c11.x && Sparse[k].y == c11.y) {
                    a11 = true;
                }
                if (Sparse[k].x == c12.x && Sparse[k].y == c12.y) {
                    a12 = true;
                }
                if (Sparse[k].x == c13.x && Sparse[k].y == c13.y) {
                    a13 = true;
                }
                if (Sparse[k].x == c21.x && Sparse[k].y == c21.y) {
                    a21 = true;
                }
                if (Sparse[k].x == c22.x && Sparse[k].y == c22.y) {
                    a22 = true;
                }
                if (Sparse[k].x == c23.x && Sparse[k].y == c23.y) {
                    a23 = true;
                }
                if (Sparse[k].x == c31.x && Sparse[k].y == c31.y) {
                    a31 = true;
                }
                if (Sparse[k].x == c32.x && Sparse[k].y == c32.y) {
                    a32 = true;
                }
                if (Sparse[k].x == c33.x && Sparse[k].y == c33.y) {
                    a33 = true;
                }
                if (k == SIZE - 1) {
                    if (!a11) { Sparse.push_back(c11); }
                    if (!a12) { Sparse.push_back(c12); }
                    if (!a13) { Sparse.push_back(c13); }
                    if (!a21) { Sparse.push_back(c21); }
                    if (!a22) { Sparse.push_back(c22); }
                    if (!a23) { Sparse.push_back(c23); }
                    if (!a31) { Sparse.push_back(c31); }
                    if (!a32) { Sparse.push_back(c32); }
                    if (!a33) { Sparse.push_back(c33); }
                }
            }
            }
        }
    }

}

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
    float determinant = C.det(4);
    //cout << determinant << endl;

    //K = B.transpose() * D * B * std::abs(C.det()) / 2.0;

    temp_B = B;
    //temp_B.Show();

    temp_B.transpose2();
    //temp_B.Show();
    //cout << B(0, 0) << "-" << temp_B(0, 0) << " ";
    //cout << C(0,1) << " ";
    //cout << IC(0,1) << " ";


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
                //nodesStructure[trplt11.get_y() + trplt11.get_x() * n]++;
            }
            if (trplt12.get_value() != 0.0) {
                triplets.push_back(trplt12);
                //nodesStructure[trplt12.get_y() + trplt12.get_x() * n]++;
            }
            if (trplt13.get_value() != 0.0) {
                triplets.push_back(trplt13);
                //nodesStructure[trplt13.get_y() + trplt13.get_x() * n]++;
            }
            if (trplt21.get_value() != 0.0) {
                triplets.push_back(trplt21);
                //nodesStructure[trplt21.get_y() + trplt21.get_x() * n]++;
            }
            if (trplt22.get_value() != 0.0) {
                triplets.push_back(trplt22);
                //nodesStructure[trplt22.get_y() + trplt22.get_x() * n]++;
            }
            if (trplt23.get_value() != 0.0) {
                triplets.push_back(trplt23);
                //nodesStructure[trplt23.get_y() + trplt23.get_x() * n]++;
            }
            if (trplt31.get_value() != 0.0) {
                triplets.push_back(trplt31);
                //nodesStructure[trplt31.get_y() + trplt31.get_x() * n]++;
            }
            if (trplt32.get_value() != 0.0) {
                triplets.push_back(trplt32);
                //nodesStructure[trplt32.get_y() + trplt32.get_x() * n]++;
            }
            if (trplt33.get_value() != 0.0) {
                triplets.push_back(trplt33);
                //nodesStructure[trplt33.get_y() + trplt33.get_x() * n]++;
            }
        }
    }
}
