//
//  MatrixFunctions.cpp
//  Estimador de Estado
//
//  Created by Noé de Lima Bezerra on 31/12/00.
//  Copyright © 2000 Noé de Lima Bezerra. All rights reserved.
//

#include "AuxFunctions.hpp"
#include <iostream>
#include <cmath>
#include <complex>
#include <new>
using namespace std;

void Transp(long double* A, int line, int row, long double* T) {
    for (int i=0; i<line; i++) {
        for (int j=0; j<row; j++) {
            T[i+line*j] = A[j+row*i];
        }
    }
    return;
} //Calcula a transposta de uma matriz

void Product(long double* A, int aline, int arow, long double* B, int bline, int brow, long double* P) {
    if (arow!=bline) {
        cout << "\nMatriz com erro de dimensao";
        return;
    }
    for (int i=0; i<aline; i++) {
        for (int j=0; j<brow; j++) {
            P[j+brow*i]=0;
            for (int k=1; k<arow; k++) {
                P[j+brow*i] = P[j+brow*i] + A[k+arow*i]*B[j+brow*k];
            }
        }
    }
    return;
} //Calcula o produto de duas matrizes

int InvertMatrix(int n, long double* X, long double* Y) {
    long double *I;
	I = new long double [(n)*(n)]; // Matriz 2 dim

    for (int i=1; i<n; i++) {
        for (int j=1; j<n; j++) {
            I[j+n*i] = X[j+n*i];
            if (i==j) {Y[j+n*i]=1;}
            else {Y[j+n*i]=0;}
        }
    } // Inicializa a matriz identidade com a matriz de entrada a ser invertida
    /*if (Det(n, I)==0) {
        //cout << "Matriz singular, nao e possivel inverter\n\n";
        return 1;
    } // Testa se a inversão é possível */
    for (int i=1; i<n; i++) {
        int k=0;
        while (I[i+n*i]==0) {
            k=k+1;
            if (k!=i) {
                for (int row=1; row<n; row++) {
                    I[row+n*i] += I[row+n*k];
                    Y[row+n*i] += Y[row+n*k];
                }
            }
        } // Evita divisão por zero somando outra linha não nula
        long double a = I[i+n*i];
        for (int row=1; row<n; row++) {
            I[row+n*i] = I[row+n*i]/a;
            Y[row+n*i] = Y[row+n*i]/a;
        } // Divide a linha pelo elemento da diagonal principal
        for (int line=1; line<n; line++) {
            if (line != i) {
                a = I[i+n*line];
                for (int row=1; row<n; row++) {
                    I[row+n*line] += -a*I[row+n*i];
                    Y[row+n*line] += -a*Y[row+n*i];
                }
            }
        } // transforma os demais elementos da coluna em zero
    }
	delete[] I;

    return 0;
} //Calcula a inversa de uma matriz

void Diag(int M, long double* sigma, long double* R) {
    M=M+1;
    for (int i=0; i<M; i++) {
        for (int j=0; j<M; j++) {
            R[j+M*i] = 0;
            if (i==j) {
                R[j+M*i] = sigma[i];
            }
        }
    }
    return;
} //Retorna uma matriz com os elementos do vetor R na diagonal e os demais elementos nulos

long double Det(int m, long double* A) {

	long double det=0;
    if (m==2) {
        det = A[3];
    }
    else {
        int n = m-1;
        for (int j=1; j<m; j++) {
            if (A[j+m]!=0) {
                long double *M;
				M = new long double [(n)*(n)]; // Matriz 2 dim

                int line=0;
                for (int l=2; l<m; l++) {
                    ++line;
                    int row=0;
                    for (int r=1; r<m; r++) {
                        if (r!=j) {
                            ++row;
                            M[row+n*line] = A[r+m*l];
                        }
                    }
                }
                det += pow(-1,j+1)*A[j+m]*Det(n, M);
				delete[] M;
            }
        }
    }
    return det;
}
