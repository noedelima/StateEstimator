//
//  SystemState.cpp
//  Estimador de Estado
//
//  Created by Noé de Lima Bezerra on 01/01/01.
//  Copyright © 2001 Noé de Lima Bezerra. All rights reserved.
//

#include "AuxFunctions.hpp"
#include <iostream>
#include <cmath>
#include <complex>
#include <new>
using namespace std;

// Esta função estima os fluxos de potência do sistema com base nas variáveis de estado
// As variáveis de estado são os módulos de tensão nas barras e os ângulos de fase

// Todos os arrays de entrada, vetores e matrizes, são enviados como endereço do elemento 0 e recebidos como ponteiros
// Os índices dos arrays multidimensionais são calculados com base na disposição em memória utilizando lógica de base numérica
// Por exemplo, o elemento Matriz[i][j][k] será Matriz[k+n*j+n*n*k], onde n é a ordem da matriz

void CalcPower(int m, long double* Ybus, long double* S, long double* V, long double* Theta) {
    long double G, B;
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            if (i==j) {
                S[j+m*i] = 0;
                S[j+m*i+m*m] = 0;
                for (int k=1; k<m; k++) {
                    G = -Ybus[k+m*i];
                    B = -Ybus[k+m*i+m*m];
                    S[j+m*i] = S[j+m*i] - V[i]*V[k]*(G*cos(Theta[i]-Theta[k])+B*sin(Theta[i]-Theta[k])); //Potência ativa
                    S[j+m*i+m*m] = S[j+m*i+m*m] - V[i]*V[k]*(G*sin(Theta[i]-Theta[k])-B*cos(Theta[i]-Theta[k])); //Potência reativa
                } //Potência líquida nas barras
            }
            else {
                G = -Ybus[j+m*i];
                B = -Ybus[j+m*i+m*m];
                S[j+m*i] = pow(V[i], 2)*G-V[i]*V[j]*(G*cos(Theta[i]-Theta[j])+B*sin(Theta[i]-Theta[j])); //Potência ativa
                S[j+m*i+m*m] = -pow(V[i], 2)*B-V[i]*V[j]*(G*sin(Theta[i]-Theta[j])-B*cos(Theta[i]-Theta[j])); //Potência reativa
            } //Fluxos de potência nas linhas
        }
    }
    return;
} // Calcula a matriz com as potências injetadas/consumidas nas barras e os fluxos de potências injetados nas linhas, tanto de potência ativa quanto de potência reativa
