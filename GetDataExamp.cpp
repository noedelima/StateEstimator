//
//  GetDataExamp.cpp
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

// Todos os arrays de entrada, vetores e matrizes, são enviados como endereço do elemento 0 e recebidos como ponteiros
// Os índices dos arrays multidimensionais são calculados com base na disposição em memória utilizando lógica de base numérica
// Por exemplo, o elemento Matriz[i][j][k] será Matriz[k+n*j+n*n*k], onde n é a ordem da matriz


//Valores de entrada para simular o exemplo 2.2 do livro Power System State Estimator, Ali Abur & Antonio Exposito
void GetYbus1(int m, long double* Ybus) {
    //cout << "Entre com as impedâncias do sistema\nCaso não haja ligação entre as barras, utilize o valor zero que será considerado como impedância infinita.\n\n";
    long double r=0, x=0;
    int i, j;
    for (i=1; i<m; i++) {
        for (j=1; j<m; j++) {
            if (i < j) {
                //cout << "Digite o vador da resistência série entre as barras " << i << " e " << j << ":\n";
                //cin >> r;
                //cout << "Digite o valor da reatância série entre as barras " << i << " e " << j << ":\n";
                //cin >> x;
                if (i==1&&j==2) {
                    r = 0.01;
                    x = 0.03;
                }
                else if (i==1&&j==3) {
                    r = 0.02;
                    x = 0.05;
                }
                else if (i==2&&j==3) {
                    r = 0.03;
                    x = 0.08;
                }
                else {
                    r = 0;
                    x = 0;
                }
                complex<long double> y(r, x);
                if (r==0 && x==0) {
                    y = complex<long double> (0,0);
                }
                else {y = pow(y, -1);}
                Ybus[j+m*i] = -real(y);
                Ybus[j+m*i+m*m] = -imag(y);
            }
            else if (i > j) {
                Ybus[j+m*i] = Ybus[i+m*j];
                Ybus[j+m*i+m*m] = Ybus[i+m*j+m*m];
            }
            if (i == j) {
                Ybus[j+m*i] = 0;
                Ybus[j+m*i+m*m] = 0;
            }
        }
    } //Elementos fora da diagonal principal
    for (i=1; i<m; i++) {
        for (j=1; j<m; j++) {
            if (i != j) {
                Ybus[i+m*i] -= Ybus[j+m*i];
                Ybus[i+m*i+m*m] -= Ybus[j+m*i+m*m];
            }
        }
    } //Elementos da Diagonal Principal
    return;
} //Solicita os valores de impedância em série nas linhas e calcula a matriz Ybus com a topologia da rede

// Função que lê as medidas de potência na rede e cria uma matriz de 3 dimensões, sendo 4 camadas de matrizes bidimensionais: camada zero com os valores de potência ativa, camada 1 com os valores de potência reativa, camada 2 com os valores de incerteza na leitura da potência ativa e camada 3 com os valores de incerteza na leitura da potência reativa

void GetMeansPower1(int m, long double *Sz) {
    //cout << "Matriz de Potências\n";
    //cout << "Entre com as Potências medidas no sistema\nCaso não haja medição, utilize o valor zero.\n\n";
    long double Pz=0, Qz=0, Rpz, Rqz;
    int i, j;
    for (i=1; i<m; i++) {
        for (j=1; j<m; j++) {
            if (i != j) {
                //cout << "Digite o fluxo de potência ativa medido da barra " << i << " para a barra " << j << ": ";
                Pz = 0; //cin >> Pz;
                if (Pz==0) {Rpz=0;}
                else {
                    //cout << "Digite a incerteza dessa medição: ";
                    cin >> Rpz;
                }
                //cout << "Digite o fluxo de potência reativa medido da barras " << i << " para a barra " << j << ": ";
                Qz = 0; //cin >> Qz;
                if (Qz==0) {Rqz=0;}
                else {
                    //cout << "Digite a incerteza dessa medição: ";
                    cin >> Rqz;
                }
                Sz[j+m*i] = Pz; //Potência ativa
                Sz[j+m*i+m*m] = Qz; //Potência reativa
                Sz[j+m*i+2*m*m] = Rpz; //Incerteza associada ao medidor de potência ativa
                Sz[j+m*i+3*m*m] = Rqz; //Incerteza associada ao medidor de potência reativa
            } //Solicita as potências medidas nas linhas
            if (i == j) {
                //cout << "Digite o balanço de potência ativa líquida medido na barra " << i << ": ";
                Pz = 0; //cin >> Pz;
                if (Pz==0) {Rpz=0;}
                else {
                    cout << "Digite a incerteza dessa medicao: ";
                    cin >> Rpz;
                }
                //cout << "Digite o balanço de potência reativa líquida medido na barra " << i << ": ";
                Qz = 0; //cin >> Qz;
                if (Qz==0) {Rqz=0;}
                else {
                    //cout << "Digite a incerteza dessa medição: ";
                    cin >> Rqz;
                }
                Sz[j+m*i] = Pz; //Potência ativa
                Sz[j+m*i+m*m] = Qz; //Potência reativa
                Sz[j+m*i+2*m*m] = Rpz; //Incerteza associada ao medidor de potência ativa
                Sz[j+m*i+3*m*m] = Rqz; //Incerteza associada ao medidor de potência reativa
            } //Solicita as potências medidas nas barras
        }
    }

    //Elemento 1-2
    Sz[2+m*1+0*m*m] = 0.888;//P
    Sz[2+m*1+1*m*m] = 0.568;//Q
    Sz[2+m*1+2*m*m] = 0.008;//∆P
    Sz[2+m*1+3*m*m] = 0.008;//∆Q
    
	//Elemento 1-3
    Sz[3+m*1+0*m*m] = 1.173;//P
    Sz[3+m*1+1*m*m] = 0.663;//Q
    Sz[3+m*1+2*m*m] = 0.008;//∆P
    Sz[3+m*1+3*m*m] = 0.008;//∆Q
    
	//Elemento 2-2
    Sz[2+m*2+0*m*m] = -0.501;//P
    Sz[2+m*2+1*m*m] = -0.286;//Q
    Sz[2+m*2+2*m*m] = 0.01;//∆P
    Sz[2+m*2+3*m*m] = 0.01;//∆Q

    return;
} //Solicita os valores de potência lido nos medidores da rede

//Função que lê os valores de módulo de tensão nas barras da rede e cria uma matriz de duas colunas, sendo uma com o valor da leitura e outra com a incerteza associada à leitura

void GetMeansVoltageMag1(int m, long double* Vz) {
    for (int i=1; i<m; i++) {
        long double vz, Rvz;
        //cout << "Entre com a medida do módulo de tensão na barra " << i << ": ";
        vz=0; //cin >> vz;
        if (vz==0) {Rvz=0;}
        else {
            cout << "Digite a incerteza dessa medicao: ";
            cin >> Rvz;
        }
            Vz[i] = vz;
        Vz[i+m] = Rvz;
    }
    
    //Barra 1
    Vz[1+0*m] = 1.006;//V
    Vz[1+1*m] = 0.004;//∆V
    //Barra 2
    Vz[2+0*m] = 0.968;//V
    Vz[2+1*m] = 0.004;//∆V
    
    return;
} //Solicita os valores de magnitude das tensões lidos nas barras

