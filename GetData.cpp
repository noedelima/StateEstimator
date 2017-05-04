//
//  GetData.cpp
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


void GetYbus(int m, long double* Ybus) {
    std::cout << "Entre com as impedancias do sistema\nCaso nao haja ligacao entre as barras, utilize o valor zero que sera considerado como impedancia infinita.\n\n";
    long double r=0, x=0;
    int i, j;
    for (i=1; i<m; i++) {
        for (j=1; j<m; j++) {
            if (i < j) {
                std::cout << "Digite o vador da resistencia serie entre as barras " << i << " e " << j << ": ";
                std::cin >> r;
                std::cout << "Digite o valor da reatancia serie entre as barras " << i << " e " << j << ": ";
                std::cin >> x;
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
                Ybus[i+m*i] = Ybus[i+m*i] - Ybus[j+m*i];
                Ybus[i+m*i+m*m] = Ybus[i+m*i+m*m] - Ybus[j+m*i+m*m];
            }
        }
    } //Elementos da Diagonal Principal
    return;
} //Solicita os valores de impedância em série nas linhas e calcula a matriz Ybus com a topologia da rede

// Função que lê as medidas de potência na rede e cria uma matriz de 3 dimensões, sendo 4 camadas de matrizes bidimensionais: camada zero com os valores de potência ativa, camada 1 com os valores de potência reativa, camada 2 com os valores de incerteza na leitura da potência ativa e camada 3 com os valores de incerteza na leitura da potência reativa

void GetMeansPower(int m, long double *Sz) {
    std::cout << "Matriz de Potencias\n";
    std::cout << "Entre com as Potencias medidas no sistema\nCaso nao haja medicao, utilize o valor zero.\n\n";
    long double Pz=0, Qz=0, Rpz, Rqz;
    int i, j;
    for (i=1; i<m; i++) {
        for (j=1; j<m; j++) {
            if (i != j) {
                std::cout << "Digite o fluxo de potencia ativa medido da barra " << i << " para a barra " << j << ": ";
                std::cin >> Pz;
                if (Pz==0) {Rpz=0;}
                else {
					std::cout << "Digite a incerteza dessa medicao: ";
					std::cin >> Rpz;
                }
				std::cout << "Digite o fluxo de potencia reativa medido da barra " << i << " para a barra " << j << ": ";
				std::cin >> Qz;
                if (Qz==0) {Rqz=0;}
                else {
					std::cout << "Digite a incerteza dessa medicao: ";
					std::cin >> Rqz;
                }
                Sz[j+m*i] = Pz; //Potência ativa
                Sz[j+m*i+m*m] = Qz; //Potência reativa
                Sz[j+m*i+2*m*m] = Rpz; //Incerteza associada ao medidor de potência ativa
                Sz[j+m*i+3*m*m] = Rqz; //Incerteza associada ao medidor de potência reativa
            } //Solicita as potências medidas nas linhas
            if (i == j) {
				std::cout << "Digite o balanco de potencia ativa liquida medido na barra " << i << ": ";
				std::cin >> Pz;
                if (Pz==0) {Rpz=0;}
                else {
					std::cout << "Digite a incerteza dessa medicao: ";
					std::cin >> Rpz;
                }
				std::cout << "Digite o balanco de potencia reativa liquida medido na barra " << i << ": ";
				std::cin >> Qz;
                if (Qz==0) {Rqz=0;}
                else {
					std::cout << "Digite a incerteza dessa medicao: ";
					std::cin >> Rqz;
                }
                Sz[j+m*i] = Pz; //Potência ativa
                Sz[j+m*i+m*m] = Qz; //Potência reativa
                Sz[j+m*i+2*m*m] = Rpz; //Incerteza associada ao medidor de potência ativa
                Sz[j+m*i+3*m*m] = Rqz; //Incerteza associada ao medidor de potência reativa
            } //Solicita as potências medidas nas barras
        }
    }
    return;
} //Solicita os valores de potência lido nos medidores da rede

//Função que lê os valores de módulo de tensão nas barras da rede e cria uma matriz de duas colunas, sendo uma com o valor da leitura e outra com a incerteza associada à leitura

void GetMeansVoltageMag(int m, long double* Vz) {
    for (int i=1; i<m; i++) {
        long double vz, Rvz;
		std::cout << "Entre com a medida do modulo de tensao na barra " << i << ": ";
		std::cin >> vz;
        if (vz==0) {Rvz=0;}
        else {
			std::cout << "Digite a incerteza dessa medicao: ";
			std::cin >> Rvz;
        }
            Vz[i] = vz;
        Vz[i+m] = Rvz;
    }
    return;
} //Solicita os valores de magnitude das tensões lidos nas barras
