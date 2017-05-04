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


//Valores de entrada para o exemplo do IEEE 14 barras com redundância de 1,5
void GetYbus14bus350(int m, long double* Ybus) {
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
                    r = 0.01938;
                    x = 0.05917;
                }
                else if (i==1&&j==5) {
                    r = 0.05403;
                    x = 0.22304;
                }
                else if (i==2&&j==3) {
                    r = 0.04699;
                    x = 0.19797;
                }
                else if (i==2&&j==4) {
                    r = 0.05811;
                    x = 0.17632;
                }
                else if (i==2&&j==5) {
                    r = 0.05695;
                    x = 0.17388;
                }
                else if (i==3&&j==4) {
                    r = 0.06701;
                    x = 0.17103;
                }
                else if (i==4&&j==5) {
                    r = 0.01335;
                    x = 0.04211;
                }
                else if (i==4&&j==7) {
                    r = 0;
                    x = 0.20912;
                }
                else if (i==4&&j==9) {
                    r = 0;
                    x = 0.55618;
                }
                else if (i==5&&j==6) {
                    r = 0;
                    x = 0.25202;
                }
                else if (i==6&&j==11) {
                    r = 0.09498;
                    x = 0.1989;
                }
                else if (i==6&&j==12) {
                    r = 0.12291;
                    x = 0.25581;
                }
                else if (i==6&&j==13) {
                    r = 0.06615;
                    x = 0.13027;
                }
                else if (i==7&&j==8) {
                    r = 0;
                    x = 0.17615;
                }
                else if (i==7&&j==9) {
                    r = 0;
                    x = 0.11001;
                }
                else if (i==9&&j==10) {
                    r = 0.03181;
                    x = 0.0845;
                }
                else if (i==9&&j==14) {
                    r = 0.12711;
                    x = 0.27038;
                }
                else if (i==10&&j==11) {
                    r = 0.08205;
                    x = 0.19207;
                }
                else if (i==12&&j==13) {
                    r = 0.22092;
                    x = 0.19988;
                }
                else if (i==13&&j==14) {
                    r = 0.17093;
                    x = 0.34802;
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
                Ybus[i+m*i] = Ybus[i+m*i] - Ybus[j+m*i];
                Ybus[i+m*i+m*m] = Ybus[i+m*i+m*m] - Ybus[j+m*i+m*m];
            }
        }
    } //Elementos da Diagonal Principal
    return;
} //Solicita os valores de impedância em série nas linhas e calcula a matriz Ybus com a topologia da rede

// Função que lê as medidas de potência na rede e cria uma matriz de 3 dimensões, sendo 4 camadas de matrizes bidimensionais: camada zero com os valores de potência ativa, camada 1 com os valores de potência reativa, camada 2 com os valores de incerteza na leitura da potência ativa e camada 3 com os valores de incerteza na leitura da potência reativa

void GetMeansPower14bus350(int m, long double *Sz) {
	//cout << "Matriz de Potências\n";
	//cout << "Entre com as Potências medidas no sistema\nCaso não haja medição, utilize o valor zero.\n\n";
	long double Pz = 0, Qz = 0, Rpz, Rqz;
	int i, j;
	for (i = 1; i<m; i++) {
		for (j = 1; j<m; j++) {
			if (i != j) {
				//cout << "Digite o fluxo de potência ativa medido da barra " << i << " para a barra " << j << ": ";
				Pz = 0; //cin >> Pz;
				if (Pz == 0) { Rpz = 0; }
				else {
					//cout << "Digite a incerteza dessa medição: ";
					cin >> Rpz;
				}
				//cout << "Digite o fluxo de potência reativa medido da barras " << i << " para a barra " << j << ": ";
				Qz = 0; //cin >> Qz;
				if (Qz == 0) { Rqz = 0; }
				else {
					//cout << "Digite a incerteza dessa medição: ";
					cin >> Rqz;
				}
				Sz[j + m*i] = Pz; //Potência ativa
				Sz[j + m*i + m*m] = Qz; //Potência reativa
				Sz[j + m*i + 2 * m*m] = Rpz; //Incerteza associada ao medidor de potência ativa
				Sz[j + m*i + 3 * m*m] = Rqz; //Incerteza associada ao medidor de potência reativa
			} //Solicita as potências medidas nas linhas
			if (i == j) {
				//cout << "Digite o balanço de potência ativa líquida medido na barra " << i << ": ";
				Pz = 0; //cin >> Pz;
				if (Pz == 0) { Rpz = 0; }
				else {
					cout << "Digite a incerteza dessa medicao: ";
					cin >> Rpz;
				}
				//cout << "Digite o balanço de potência reativa líquida medido na barra " << i << ": ";
				Qz = 0; //cin >> Qz;
				if (Qz == 0) { Rqz = 0; }
				else {
					//cout << "Digite a incerteza dessa medição: ";
					cin >> Rqz;
				}
				Sz[j + m*i] = Pz; //Potência ativa
				Sz[j + m*i + m*m] = Qz; //Potência reativa
				Sz[j + m*i + 2 * m*m] = Rpz; //Incerteza associada ao medidor de potência ativa
				Sz[j + m*i + 3 * m*m] = Rqz; //Incerteza associada ao medidor de potência reativa
			} //Solicita as potências medidas nas barras
		}
	}

	//Elemento 1-1
	Sz[1 + m * 1 + 0 * m*m] = 1;//P
	Sz[1 + m * 1 + 1 * m*m] = 1;//Q
	Sz[1 + m * 1 + 2 * m*m] = 0.008;//∆P
	Sz[1 + m * 1 + 3 * m*m] = 0.008;//∆Q

			//Elemento 1-2
	Sz[2 + m * 1 + 0 * m*m] = 1;//P
	Sz[2 + m * 1 + 1 * m*m] = 1;//Q
	Sz[2 + m * 1 + 2 * m*m] = 0.008;//∆P
	Sz[2 + m * 1 + 3 * m*m] = 0.008;//∆Q

				//Elemento 2-1
	Sz[1 + m * 2 + 0 * m*m] = -1;//P
	Sz[1 + m * 2 + 1 * m*m] = -1;//Q
	Sz[1 + m * 2 + 2 * m*m] = 0.008;//∆P
	Sz[1 + m * 2 + 3 * m*m] = 0.008;//∆Q
									
			//Elemento 1-5
	Sz[5 + m * 1 + 0 * m*m] = 1;//P
	Sz[5 + m * 1 + 1 * m*m] = 1;//Q
	Sz[5 + m * 1 + 2 * m*m] = 0.008;//∆P
	Sz[5 + m * 1 + 3 * m*m] = 0.008;//∆Q


	//Elemento 2-2
	Sz[2 + m * 2 + 0 * m*m] = 1;//P
	Sz[2 + m * 2 + 1 * m*m] = 1;//Q
	Sz[2 + m * 2 + 2 * m*m] = 0.008;//∆P
	Sz[2 + m * 2 + 3 * m*m] = 0.008;//∆Q

			//Elemento 2-3
	Sz[3 + m * 2 + 0 * m*m] = 1;//P
	Sz[3 + m * 2 + 1 * m*m] = 1;//Q
	Sz[3 + m * 2 + 2 * m*m] = 0.008;//∆P
	Sz[3 + m * 2 + 3 * m*m] = 0.008;//∆Q

				//Elemento 3-2
	Sz[2 + m * 3 + 0 * m*m] = -1;//P
	Sz[2 + m * 3 + 1 * m*m] = -1;//Q
	Sz[2 + m * 3 + 2 * m*m] = 0.008;//∆P
	Sz[2 + m * 3 + 3 * m*m] = 0.008;//∆Q

			//Elemento 2-4
	Sz[4 + m * 2 + 0 * m*m] = 1;//P
	Sz[4 + m * 2 + 1 * m*m] = 1;//Q
	Sz[4 + m * 2 + 2 * m*m] = 0.008;//∆P
	Sz[4 + m * 2 + 3 * m*m] = 0.008;//∆Q

				//Elemento 4-2
	Sz[2 + m * 4 + 0 * m*m] = -1;//P
	Sz[2 + m * 4 + 1 * m*m] = -1;//Q
	Sz[2 + m * 4 + 2 * m*m] = 0.008;//∆P
	Sz[2 + m * 4 + 3 * m*m] = 0.008;//∆Q

			//Elemento 2-5
	Sz[5 + m * 2 + 0 * m*m] = 1;//P
	Sz[5 + m * 2 + 1 * m*m] = 1;//Q
	Sz[5 + m * 2 + 2 * m*m] = 0.008;//∆P
	Sz[5 + m * 2 + 3 * m*m] = 0.008;//∆Q

				//Elemento 5-2
	Sz[2 + m * 5 + 0 * m*m] = -1;//P
	Sz[2 + m * 5 + 1 * m*m] = -1;//Q
	Sz[2 + m * 5 + 2 * m*m] = 0.008;//∆P
	Sz[2 + m * 5 + 3 * m*m] = 0.008;//∆Q


	//Elemento 3-3
	Sz[3 + m * 3 + 0 * m*m] = 1;//P
	Sz[3 + m * 3 + 1 * m*m] = 1;//Q
	Sz[3 + m * 3 + 2 * m*m] = 0.008;//∆P
	Sz[3 + m * 3 + 3 * m*m] = 0.008;//∆Q

			//Elemento 3-4
	Sz[4 + m * 3 + 0 * m*m] = 1;//P
	Sz[4 + m * 3 + 1 * m*m] = 1;//Q
	Sz[4 + m * 3 + 2 * m*m] = 0.008;//∆P
	Sz[4 + m * 3 + 3 * m*m] = 0.008;//∆Q


	//Elemento 4-4
	Sz[4 + m * 4 + 0 * m*m] = 1;//P
	Sz[4 + m * 4 + 1 * m*m] = 1;//Q
	Sz[4 + m * 4 + 2 * m*m] = 0.008;//∆P
	Sz[4 + m * 4 + 3 * m*m] = 0.008;//∆Q

			//Elemento 4-5
	Sz[5 + m * 4 + 0 * m*m] = 1;//P
	Sz[5 + m * 4 + 1 * m*m] = 1;//Q
	Sz[5 + m * 4 + 2 * m*m] = 0.008;//∆P
	Sz[5 + m * 4 + 3 * m*m] = 0.008;//∆Q

				//Elemento 5-4
	Sz[4 + m * 5 + 0 * m*m] = -1;//P
	Sz[4 + m * 5 + 1 * m*m] = -1;//Q
	Sz[4 + m * 5 + 2 * m*m] = 0.008;//∆P
	Sz[4 + m * 5 + 3 * m*m] = 0.008;//∆Q

			//Elemento 4-7
	Sz[7 + m * 4 + 0 * m*m] = 1;//P
	Sz[7 + m * 4 + 1 * m*m] = 1;//Q
	Sz[7 + m * 4 + 2 * m*m] = 0.008;//∆P
	Sz[7 + m * 4 + 3 * m*m] = 0.008;//∆Q

				//Elemento 7-4
	Sz[4 + m * 7 + 0 * m*m] = -1;//P
	Sz[4 + m * 7 + 1 * m*m] = -1;//Q
	Sz[4 + m * 7 + 2 * m*m] = 0.008;//∆P
	Sz[4 + m * 7 + 3 * m*m] = 0.008;//∆Q

			//Elemento 4-9
	Sz[9 + m * 4 + 0 * m*m] = 1;//P
	Sz[9 + m * 4 + 1 * m*m] = 1;//Q
	Sz[9 + m * 4 + 2 * m*m] = 0.008;//∆P
	Sz[9 + m * 4 + 3 * m*m] = 0.008;//∆Q

				//Elemento 9-4
	Sz[4 + m * 9 + 0 * m*m] = -1;//P
	Sz[4 + m * 9 + 1 * m*m] = -1;//Q
	Sz[4 + m * 9 + 2 * m*m] = 0.008;//∆P
	Sz[4 + m * 9 + 3 * m*m] = 0.008;//∆Q


	//Elemento 5-5
	Sz[5 + m * 5 + 0 * m*m] = 1;//P
	Sz[5 + m * 5 + 1 * m*m] = 1;//Q
	Sz[5 + m * 5 + 2 * m*m] = 0.008;//∆P
	Sz[5 + m * 5 + 3 * m*m] = 0.008;//∆Q

			//Elemento 5-6
	Sz[6 + m * 5 + 0 * m*m] = 1;//P
	Sz[6 + m * 5 + 1 * m*m] = 1;//Q
	Sz[6 + m * 5 + 2 * m*m] = 0.008;//∆P
	Sz[6 + m * 5 + 3 * m*m] = 0.008;//∆Q

				//Elemento 6-5
	Sz[5 + m * 6 + 0 * m*m] = -1;//P
	Sz[5 + m * 6 + 1 * m*m] = -1;//Q
	Sz[5 + m * 6 + 2 * m*m] = 0.008;//∆P
	Sz[5 + m * 6 + 3 * m*m] = 0.008;//∆Q


	//Elemento 6-6
	Sz[6 + m * 6 + 0 * m*m] = 1;//P
	Sz[6 + m * 6 + 1 * m*m] = 1;//Q
	Sz[6 + m * 6 + 2 * m*m] = 0.008;//∆P
	Sz[6 + m * 6 + 3 * m*m] = 0.008;//∆Q

			//Elemento 6-11
	Sz[11 + m * 6 + 0 * m*m] = 1;//P
	Sz[11 + m * 6 + 1 * m*m] = 1;//Q
	Sz[11 + m * 6 + 2 * m*m] = 0.008;//∆P
	Sz[11 + m * 6 + 3 * m*m] = 0.008;//∆Q

				//Elemento 11-6
	Sz[6 + m * 11 + 0 * m*m] = -1;//P
	Sz[6 + m * 11 + 1 * m*m] = -1;//Q
	Sz[6 + m * 11 + 2 * m*m] = 0.008;//∆P
	Sz[6 + m * 11 + 3 * m*m] = 0.008;//∆Q

			//Elemento 6-12
	Sz[12 + m * 6 + 0 * m*m] = 1;//P
	Sz[12 + m * 6 + 1 * m*m] = 1;//Q
	Sz[12 + m * 6 + 2 * m*m] = 0.008;//∆P
	Sz[12 + m * 6 + 3 * m*m] = 0.008;//∆Q

				//Elemento 12-6
	Sz[6 + m * 12 + 0 * m*m] = -1;//P
	Sz[6 + m * 12 + 1 * m*m] = -1;//Q
	Sz[6 + m * 12 + 2 * m*m] = 0.008;//∆P
	Sz[6 + m * 12 + 3 * m*m] = 0.008;//∆Q

			//Elemento 6-13
	Sz[13 + m * 6 + 0 * m*m] = 1;//P
	Sz[13 + m * 6 + 1 * m*m] = 1;//Q
	Sz[13 + m * 6 + 2 * m*m] = 0.008;//∆P
	Sz[13 + m * 6 + 3 * m*m] = 0.008;//∆Q


	//Elemento 7-7
	Sz[7 + m * 7 + 0 * m*m] = 1;//P
	Sz[7 + m * 7 + 1 * m*m] = 1;//Q
	Sz[7 + m * 7 + 2 * m*m] = 0.008;//∆P
	Sz[7 + m * 7 + 3 * m*m] = 0.008;//∆Q

			//Elemento 7-8
	Sz[8 + m * 7 + 0 * m*m] = 1;//P
	Sz[8 + m * 7 + 1 * m*m] = 1;//Q
	Sz[8 + m * 7 + 2 * m*m] = 0.008;//∆P
	Sz[8 + m * 7 + 3 * m*m] = 0.008;//∆Q

				//Elemento 8-7
	Sz[7 + m * 8 + 0 * m*m] = -1;//P
	Sz[7 + m * 8 + 1 * m*m] = -1;//Q
	Sz[7 + m * 8 + 2 * m*m] = 0.008;//∆P
	Sz[7 + m * 8 + 3 * m*m] = 0.008;//∆Q

			//Elemento 7-9
	Sz[9 + m * 7 + 0 * m*m] = 1;//P
	Sz[9 + m * 7 + 1 * m*m] = 1;//Q
	Sz[9 + m * 7 + 2 * m*m] = 0.008;//∆P
	Sz[9 + m * 7 + 3 * m*m] = 0.008;//∆Q


	//Elemento 8-8
	Sz[8 + m * 8 + 0 * m*m] = 1;//P
	Sz[8 + m * 8 + 1 * m*m] = 1;//Q
	Sz[8 + m * 8 + 2 * m*m] = 0.008;//∆P
	Sz[8 + m * 8 + 3 * m*m] = 0.008;//∆Q


	//Elemento 9-9
	Sz[9 + m * 9 + 0 * m*m] = 1;//P
	Sz[9 + m * 9 + 1 * m*m] = 1;//Q
	Sz[9 + m * 9 + 2 * m*m] = 0.008;//∆P
	Sz[9 + m * 9 + 3 * m*m] = 0.008;//∆Q

			//Elemento 9-10
	Sz[10 + m * 9 + 0 * m*m] = 1;//P
	Sz[10 + m * 9 + 1 * m*m] = 1;//Q
	Sz[10 + m * 9 + 2 * m*m] = 0.008;//∆P
	Sz[10 + m * 9 + 3 * m*m] = 0.008;//∆Q

				//Elemento 10-9
	Sz[9 + m * 10 + 0 * m*m] = -1;//P
	Sz[9 + m * 10 + 1 * m*m] = -1;//Q
	Sz[9 + m * 10 + 2 * m*m] = 0.008;//∆P
	Sz[9 + m * 10 + 3 * m*m] = 0.008;//∆Q

			//Elemento 9-14
	Sz[14 + m * 9 + 0 * m*m] = 1;//P
	Sz[14 + m * 9 + 1 * m*m] = 1;//Q
	Sz[14 + m * 9 + 2 * m*m] = 0.008;//∆P
	Sz[14 + m * 9 + 3 * m*m] = 0.008;//∆Q


	//Elemento 10-10
	Sz[10 + m * 10 + 0 * m*m] = 1;//P
	Sz[10 + m * 10 + 1 * m*m] = 1;//Q
	Sz[10 + m * 10 + 2 * m*m] = 0.008;//∆P
	Sz[10 + m * 10 + 3 * m*m] = 0.008;//∆Q

			//Elemento 10-11
	Sz[11 + m * 10 + 0 * m*m] = 1;//P
	Sz[11 + m * 10 + 1 * m*m] = 1;//Q
	Sz[11 + m * 10 + 2 * m*m] = 0.008;//∆P
	Sz[11 + m * 10 + 3 * m*m] = 0.008;//∆Q

				//Elemento 11-10
	Sz[10 + m * 11 + 0 * m*m] = -1;//P
	Sz[10 + m * 11 + 1 * m*m] = -1;//Q
	Sz[10 + m * 11 + 2 * m*m] = 0.008;//∆P
	Sz[10 + m * 11 + 3 * m*m] = 0.008;//∆Q


	//Elemento 11-11
	Sz[11 + m * 11 + 0 * m*m] = 1;//P
	Sz[11 + m * 11 + 1 * m*m] = 1;//Q
	Sz[11 + m * 11 + 2 * m*m] = 0.008;//∆P
	Sz[11 + m * 11 + 3 * m*m] = 0.008;//∆Q

	//Elemento 12-12
	Sz[12 + m * 12 + 0 * m*m] = 1;//P
	Sz[12 + m * 12 + 1 * m*m] = 1;//Q
	Sz[12 + m * 12 + 2 * m*m] = 0.008;//∆P
	Sz[12 + m * 12 + 3 * m*m] = 0.008;//∆Q

			//Elemento 12-13
	Sz[13 + m * 12 + 0 * m*m] = 1;//P
	Sz[13 + m * 12 + 1 * m*m] = 1;//Q
	Sz[13 + m * 12 + 2 * m*m] = 0.008;//∆P
	Sz[13 + m * 12 + 3 * m*m] = 0.008;//∆Q

				//Elemento 13-12
	Sz[12 + m * 13 + 0 * m*m] = -1;//P
	Sz[12 + m * 13 + 1 * m*m] = -1;//Q
	Sz[12 + m * 13 + 2 * m*m] = 0.008;//∆P
	Sz[12 + m * 13 + 3 * m*m] = 0.008;//∆Q


	//Elemento 13-13
	Sz[13 + m * 13 + 0 * m*m] = 1;//P
	Sz[13 + m * 13 + 1 * m*m] = 1;//Q
	Sz[13 + m * 13 + 2 * m*m] = 0.008;//∆P
	Sz[13 + m * 13 + 3 * m*m] = 0.008;//∆Q

			//Elemento 13-14
	Sz[14 + m * 13 + 0 * m*m] = 1;//P
	Sz[14 + m * 13 + 1 * m*m] = 1;//Q
	Sz[14 + m * 13 + 2 * m*m] = 0.008;//∆P
	Sz[14 + m * 13 + 3 * m*m] = 0.008;//∆Q

				//Elemento 14-13
	Sz[13 + m * 14 + 0 * m*m] = -1;//P
	Sz[13 + m * 14 + 1 * m*m] = -1;//Q
	Sz[13 + m * 14 + 2 * m*m] = 0.008;//∆P
	Sz[13 + m * 14 + 3 * m*m] = 0.008;//∆Q


	//Elemento 14-14
	Sz[14 + m * 14 + 0 * m*m] = 1;//P
	Sz[14 + m * 14 + 1 * m*m] = 1;//Q
	Sz[14 + m * 14 + 2 * m*m] = 0.008;//∆P
	Sz[14 + m * 14 + 3 * m*m] = 0.008;//∆Q


	return;
} //Solicita os valores de potência lido nos medidores da rede

//Função que lê os valores de módulo de tensão nas barras da rede e cria uma matriz de duas colunas, sendo uma com o valor da leitura e outra com a incerteza associada à leitura

void GetMeansVoltageMag14bus350(int m, long double* Vz) {
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
    
	/*
    //Barra 1
    Vz[1+0*m] = 1.060;//V
    Vz[1+1*m] = 0.004;//∆V
    
    //Barra 2
    Vz[2+0*m] = 1.045;//V
    Vz[2+1*m] = 0.004;//∆V
    
    //Barra 3
    Vz[3+0*m] = 1.010;//V
    Vz[3+1*m] = 0.004;//∆V
	*/

    return;
} //Solicita os valores de magnitude das tensões lidos nas barras

