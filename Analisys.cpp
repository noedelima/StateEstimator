//
//  Analisys.cpp
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



int ResidualAnalisys(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S) {
	/******************************************************************
	 *      Cálculo da matriz de ganho G por partes                   *
	 *      Gaa = Haa'*Inv(Raa)*Haa                                   *
	 *      G = Ht*Rn*H, Onde Rn = R^(-1) e Ht = H' (transposta)      *
	 *      Logo, Haat = Haa', RaanHaa = Inv(Raa)*Haa e Gaa = HtRnHaa *
	 ******************************************************************/


	int m = n + 1, option=0;

	int M = DimensionH(m, Sz, Vz);
	int A = DimensionHaa(m, Sz);
	std::cout << "\n\nForam detectadas " << M << " medicoes, sendo " << A << " medicoes de potencia ativa\n";

	long double *Haa, *Raa, *Raan, *I, *sigma, *sigmaaa, *raa, *zaa, *Gaa, *Gaan, *Haat, *HtRnaa, *HaaGaan, *HaaGaanHaat, *Eaa, *rnormaa, *Ynorm, *Gama;
	Haa = new long double[(A + 1)*(n)]; // Matriz de 2 dim
	Raa = new long double[(A + 1)*(A + 1)]; // Matriz de 2 dim
	Raan = new long double[(A + 1)*(A + 1)]; // Matriz de 2 dim
	I = new long double[(A + 1)*(A + 1)]; // Matriz de 2 dim
	sigma = new long double[M + 1]; // Vetor
	sigmaaa = new long double[A + 1];// Vetor
	Gaa = new long double[(n)*(n)]; // Matriz 2 dim
	Gaan = new long double[(n)*(n)]; // Matriz 2 dim
	Haat = new long double[(n)*(A + 1)]; // Matriz 2 dim
	HtRnaa = new long double[(n)*(A + 1)]; // Matriz 2 dim
	HaaGaan = new long double[(A + 1)*(n)]; // Matriz 2 dim
	HaaGaanHaat = new long double[(A + 1)*(A + 1)]; // Matriz 2 dim
	Eaa = new long double[(A + 1)*(A + 1)]; // Matriz 2 dim
	Gama = new long double[(A + 1)*(A + 1)]; // Matriz 2 dim
	zaa = new long double[(3)*(A + 1)]; // Matriz de 2 dim
	Ynorm = new long double[(2)*(m)*(m)]; // Matriz de 3 dim
	rnormaa = new long double[A + 1]; // Vetor
	raa = new long double[A + 1]; // Vetor

	DimensionHaaObs(m, A+1, Sz, zaa);

	//Normaliza Ybus
	for (int i = 1; i < m*m; i++) {
		if (Ybus[i] == 0) {
			Ynorm[i] = 0;
		}
		else if (Ybus[i] < 0) {
			Ynorm[i] = 0;
		}
		else {
			Ynorm[i] = 0;
		}
	} //Condutância
	for (int i = 1 + m*m; i < 2 * m*m; i++) {
		if (Ybus[i]==0) {
			Ynorm[i] = 0;
		}
		else if (Ybus[i] < 0) {
			Ynorm[i] = -1;
		} 
		else {
			Ynorm[i] = 1;
		}
	} //Susceptância


	Residual(m, Sz, Vz, sigma);
	for (int i = 1; i <= A; i++) { sigmaaa[i] = (sigma[i] * sigma[i]); }
	Diag(A, sigmaaa, Raa); //Cria a matriz de resíduos Raa em Raan
	for (int i = 1; i <= A; i++) { sigmaaa[i] = 1 / sigma[i]; }
	Diag(A, sigmaaa, Raan); //Cria a matriz de resíduos Raa^-1 em Raan

	for (int i = 1; i <= A; i++) { sigmaaa[i] = 1; }
	Diag(A, sigmaaa, I); //Cria a matriz de resíduos Raa em Raan


	//Matriz Gaa^-1
	DecoupledJacobianObs(m, Ynorm, V, Theta, Sz, Haa);
	Transp(Haa, A + 1, n, Haat);
	Product(Haat, n, A + 1, Haa, A + 1, n, Gaa);
	std::cout << "\nDeseja analisar a observabilidade da rede?\nObserve que para sistemas com muitas barras o programa pode travar!\nDigite 1 para sim, ou qualquer valor para prosseguir sem verificar: ";
	std::cin >> option;
	if (option == 1) { if (Det(n, Gaa) == 0) { std::cout << "\nDet(Gaa) = 0\n"; return 1; } }
	std::cout << "\n\n";
	InvertMatrix(n, Gaa, Gaan); //Calcula Inv(Gaa)

	//Calcular E = I - H*(G^-1)*Ht
	Product(Haa, A + 1, n, Gaan, n, n, HaaGaan);
	Product(HaaGaan, A + 1, n, Haat, n, A + 1, HaaGaanHaat);
	for (int i = 0; i < (A + 1); i++) {
		for (int j = 0; j < (A + 1); j++) {
			Eaa[j + (A + 1)*i] = I[j + (A + 1)*i] - HaaGaanHaat[j + (A + 1)*i];
		}
	} //Calcula a matriz de Covariância dos Resíduos E em Eaa
	for (int i = 0; i < (A + 1); i++) {
		for (int j = 0; j < (A + 1); j++) {
			Gama[j + (A + 1)*i] = sqrt(pow(Eaa[j + (A + 1)*i], 2)) / sqrt(Eaa[i + (A + 1)*i] * Eaa[j + (A + 1)*j]);
		}
	} //Calcula a matriz de Coeficientes de Correlação Gama

	for (int i = 1; i < A + 1; i++) {
		if (zaa[i + (A + 1) * 1] > zaa[i + (A + 1) * 2]) { sigma[i] = -1; }
		else { sigma[i] = 1; }
	}
	//Calcula o vetor de resíduos e resíduos normalizados rnorm[i] = |r[i]|/sqrt(E[i][i])
	Product(Eaa, A+1, A+1, sigma, A+1, 1, raa); // Utiliza a fórmula r = E*z, onde z é um vetor de elementos unitários. Aqui foi utilizado sigma, pois armazena os valores unitários utilizados para gerar I
	for (int i = 1; i < A + 1; i++) {rnormaa[i] = sqrt(pow(raa[i], 2)) / sqrt(Eaa[i + (A + 1)*i]);} // Calcula os resíduos normalizados																									// Apresentação dos parâmetros de resíduo

	/*
	// Apresentação dos parâmetros de resíduos
	std::cout << "\nMatriz de Covariancia dos Residuos (E):\n\n";
	for (int i = 1; i < A + 1; i++) {
		std::cout << "\nLinha " << i << ":";
		for (int j = 1; j < A + 1; j++) {
			std::cout << "			" << Eaa[j + (A + 1)*i];
		}
	} //Apresenta a matriz de covariância
	std::cout << "\n\n";
	std::cout << "\nMatriz de Coeficientes de Correlacao (Gama):\n\n";
	for (int i = 1; i < A + 1; i++) {
		std::cout << "\nLinha " << i << ":";
		for (int j = 1; j < A + 1; j++) {
			std::cout << "			" << Gama[j + (A + 1)*i];
		}
	} //Apresenta a matriz de correlação
	std::cout << "\n\n";
	std::cout << "\nVetor de Residuos:\n";
	for (int i = 1; i < A + 1; i++) {
		std::cout << "P(" << zaa[i + (A + 1) * 1] << "-" << zaa[i + (A + 1) * 2] << ") ===> " << raa[i] << "\n";
	}
	std::cout << "\nVetor de Residuos Normalizados:\n";
	for (int i = 1; i < A + 1; i++) {
		std::cout << "P(" << zaa[i + (A + 1) * 1] << "-" << zaa[i + (A + 1) * 2] << ") ===> " << rnormaa[i] << "\n";
	}
	// Fim da apresentação de resíduos */


	//Detecção de medidas críticas e conjuntos críticos
	for (int i = 1; i < A + 1; i++) {
		if (pow(raa[i], 2) < pow(10, -6)) {
			if (zaa[i + (A + 1) * 1] != zaa[i + (A + 1) * 2]) {
				std::cout << "O medidor do fluxo de potencia da barra " << zaa[i + (A + 1) * 1] << " para a barra " << zaa[i + (A + 1) * 2] << " e critico!";
			}
			else {
				std::cout << "\nO medidor de potencia injetada na barra " << zaa[i + (A + 1) * 1] << " e critico!\n";
			}
		}
	} // Apresenta lista de Cmed - medidas críticas

	//Detecção de pares críticos
	for (int i = 1; i < A + 1; i++) {
		for (int j = 1; j < A + 1; j++) {
			if (j > i) {
				if (pow(rnormaa[i] - rnormaa[j], 2) < 0.00001 && Gama[j + (A + 1)*i] > 0.99) { // incluído tolerância de e-10 para compensar erros de cálculo na comparação
					std::cout << "\nO medidor " << zaa[i + (A + 1) * 1] << "-" << zaa[i + (A + 1) * 2] << " forma uma dupla critica com o medidor " << zaa[j + (A + 1) * 1] << "-" << zaa[j + (A + 1) * 2] << "\n";
				}
			}
		}
	} // Apresenta lista de duplas de Cconj - conjuntos críticos

	//Detecção de trios críticos
	for (int i = 1; i < A + 1; i++) {
		for (int j = 1; j < A + 1; j++) {
			for (int k = 1; k < A + 1; k++) {
				if (k > j && j > i) {
					if (pow(rnormaa[i] - rnormaa[j], 2) < 0.00001 && pow(rnormaa[i] - rnormaa[k], 2) < 0.00001 && (Gama[j + (A + 1)*i] > 0.99) || Gama[k + (A + 1)*i] > 0.99 || Gama[k + (A + 1)*j] > 0.99) { // incluído tolerância de e-10 para compensar erros de cálculo na comparação
						std::cout << "\nO medidor " << zaa[i + (A + 1) * 1] << "-" << zaa[i + (A + 1) * 2] << " forma um trio de medidas criticas com os medidores " << zaa[j + (A + 1) * 1] << "-" << zaa[j + (A + 1) * 2] << " e " << zaa[k + (A + 1) * 1] << "-" << zaa[k + (A + 1) * 2] << "\n";
					}
				}
			}
		}
	} // Apresenta a lista de trios de Cconj - conjuntos críticos

	// Grau de redundância das medidas do sistema
	// Utilizaremos a fórmula n = (m-n)/(max-n), onde max = (B+2L)/(B-1), e B=barras, L=Linhas
	float red = 0;
	for (int i = 1; i < m; i++) {
		for (int j = 1; j < m; j++) {
			if (i != j && Ynorm[j + m*i + m*m * 1] != 0) { ++red; }
		}
	} //Calcula linha = 2L
	red += n;
	red = (2 * n - 1) * red / (n - 1);
	red = 1 + (2 * A - 2 * n + 1) / (red - 2 * n + 1);
	std::cout << "\n\nA rede apresenta nivel de redundancia igual a " << 100*red << "%, sendo a redundancia global = " << 100*A/n << "%\n";

	delete[] Haa; // Matriz de 2 dim
	delete[] Raa; // Matriz de 2 dim
	delete[] Raan; // Matriz de 2 dim
	delete[] I; // Matriz de 2 dim
	delete[] sigma; // Vetor
	delete[] sigmaaa;// Vetor
	delete[] Gaa; // Matriz 2 dim
	delete[] Gaan; // Matriz 2 dim
	delete[] Haat; // Matriz 2 dim
	delete[] HtRnaa; // Matriz 2 dim
	delete[] HaaGaan; // Matriz 2 dim
	delete[] HaaGaanHaat; // Matriz 2 dim
	delete[] Eaa; // Matriz 2 dim
	delete[] Gama; // Matriz 2 dim
	delete[] zaa; // Matriz de 2 dim
	delete[] Ynorm; // Matriz de 3 dim
	delete[] rnormaa; // Vetor
	delete[] raa; // Vetor

	return 0;
}
