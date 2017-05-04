//
//  AuxFunctions.hpp
//  Estimador de Estado
//
//  Created by Noé de Lima Bezerra on 31/12/00.
//  Copyright © 2000 Noé de Lima Bezerra. All rights reserved.
//

#ifndef AuxFunctions_hpp
#define AuxFunctions_hpp

#include <stdio.h>

//Funções de entrada de dados
void GetYbus(int n, long double* Ybus);
void GetMeansPower(int m, long double* Sz);
void GetMeansVoltageMag(int m, long double* Vz);

//Função de estimação do estado
void CalcPower(int m, long double* Ybus, long double* S, long double* V, long double* Theta);

//Funções para calcular o Jacobiano
int DimensionH(int m, long double* Sz, long double* Vz);
int DimensionHaa(int m, long double* Sz);
int DimensionHrr(int m, long double* Sz, long double* Vz);
int DimensionHaaObs(int m, int layer, long double* Sz, long double* zaa);
int Jacobian(int m, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* H);
void Residual(int m, long double* Sz, long double* Vz, long double* sigma);
long double Max(int m, long double* Sz, long double* S, long double* Vz, long double* V, long double* zh);
int DecoupledJacobianAA(int m, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Haa);
int DecoupledJacobianRR(int m, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* Hrr);
int DecoupledJacobianObs(int m, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Haa);

//Funções para solução do sistema
int NewtonRaphson(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S, long double e);
int FastNewtonRaphson(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S, long double e);
int FastDecoupled(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S, long double e);
int Decoupled(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S, long double e);

//Funções matriciais
void Transp(long double* A, int line, int row, long double* T);
void Product(long double* A, int aline, int arow, long double* B, int bline, int brow, long double* P);
int InvertMatrix(int n, long double* X, long double* Y);
void Diag(int M, long double* sigma, long double* R);
long double Det(int m, long double* A);

//Funções de entrada dos exemplos
//Exemplo 2.6 do livro do Ali Abur
void GetYbus1(int n, long double* Ybus);
void GetMeansPower1(int m, long double* Sz);
void GetMeansVoltageMag1(int m, long double* Vz);
//IEEE 14-bus 150%
void GetYbus14bus150(int n, long double* Ybus);
void GetMeansPower14bus150(int m, long double* Sz);
void GetMeansVoltageMag14bus150(int m, long double* Vz);
//IEEE 14-bus 250%
void GetYbus14bus250(int n, long double* Ybus);
void GetMeansPower14bus250(int m, long double* Sz);
void GetMeansVoltageMag14bus250(int m, long double* Vz);
//IEEE 14-bus 350%
void GetYbus14bus350(int n, long double* Ybus);
void GetMeansPower14bus350(int m, long double* Sz);
void GetMeansVoltageMag14bus350(int m, long double* Vz);
//Utilizar para teste

//Funções de Análise de Resíduos e Observabilidade
int ResidualAnalisys(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S);

#endif /* AuxFunctions_hpp */
