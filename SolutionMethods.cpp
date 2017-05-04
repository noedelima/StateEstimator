//
//  SolutionMethods.cpp
//  Estimador de Estado
//
//  Created by Noé de Lima Bezerra on 14/11/16.
//  Copyright © 2016 Noé de Lima Bezerra. All rights reserved.
//

#include "AuxFunctions.hpp"
#include <iostream>
#include <cmath>
#include <complex>
#include <new>
using namespace std;

int NewtonRaphson(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S, long double e) {
    /******************************************************************
     *      Cálculo da matriz de ganho G                              *
     *      G = H'*Inv(R)*H                                           *
     *      Aqui, R já está invertida, portanto, temos:               *
     *      G = Ht*Rn*H, Onde Rn = R^(-1) e Ht = H' (transposta)      *
     *      Logo, Ht = H', RnH = Inv(R)*H e G = HtRnH                 *
     ******************************************************************/
    

    int m=n+1, countit=0, conv=0;
    int M = DimensionH(m, Sz, Vz);
    //std::cout << "\n\nForam detectadas " << M << " medicoes\n";

    long double *H, *Rn, *sigma;
	H = new long double [(M + 1)*(2 * n + 1)]; // Matriz 2 dim
	Rn = new long double [(M + 1)*(M + 1)]; // Matriz 2 dim
	sigma = new long double [M + 1]; // Matriz 2 dim

    Residual(m, Sz, Vz, sigma);
    for (int i=1; i<=M; i++) {sigma[i] = 1/(sigma[i]*sigma[i]);} //Transforma o vetor sigma em sigma.^(-2)
    Diag(M, sigma, Rn); //Cria a matriz de resíduos R^-1 em Rn
    
    
    long double *G, *GG, *Gn, *Ht, *HtRn, max=10, *zh, *t, *dx;
	G = new long double [(2 * n + 1)*(2 * n + 1)]; // Matriz 2 dim
	GG = new long double [(2 * n)*(2 * n)]; // Matriz 2 dim
	Gn = new long double [(2 * n)*(2 * n)]; // Matriz 2 dim
	Ht = new long double [(2 * n + 1)*(M + 1)]; // Matriz 2 dim
	HtRn = new long double [(2 * n + 1)*(M + 1)]; // Matriz 2 dim
	zh = new long double [M + 1]; // Vetor
	t = new long double [2 * n + 1]; // Vetor
	dx = new long double [2 * n + 1]; // Vetor
    
    //Utilizando a solução por Newton-Raphson
    while (max > e) {
        countit = countit+1; //Conta a iteração

        //Calcula a matriz inversa dentro da iteração
        Jacobian(m, Ybus, V, Theta, Sz, Vz, H);
        Transp(H, M+1, 2*n+1, Ht);
        Product(Ht, 2*n+1, M+1, Rn, M+1, M+1, HtRn);
        Product(HtRn, 2*n+1, M+1, H, M+1, 2*n+1, G);

        for (int i=1; i<2*n; i++) {
            for (int j=1; j<2*n; j++) {
                GG[j+(2*n)*i] = G[(j+1)+(2*n+1)*(i+1)];
            }
        } //Redimensiona a matriz G na matriz GG para retirar os resíduos de memória antes de inverter

		conv = InvertMatrix(2*n, GG, Gn); //Calcula Inv(G)
		if (conv!=0) {return 2;}

        CalcPower(m, Ybus, S, V, Theta); //Estima o estado do sistema
        Max(m, Sz, S, Vz, V, zh); //Calcula (z-h(X))
        Product(HtRn, 2*n+1, M+1, zh, M+1, 1, t); //Calcula t = H'Inv(R)(z-h(x))
        for (int i=1; i<2*n; i++) {t[i]=t[i+1];} //Desloca t para baixo para excluir a barra 1
        Product(Gn, 2*n, 2*n, t, 2*n, 1, dx+1); //Calcula dx = Inv(G)*t ***Note que o resultado é colocado a partir do segundo endereço de dx (utilizando dx+1), para não impactar nas referências das barras. O deslocamento retirando a barra 1 visa corrigir os erros da inversão da matriz G.
        for (int i=2; i<m; i++) {
            Theta[i] = Theta[i]+dx[i];
        } //Calcula o novo estado do sistema x = x+dx para os ângulos
        for (int i=1; i<m; i++) {
            V[i] = V[i]+dx[i+n];
        } //Calcula o novo estado do sistema x = x+dx para as tenções
        for (int i=2; i<=2*n; i++) {
            max=0;
            long double comp = sqrt(pow(dx[i], 2));
            if (comp > max) {max = comp;}
        }
        //cout << "Iteração " << countit << ", tolerância " << e << " e |∆x|máx = " << max << "\n";
        if (countit>=100) {std::cout << "\nO programa nao convergiu em 100 iteracoes\n"; return 1;}
    }
    std::cout << "\nA solucao convergiu apos " << countit << " iteracoes\n";

	delete[] H; // deleta Matriz 2 dim da memória
	delete[] Rn; // deleta Matriz 2 dim da memória
	delete[] sigma; // deleta Matriz 2 dim da memória
	delete[] G; // deleta Matriz 2 dim da memória
	delete[] GG; // deleta Matriz 2 dim da memória
	delete[] Gn; // deleta Matriz 2 dim da memória
	delete[] Ht; // deleta Matriz 2 dim da memória
	delete[] HtRn; // deleta Matriz 2 dim da memória
	delete[] zh; // deleta Vetor da memória
	delete[] t; // deleta Vetor da memória
	delete[] dx; // deleta Vetor da memória

    return 0;
}


int FastNewtonRaphson(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S, long double e) {
    /******************************************************************
     *      Cálculo da matriz de ganho G                              *
     *      G = H'*Inv(R)*H                                           *
     *      Aqui, R já está invertida, portanto, temos:               *
     *      G = Ht*Rn*H, Onde Rn = R^(-1) e Ht = H' (transposta)      *
     *      Logo, Ht = H', RnH = Inv(R)*H e G = HtRnH                 *
     *      Neste método a matriz G será invertida apenas uma vez!    *
     ******************************************************************/
    
    
    int m=n+1, countit=0, conv=0;
    int M = DimensionH(m, Sz, Vz);
    //std::cout << "\n\nForam detectadas " << M << " medicoes\n";

    long double *H, *Rn, *sigma;
	H = new long double [(M + 1)*(2 * n + 1)]; // Matriz 2 dim
	Rn = new long double [(M + 1)*(M + 1)]; // Matriz 2 dim
	sigma = new long double [M + 1]; // Vetor

    Residual(m, Sz, Vz, sigma);
    for (int i=1; i<=M; i++) {sigma[i] = 1/(sigma[i]*sigma[i]);} //Transforma o vetor sigma em sigma.^(-2)
    Diag(M, sigma, Rn); //Cria a matriz de resíduos R^-1 em Rn
    
    
    long double *G, *GG, *Gn, *Ht, *HtRn, max=10, *zh, *t, *dx;
	G = new long double [(2 * n + 1)*(2 * n + 1)]; // Matriz 2 dim
	GG = new long double [(2 * n)*(2 * n)]; // Matriz 2 dim
	Gn = new long double [(2 * n)*(2 * n)]; // Matriz 2 dim
	Ht = new long double [(2 * n + 1)*(M + 1)]; // Matriz 2 dim
	HtRn = new long double [(2 * n + 1)*(M + 1)]; // Matriz 2 dim
	zh = new long double [M + 1]; // Vetor
	t = new long double [2 * n + 1]; // Vetor
	dx = new long double [2 * n + 1]; // Vetor

    //Calcula a matriz inversa fora da iteração
    Jacobian(m, Ybus, V, Theta, Sz, Vz, H);
    Transp(H, M+1, 2*n+1, Ht);
    Product(Ht, 2*n+1, M+1, Rn, M+1, M+1, HtRn);
    Product(HtRn, 2*n+1, M+1, H, M+1, 2*n+1, G);
    for (int i=1; i<2*n; i++) {
        for (int j=0;j<2*n; j++) {
            GG[j+(2*n)*i] = G[(j+1)+(2*n+1)*(i+1)];
        }
    } //Redimensiona a matriz G na matriz GG para retirar os resíduos de memória antes de inverter
    conv = InvertMatrix(2*n, GG, Gn); //Calcula Inv(G)
    if (conv!=0) {return 2;}
    
    //Utilizando a solução por Newton-Raphson
    while (max > e) {
        countit = countit+1; //Conta a iteração
        
        CalcPower(m, Ybus, S, V, Theta); //Estima o estado do sistema
        Max(m, Sz, S, Vz, V, zh); //Calcula (z-h(X))
        Product(HtRn, 2*n+1, M+1, zh, M+1, 1, t); //Calcula t = H'Inv(R)(z-h(x))
        for (int i=1; i<2*n; i++) {t[i]=t[i+1];} //Desloca t para baixo para excluir a barra 1
        Product(Gn, 2*n, 2*n, t, 2*n, 1, dx+1); //Calcula dx = Inv(G)*t ***Note que o resultado é colocado a partir do segundo endereço de dx (utilizando dx+1), para não impactar nas referências das barras. O deslocamento retirando a barra 1 visa corrigir os erros da inversão da matriz G.
        for (int i=2; i<m; i++) {
            Theta[i] = Theta[i]+dx[i];
        } //Calcula o novo estado do sistema x = x+dx para os ângulos
        for (int i=1; i<m; i++) {
            V[i] = V[i]+dx[i+n];
        } //Calcula o novo estado do sistema x = x+dx para as tenções
        for (int i=2; i<=2*n; i++) {
            max=0;
            long double comp = sqrt(pow(dx[i], 2));
            if (comp > max) {max = comp;}
        }
        //cout << "Iteração " << countit << ", tolerância " << e << " e |∆x|máx = " << max << "\n";
        if (countit>=100) {std::cout << "\nO programa nao convergiu em 100 iteracoes\n"; return 1;}
    }
    std::cout << "\nA solucao convergiu apos " << countit << " iteracoes\n";

	delete[] H; // deleta Matriz 2 dim da memória
	delete[] Rn; // deleta Matriz 2 dim da memória
	delete[] sigma; // deleta Matriz 2 dim da memória
	delete[] G; // deleta Matriz 2 dim da memória
	delete[] GG; // deleta Matriz 2 dim da memória
	delete[] Gn; // deleta Matriz 2 dim da memória
	delete[] Ht; // deleta Matriz 2 dim da memória
	delete[] HtRn; // deleta Matriz 2 dim da memória
	delete[] zh; // deleta Vetor da memória
	delete[] t; // deleta Vetor da memória
	delete[] dx; // deleta Vetor da memória

    return 0;
}


int Decoupled(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S, long double e) {
    /******************************************************************
     *      Cálculo da matriz de ganho G por partes                   *
     *      Gaa = Haa'*Inv(Raa)*Haa                                   *
     *      Grr = Hrr'*Inv(Rrr)*Hrr                                   *
     *      G = Ht*Rn*H, Onde Rn = R^(-1) e Ht = H' (transposta)      *
     *      Logo, Haat = Haa', RaanHaa = Inv(Raa)*Haa e Gaa = HtRnHaa *
     *      Neste método a matriz Gaa será invertida a cada iteração  *
     *      Neste método a matriz Grr será invertida a cada iteração  *
     ******************************************************************/
    
    
    int m=n+1, conv=0;
    double countit=0;
    
    int M = DimensionH(m, Sz, Vz);
    //std::cout << "\n\nForam detectadas " << M << " medicoes\n";
    int A = DimensionHaa(m, Sz);
    int R = DimensionHrr(m, Sz, Vz);
    
	long double *Haa, *Hrr, *Raan, *Rrrn, *sigma, *sigmaaa, *sigmarr, *Gaa, *Gaan, *Grr, *Grrn, *Haat, *Hrrt, *HtRnaa, *HtRnrr, max = 10, *zh, *zhaa, *zhrr, *taa, *trr, *dtheta, *dV;
	Haa = new long double[(A + 1)*(n)]; // Matriz de 2 dim
	Hrr = new long double[(R + 1)*(m)]; // Matriz de 2 dim
	Raan = new long double[(A + 1)*(A + 1)]; // Matriz de 2 dim
	Rrrn = new long double[(R + 1)*(R + 1)]; // Matriz de 2 dim
	sigma = new long double[M + 1]; // Vetor
	sigmaaa = new long double[A + 1];// Vetor
	sigmarr = new long double[R + 1]; // Vetor
	Gaa = new long double[(n)*(n)]; // Matriz 2 dim
	Gaan = new long double[(n)*(n)]; // Matriz 2 dim
	Grr = new long double[(m)*(m)]; // Matriz 2 dim
	Grrn = new long double[(m)*(m)]; // Matriz 2 dim
	Haat = new long double[(n)*(A + 1)]; // Matriz 2 dim
	Hrrt = new long double[(m)*(R + 1)]; // Matriz 2 dim
	HtRnaa = new long double[(n)*(A + 1)]; // Matriz 2 dim
	HtRnrr = new long double[(m)*(R + 1)]; // Matriz 2 dim
	zh = new long double[M + 1]; // Vetor
	zhaa = new long double[A + 1]; // Vetor
	zhrr = new long double[R + 1]; // Vetor
	taa = new long double[n]; // Vetor
	trr = new long double[m]; // Vetor
	dtheta = new long double[m]; // Vetor
	dV = new long double[m]; // Vetor

    Residual(m, Sz, Vz, sigma);
    for (int i=1; i<=A; i++) {sigmaaa[i] = 1/(sigma[i]*sigma[i]);}
    for (int i=1; i<=R; i++) {sigmarr[i] = 1/(sigma[i+A]*sigma[i+A]);}
    
    Diag(A, sigmaaa, Raan); //Cria a matriz de resíduos Raa^-1 em Raan
    Diag(R, sigmarr, Rrrn); //Cria a matriz de resíduos Rrr^-1 em Rrrn
    
    //Utilizando a solução pelo método desacoplado
    while (max > e) {
        countit = countit+0.5; //Conta meia iteração
        
        //Calcula as matrizes inversas fora da iteração
        //Matriz Gaa^-1
        DecoupledJacobianAA(m, Ybus, V, Theta, Sz, Haa);
        Transp(Haa, A+1, n, Haat);
        Product(Haat, n, A+1, Raan, A+1, A+1, HtRnaa);
        Product(HtRnaa, n, A+1, Haa, A+1, n, Gaa);
        conv = InvertMatrix(n, Gaa, Gaan); //Calcula Inv(Gaa)
        if (conv!=0) {return 2;}
        
        //Matriz Grr^-1
        DecoupledJacobianRR(m, Ybus, V, Theta, Sz, Vz, Hrr);
        Transp(Hrr, R+1, m, Hrrt);
        Product(Hrrt, m, R+1, Rrrn, R+1, R+1, HtRnrr);
        Product(HtRnrr, m, R+1, Hrr, R+1, m, Grr);
        conv = InvertMatrix(m, Grr, Grrn); //Calcula Inv(Grr)
        if (conv!=0) {return 3;}
        
        CalcPower(m, Ybus, S, V, Theta); //Estima o estado do sistema
        Max(m, Sz, S, Vz, V, zh); //Calcula (z-h(X))
        for (int i=1; i<=A; i++) {zhaa[i]=zh[i];}
        for (int i=1; i<=R; i++) {zhrr[i]=zh[i+A];}
        
        
        Product(HtRnaa, n, A+1, zhaa, A+1, 1, taa); //Calcula taa = Haa'Inv(Raa)(z-h(x)aa)
        Product(Gaan, n, n, taa, n, 1, dtheta+1); //Calcula dtheta = Inv(Gaa)*taa ***Note que o resultado é colocado a partir do segundo endereço de dtheta (utilizando dtheta+1), para não impactar nas referências das barras. O deslocamento retirando a barra 1 visa corrigir os erros de referência nos ângulos
        for (int i=2; i<m; i++) {
            Theta[i] += dtheta[i];
        } //Calcula o novo estado do sistema x = x+dx pata o ângulo de fase nas barras
        
        
        max=0;
        for (int i=2; i<m; i++) {
            long double comp = sqrt(pow(dtheta[i], 2));
            if (comp > max) {max = comp;}
        }
        for (int i=1; i<m; i++) {
            long double comp = sqrt(pow(dV[i], 2));
            if (comp > max) {max = comp;}
        }
        //cout << "Iteração " << countit << ", tolerância " << e << " e |∆x|máx = " << max << "\n";
        
        if (max>e) {
            countit = countit+0.5; //Conta meia iteração
            
            Product(HtRnrr, m, R+1, zhrr, R+1, 1, trr); //Calcula trr = Hrr'Inv(Rrr)(z-h(x)rr)
            Product(Grrn, m, m, trr, m, 1, dV); //Calcula dV = Inv(Grr)*trr
            for (int i=1; i<m; i++) {
                V[i] += dV[i];
            } //Calcula o novo estado do sistema x = x+dx pata o ângulo de fase nas barras
            
            
            max=0;
            for (int i=2; i<m; i++) {
                long double comp = sqrt(pow(dtheta[i], 2));
                if (comp > max) {max = comp;}
            }
            for (int i=1; i<m; i++) {
                long double comp = sqrt(pow(dV[i], 2));
                if (comp > max) {max = comp;}
            }
            //cout << "Iteração " << countit << ", tolerância " << e << " e |∆x|máx = " << max << "\n";
        }
        if (countit>=100) {std::cout << "\nO programa nao convergiu em 100 iteracoes\n"; return 1;}
    }
    std::cout << "\nA solucao convergiu apos " << countit << " iteracoes\n";

	delete[] Haa; // deleta Matriz de 2 dim da memória
	delete[] Hrr; // deleta Matriz de 2 dim da memória
	delete[] Raan; // deleta Matriz de 2 dim da memória
	delete[] Rrrn; // deleta Matriz de 2 dim da memória
	delete[] sigma; // deleta Vetor da memória
	delete[] sigmaaa;// deleta Vetor da memória
	delete[] sigmarr; // deleta Vetor da memória
	delete[] Gaa; // deleta Matriz 2 dim da memória
	delete[] Gaan; // deleta Matriz 2 dim da memória
	delete[] Grr; // deleta Matriz 2 dim da memória
	delete[] Grrn; // deleta Matriz 2 dim da memória
	delete[] Haat; // deleta Matriz 2 dim da memória
	delete[] Hrrt; // deleta Matriz 2 dim da memória
	delete[] HtRnaa; // deleta Matriz 2 dim da memória
	delete[] HtRnrr; // deleta Matriz 2 dim da memória
	delete[] zh; // deleta Vetor da memória
	delete[] zhaa; // deleta Vetor da memória
	delete[] zhrr; // deleta Vetor da memória
	delete[] taa; // deleta Vetor da memória
	delete[] trr; // deleta Vetor da memória
	delete[] dtheta; // deleta Vetor da memória
	delete[] dV; // deleta Vetor da memória
	
    return 0;
}


int FastDecoupled(int n, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* S, long double e) {
    /******************************************************************
     *      Cálculo da matriz de ganho G por partes                   *
     *      Gaa = Haa'*Inv(Raa)*Haa                                   *
     *      Grr = Hrr'*Inv(Rrr)*Hrr                                   *
     *      G = Ht*Rn*H, Onde Rn = R^(-1) e Ht = H' (transposta)      *
     *      Logo, Haat = Haa', RaanHaa = Inv(Raa)*Haa e Gaa = HtRnHaa *
     *      Neste método a matriz Gaa será invertida apenas uma vez!  *
     *      Neste método a matriz Grr será invertida apenas uma vez!  *
     ******************************************************************/
    
    
	int m = n + 1, conv = 0;
	double countit = 0;

	int M = DimensionH(m, Sz, Vz);
	//std::cout << "\n\nForam detectadas " << M << " medicoes\n";
	int A = DimensionHaa(m, Sz);
	int R = DimensionHrr(m, Sz, Vz);

	long double *Haa, *Hrr, *Raan, *Rrrn, *sigma, *sigmaaa, *sigmarr, *Gaa, *Gaan, *Grr, *Grrn, *Haat, *Hrrt, *HtRnaa, *HtRnrr, max = 10, *zh, *zhaa, *zhrr, *taa, *trr, *dtheta, *dV;
	Haa = new long double[(A + 1)*(n)]; // Matriz de 2 dim
	Hrr = new long double[(R + 1)*(m)]; // Matriz de 2 dim
	Raan = new long double[(A + 1)*(A + 1)]; // Matriz de 2 dim
	Rrrn = new long double[(R + 1)*(R + 1)]; // Matriz de 2 dim
	sigma = new long double[M + 1]; // Vetor
	sigmaaa = new long double[A + 1];// Vetor
	sigmarr = new long double[R + 1]; // Vetor
	Gaa = new long double[(n)*(n)]; // Matriz 2 dim
	Gaan = new long double[(n)*(n)]; // Matriz 2 dim
	Grr = new long double[(m)*(m)]; // Matriz 2 dim
	Grrn = new long double[(m)*(m)]; // Matriz 2 dim
	Haat = new long double[(n)*(A + 1)]; // Matriz 2 dim
	Hrrt = new long double[(m)*(R + 1)]; // Matriz 2 dim
	HtRnaa = new long double[(n)*(A + 1)]; // Matriz 2 dim
	HtRnrr = new long double[(m)*(R + 1)]; // Matriz 2 dim
	zh = new long double[M + 1]; // Vetor
	zhaa = new long double[A + 1]; // Vetor
	zhrr = new long double[R + 1]; // Vetor
	taa = new long double[n]; // Vetor
	trr = new long double[m]; // Vetor
	dtheta = new long double[m]; // Vetor
	dV = new long double[m]; // Vetor

	Residual(m, Sz, Vz, sigma);
	for (int i = 1; i <= A; i++) { sigmaaa[i] = 1 / (sigma[i] * sigma[i]); }
	for (int i = 1; i <= R; i++) { sigmarr[i] = 1 / (sigma[i + A] * sigma[i + A]); }

	Diag(A, sigmaaa, Raan); //Cria a matriz de resíduos Raa^-1 em Raan
	Diag(R, sigmarr, Rrrn); //Cria a matriz de resíduos Rrr^-1 em Rrrn

        
    //Calcula as matrizes inversas fora da iteração
    //Matriz Gaa^-1
    DecoupledJacobianAA(m, Ybus, V, Theta, Sz, Haa);
    Transp(Haa, A+1, n, Haat);
    Product(Haat, n, A+1, Raan, A+1, A+1, HtRnaa);
    Product(HtRnaa, n, A+1, Haa, A+1, n, Gaa);
    conv = InvertMatrix(n, Gaa, Gaan); //Calcula Inv(Gaa)
    if (conv!=0) {return 2;}
    
    //Matriz Grr^-1
    DecoupledJacobianRR(m, Ybus, V, Theta, Sz, Vz, Hrr);
    Transp(Hrr, R+1, m, Hrrt);
    Product(Hrrt, m, R+1, Rrrn, R+1, R+1, HtRnrr);
    Product(HtRnrr, m, R+1, Hrr, R+1, m, Grr);
    conv = InvertMatrix(m, Grr, Grrn); //Calcula Inv(Grr)
    if (conv!=0) {return 3;}
    
	//Utilizando a solução pelo método desacoplado
    while (max > e) {
        countit = countit+0.5; //Conta meia iteração
        
        CalcPower(m, Ybus, S, V, Theta); //Estima o estado do sistema
        Max(m, Sz, S, Vz, V, zh); //Calcula (z-h(X))
        for (int i=1; i<=A; i++) {zhaa[i]=zh[i];}
        for (int i=1; i<=R; i++) {zhrr[i]=zh[i+A];}
        
        
        Product(HtRnaa, n, A+1, zhaa, A+1, 1, taa); //Calcula taa = Haa'Inv(Raa)(z-h(x)aa)
        Product(Gaan, n, n, taa, n, 1, dtheta+1); //Calcula dtheta = Inv(Gaa)*taa ***Note que o resultado é colocado a partir do segundo endereço de dtheta (utilizando dtheta+1), para não impactar nas referências das barras. O deslocamento retirando a barra 1 visa corrigir os erros de referência nos ângulos
        for (int i=2; i<m; i++) {
            Theta[i] += dtheta[i];
        } //Calcula o novo estado do sistema x = x+dx pata o ângulo de fase nas barras
        
        
        max=0;
        for (int i=2; i<m; i++) {
            long double comp = sqrt(pow(dtheta[i], 2));
            if (comp > max) {max = comp;}
        }
        for (int i=1; i<m; i++) {
            long double comp = sqrt(pow(dV[i], 2));
            if (comp > max) {max = comp;}
        }
        //cout << "Iteração " << countit << ", tolerância " << e << " e |∆x|máx = " << max << "\n";
        
        if (max>e) {
            countit = countit+0.5; //Conta meia iteração
            
            Product(HtRnrr, m, R+1, zhrr, R+1, 1, trr); //Calcula trr = Hrr'Inv(Rrr)(z-h(x)rr)
            Product(Grrn, m, m, trr, m, 1, dV); //Calcula dV = Inv(Grr)*trr
            for (int i=1; i<m; i++) {
                V[i] += dV[i];
            } //Calcula o novo estado do sistema x = x+dx pata o ângulo de fase nas barras
            
            
            max=0;
            for (int i=2; i<m; i++) {
                long double comp = sqrt(pow(dtheta[i], 2));
                if (comp > max) {max = comp;}
            }
            for (int i=1; i<m; i++) {
                long double comp = sqrt(pow(dV[i], 2));
                if (comp > max) {max = comp;}
            }
            //cout << "Iteração " << countit << ", tolerância " << e << " e |∆x|máx = " << max << "\n";
        }
        if (countit>=100) {std::cout << "\nO programa nao convergiu em 100 iteracoes\n"; return 1;}
    }
    std::cout << "\nA solucao convergiu apos " << countit << " iteracoes\n";
	
	delete[] Haa; // deleta Matriz de 2 dim da memória
	delete[] Hrr; // deleta Matriz de 2 dim da memória
	delete[] Raan; // deleta Matriz de 2 dim da memória
	delete[] Rrrn; // deleta Matriz de 2 dim da memória
	delete[] sigma; // deleta Vetor da memória
	delete[] sigmaaa;// deleta Vetor da memória
	delete[] sigmarr; // deleta Vetor da memória
	delete[] Gaa; // deleta Matriz 2 dim da memória
	delete[] Gaan; // deleta Matriz 2 dim da memória
	delete[] Grr; // deleta Matriz 2 dim da memória
	delete[] Grrn; // deleta Matriz 2 dim da memória
	delete[] Haat; // deleta Matriz 2 dim da memória
	delete[] Hrrt; // deleta Matriz 2 dim da memória
	delete[] HtRnaa; // deleta Matriz 2 dim da memória
	delete[] HtRnrr; // deleta Matriz 2 dim da memória
	delete[] zh; // deleta Vetor da memória
	delete[] zhaa; // deleta Vetor da memória
	delete[] zhrr; // deleta Vetor da memória
	delete[] taa; // deleta Vetor da memória
	delete[] trr; // deleta Vetor da memória
	delete[] dtheta; // deleta Vetor da memória
	delete[] dV; // deleta Vetor da memória
	
    return 0;
}
