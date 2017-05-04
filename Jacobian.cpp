//
//  Jacobian.cpp
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

// Todos os arrays de entrada, vetores e matrizes, são enviados como endereço do elemento 0 e recebidos como ponteiros
// Os índices dos arrays multidimensionais são calculados com base na disposição em memória utilizando lógica de base numérica
// Por exemplo, o elemento Matriz[i][j][k] será Matriz[k+n*j+n*n*k], onde n é a ordem da matriz
// Para matriz não "quadrada", com dimensão Matriz[a][b][c], o elemento Matriz[i][j][k] será Matriz[(k)+(c*j)+(c*b*i)]


// Funções de derivadas para os elementos internos da matriz jacobiana

// Elementos correspondentes às potências nas barras

long double dPiTi(int m, int i, long double* Ybus, long double* V, long double* Theta) {
    long double G, B, dPiTi=0;
    for (int k=1; k<m; k++) {
        G = Ybus[k+m*i]; //Parte real
        B = Ybus[k+m*i+m*m]; //Parte imaginária
        dPiTi = dPiTi + V[i]*V[k]*(-G*sin(Theta[i]-Theta[k])+B*cos(Theta[i]-Theta[k]));
    }
    B = Ybus[i+m*i+m*m];
    dPiTi = dPiTi - pow(V[i], 2)*B;
    return dPiTi;
} //Calcula a derivada da potência ativa na barra i em relação a teta na barra i

long double dQiTi(int m, int i, long double* Ybus, long double* V, long double* Theta) {
    long double G, B, dQiTi=0;
    for (int k=1; k<m; k++) {
        G = Ybus[k+m*i]; //Parte real
        B = Ybus[k+m*i+m*m]; //Parte imaginária
        dQiTi = dQiTi + V[i]*V[k]*(G*cos(Theta[i]-Theta[k])+B*sin(Theta[i]-Theta[k]));
    }
    G = Ybus[i+m*i];
    dQiTi = dQiTi - pow(V[i], 2)*G;
    return dQiTi;
} //Calcula a derivada da potência reativa na barra i em relação a teta na barra i

long double dPiTj(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double G, B, dPiTj=0;
    G = Ybus[j+m*i]; //Parte real
    B = Ybus[j+m*i+m*m]; //Parte imaginária
    dPiTj = V[i]*V[j]*(G*sin(Theta[i]-Theta[j])-B*cos(Theta[i]-Theta[j]));
    return dPiTj;
} //Calcula a derivada da potência ativa na barra i em relação a teta na barra j

long double dQiTj(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double G, B, dQiTj=0;
    G = Ybus[j+m*i]; //Parte real
    B = Ybus[j+m*i+m*m]; //Parte imaginária
    dQiTj = V[i]*V[j]*(-G*cos(Theta[i]-Theta[j])-B*sin(Theta[i]-Theta[j]));
    return dQiTj;
} //Calcula a derivada da potência reativa na barra i em relação a teta na barra j

long double dPiVi(int m, int i, long double* Ybus, long double* V, long double* Theta) {
    long double G, B, dPiVi=0;
    for (int k=1; k<m; k++) {
        G = Ybus[k+m*i]; //Parte real
        B = Ybus[k+m*i+m*m]; //Parte imaginária
        dPiVi = dPiVi + V[k]*(G*cos(Theta[i]-Theta[k])+B*sin(Theta[i]-Theta[k]));
    }
    G = Ybus[i+m*i];
    dPiVi = dPiVi + V[i]*G;
    return dPiVi;
} //Calcula a derivada da potência ativa na barra i em relação a tensão na barra i

long double dQiVi(int m, int i, long double* Ybus, long double* V, long double* Theta) {
    long double G, B, dQiVi=0;
    for (int k=1; k<m; k++) {
        G = Ybus[k+m*i]; //Parte real
        B = Ybus[k+m*i+m*m]; //Parte imaginária
        dQiVi = dQiVi + V[k]*(G*sin(Theta[i]-Theta[k])-B*cos(Theta[i]-Theta[k]));
    }
    B = Ybus[i+m*i+m*m];
    dQiVi = dQiVi - V[i]*B;
    return dQiVi;
} //Calcula a derivada da potência reativa na barra i em relação a tensão na barra i

long double dPiVj(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double G, B, dPiVj=0;
    G = Ybus[j+m*i]; //Parte real
    B = Ybus[j+m*i+m*m]; //Parte imaginária
    dPiVj = V[i]*(G*cos(Theta[i]-Theta[j])+B*sin(Theta[i]-Theta[j]));
    return dPiVj;
} //Calcula a derivada da potência ativa na barra i em relação a tesão na barra j

long double dQiVj(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double G, B, dQiVj=0;
    G = Ybus[j+m*i]; //Parte real
    B = Ybus[j+m*i+m*m]; //Parte imaginária
    dQiVj = V[i]*(G*sin(Theta[i]-Theta[j])-B*cos(Theta[i]-Theta[j]));
    return dQiVj;
} //Calcula a derivada da potência reativa na barra i em relação a tesão na barra j

// Elementos correspondentes às potências nas linhas

long double dPijTi(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double dPijTi, g, b;
    g = -Ybus[j+m*i]; //Parte real
    b = -Ybus[j+m*i+m*m]; //Parte imaginária
    dPijTi = V[i]*V[j]*(g*sin(Theta[i]-Theta[j])-b*cos(Theta[i]-Theta[j]));
    return dPijTi;
} //Calcula a derivada da potência ativa injetada da barra i para a barra j em relação a teta na barra i

long double dQijTi(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double dQijTi, g, b;
    g = -Ybus[j+m*i]; //Parte real
    b = -Ybus[j+m*i+m*m]; //Parte imaginária
    dQijTi = -V[i]*V[j]*(g*cos(Theta[i]-Theta[j])-b*sin(Theta[i]-Theta[j]));
    return dQijTi;
} //Calcula a derivada da potência reativa injetada da barra i para a barra j em relação a teta na barra i

long double dPijTj(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double dPijTj, g, b;
    g = -Ybus[j+m*i]; //Parte real
    b = -Ybus[j+m*i+m*m]; //Parte imaginária
    dPijTj = -V[i]*V[j]*(g*sin(Theta[i]-Theta[j])-b*cos(Theta[i]-Theta[j]));
    return dPijTj;
} //Calcula a derivada da potência ativa injetada da barra i para a barra j em relação a teta na barra j

long double dQijTj(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double dQijTj, g, b;
    g = -Ybus[j+m*i]; //Parte real
    b = -Ybus[j+m*i+m*m]; //Parte imaginária
    dQijTj = V[i]*V[j]*(g*cos(Theta[i]-Theta[j])-b*sin(Theta[i]-Theta[j]));
    return dQijTj;
} //Calcula a derivada da potência reativa injetada da barra i para a barra j em relação a teta na barra j

long double dPijVi(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double dPijVi, g, b;
    g = -Ybus[j+m*i]; //Parte real
    b = -Ybus[j+m*i+m*m]; //Parte imaginária
    dPijVi = -V[j]*(g*cos(Theta[i]-Theta[j])+b*sin(Theta[i]-Theta[j])) + 2*g*V[i];
    return dPijVi;
} //Calcula a derivada da potência ativa injetada da barra i para a barra j em relação a tensão na barra i

long double dQijVi(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double dQijVi, g, b;
    g = -Ybus[j+m*i]; //Parte real
    b = -Ybus[j+m*i+m*m]; //Parte imaginária
    dQijVi = -V[j]*(g*sin(Theta[i]-Theta[j])-b*cos(Theta[i]-Theta[j])) - 2*b*V[i];
    return dQijVi;
} //Calcula a derivada da potência reativa injetada da barra i para a barra j em relação a tensão na barra i

long double dPijVj(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double dPijVj, g, b;
    g = -Ybus[j+m*i]; //Parte real
    b = -Ybus[j+m*i+m*m]; //Parte imaginária
    dPijVj = -V[i]*(g*cos(Theta[i]-Theta[j])+b*sin(Theta[i]-Theta[j]));
    return dPijVj;
} //Calcula a derivada da potência ativa injetada da barra i para a barra j em relação a tensão na barra j

long double dQijVj(int m, int i, int j, long double* Ybus, long double* V, long double* Theta) {
    long double dQijVj, g, b;
    g = -Ybus[j+m*i]; //Parte real
    b = -Ybus[j+m*i+m*m]; //Parte imaginária
    dQijVj = -V[i]*(g*sin(Theta[i]-Theta[j])-b*cos(Theta[i]-Theta[j]));
    return dQijVj;
} //Calcula a derivada da potência reativa injetada da barra i para a barra j em relação a tensão na barra j

// Verificando a quantidade de linhas da matriz jacobiana (número de medidores)

int DimensionH(int m, long double* Sz, long double* Vz) {
    int k=0;
    long double P, Q;
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            P = Sz[j+m*i]; //Medição de potência ativa
            Q = Sz[j+m*i+m*m]; //Medição de potência reativa
            if (P!=0) {k=k+1;}
            if (Q!=0) {k=k+1;}
        }
        if (Vz[i]!=0) {
            k=k+1;
        }
    }
    return k;
} //Calcula o número de elementos diferentes de 0 na matriz de medição de potência e no vetor de tensões medidas


int DimensionHaa(int m, long double* Sz) {
    int k=0;
    long double P;
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            P = Sz[j+m*i]; //Medição de potência ativa
            if (P!=0) {k=k+1;}
        }
    }
    return k;
} //Calcula o número de elementos diferentes de 0 na matriz de medição de potência apenas para medidas de potência ativa


int DimensionHrr(int m, long double* Sz, long double* Vz) {
    int k=0;
    long double Q;
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            Q = Sz[j+m*i+m*m]; //Medição de potência reativa
            if (Q!=0) {k=k+1;}
        }
        if (Vz[i]!=0) {
            k=k+1;
        }
    }
    return k;
} //Calcula o número de elementos diferentes de 0 na matriz de medição de potência apenas para medidas de potência reativa


int DimensionHaaObs(int m, int layer, long double* Sz, long double* zaa) {
	int k = 0;
	long double P;
	for (int i = 1; i<m; i++) {
		for (int j = 1; j<m; j++) {
			P = Sz[j + m*i];
			if (i != j && P != 0) {
				k = k + 1;
				zaa[k + 0 * layer] = P;
				zaa[k + 1 * layer] = i;
				zaa[k + 2 * layer] = j;
			}
		}
	} //Medição de fluxo de potência ativa
	for (int i = 1; i<m; i++) {
		for (int j = 1; j<m; j++) {
			P = Sz[j + m*i];
			if (i == j && P != 0) {
				k = k + 1;
				zaa[k + 0 * layer] = P;
				zaa[k + 1 * layer] = i;
				zaa[k + 2 * layer] = j;
			}
		}
	} //Medição de injeção de potência ativa
	return k;
} //Calcula o número de elementos diferentes de 0 na matriz de medição de potência apenas para medidas de potência ativa


// Calculo da matriz jacobiana

int Jacobian(int m, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* H) {
    int line=0, n=2*m-1;
    long double P, Q;
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            P = Sz[j+m*i];
            if (i!=j && P!=0) {
                line = line+1;
                for (int row=1; row<m; row++) {
                    if (row==i) {
                        H[row+n*line] = dPijTi(m, i, j, Ybus, V, Theta);
                        //std::cout << "Elemento P-T-iji (" << line << "," << row << "): " << H[row+n*line] << " Função: dPijTi\n";
                        H[row+m-1+n*line] = dPijVi(m, i, j, Ybus, V, Theta);
						//std::cout << "Elemento P-V-iji (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dPijVi\n";
                    }
                    else if (row==j) {
                        H[row+n*line] = dPijTj(m, i, j, Ybus, V, Theta);
						//std::cout << "Elemento P-T-ijj (" << line << "," << row << "): " << H[row+n*line] << " Função: dPijTj\n";
                        H[row+m-1+n*line] = dPijVj(m, i, j, Ybus, V, Theta);
						//std::cout << "Elemento P-V-ijj (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dPijVj\n";
                    }
                    else {
                        H[row+n*line] = 0;
						//std::cout << "Elemento P-T-ijk (" << line << "," << row << "): " << H[row+n*line] << " Função: 0\n";
                        H[row+m-1+n*line] = 0;
						//std::cout << "Elemento P-V-ijk (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: 0\n";
                    }
                }
            }
        }
    } //Para os medidores de fluxo de potência ativa nas linhas
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
                P = Sz[j+m*i];
            if (i==j && P!=0) {
                    line = line+1;
                for (int row=1; row<m; row++) {
                    if (row==i) {
                        H[row+n*line] = dPiTi(m, i, Ybus, V, Theta);
						//std::cout << "Elemento P-T-iii (" << line << "," << row << "): " << H[row+n*line] << " Função: dPiTi\n";
                        H[row+m-1+n*line] = dPiVi(m, i, Ybus, V, Theta);
						//std::cout << "Elemento P-V-iii (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dPiVi\n";
                    }
                    else {
                        H[row+n*line] = dPiTj(m, i, row, Ybus, V, Theta);
						//std::cout << "Elemento P-T-iij (" << line << "," << row << "): " << H[row+n*line] << " Função: dPiTj\n";
                        H[row+m-1+n*line] = dPiVj(m, i, row, Ybus, V, Theta);
						//std::cout << "Elemento P-V-iij (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dPiVj\n";
                    }
                }
            }
        }
    } //Para os medidores de potência ativa nas barras
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            Q = Sz[j+m*i+m*m];
            if (i!=j && Q!=0) {
                line = line+1;
                for (int row=1; row<m; row++) {
                    if (row==i) {
                        H[row+n*line] = dQijTi(m, i, j, Ybus, V, Theta);
						//std::cout << "Elemento Q-T-iji (" << line << "," << row << "): " << H[row+n*line] << " Função: dQijTi\n";
                        H[row+m-1+n*line] = dQijVi(m, i, j, Ybus, V, Theta);
						//std::cout << "Elemento Q-V-iji (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dQijVi\n";
                    }
                    else if (row==j) {
                        H[row+n*line] = dQijTj(m, i, j, Ybus, V, Theta);
						//std::cout << "Elemento Q-T-ijj (" << line << "," << row << "): " << H[row+n*line] << " Função: dQijTj\n";
                        H[row+m-1+n*line] = dQijVj(m, i, j, Ybus, V, Theta);
						//std::cout << "Elemento Q-V-ijj (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dQijVj\n";
                    }
                    else {
                        H[row+n*line] = 0;
						//std::cout << "Elemento Q-T-ijk (" << line << "," << row << "): " << H[row+n*line] << " Função: 0\n";
                        H[row+m-1+n*line] = 0;
						//std::cout << "Elemento Q-V-ijk (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: 0\n";
                    }
                }
            }
        }
    } //Para os medidores de fluxo de potência reativa nas linhas
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            Q = Sz[j+m*i+m*m];
            if (i==j && Q!=0) {
                line = line+1;
                for (int row=1; row<m; row++) {
                    if (row==i) {
                        H[row+n*line] = dQiTi(m, i, Ybus, V, Theta);
						//std::cout << "Elemento Q-T-iii (" << line << "," << row << "): " << H[row+n*line] << " Função: dQiTi\n";
                        H[row+m-1+n*line] = dQiVi(m, i, Ybus, V, Theta);
						//std::cout << "Elemento Q-V-iii (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dQiVi\n";
                    }
                    else {
                        H[row+n*line] = dQiTj(m, i, row, Ybus, V, Theta);
						//std::cout << "Elemento Q-T-iij (" << line << "," << row << "): " << H[row+n*line] << " Função: dQiTj\n";
                        H[row+m-1+n*line] = dQiVj(m, i, row, Ybus, V, Theta);
						//std::cout << "Elemento Q-V-iij (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dQiVj\n";
                    }
                }
            }
        }
    } //Para os medidores de potência reativa nas barras
    for (int i=1; i<m; i++) {
        if (Vz[i]!=0) {
            line = line+1;
            for (int row=1; row<m; row++) {
                H[row+n*line]=0;
                H[row+m-1+n*line]=0;
				//std::cout << "Elemento V-T-ij (" << line << "," << row << "): " << H[row+n*line] << "\n";
                if (row==i) {
                    H[row+m-1+n*line] = 1;
					//std::cout << "Elemento V-V-ii (" << line << "," << row << "): " << H[row+m-1+n*line] << "\n";
                }
            }
        }
    } //Para os medidores de módulo de tensão nas barras
    return 0;
}


void Residual(int m, long double* Sz, long double* Vz, long double* sigma) {
    int line=0;
    long double P, Q;
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            P = Sz[j+m*i];
            if (i!=j && P!=0) {
                line = line+1;
                sigma[line] = Sz[j+m*i+m*m*2];
            }
        }
    } //Para os medidores de fluxo de potência ativa nas linhas
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            P = Sz[j+m*i];
            if (i==j && P!=0) {
                line = line+1;
                sigma[line] = Sz[j+m*i+m*m*2];
            }
        }
    } //Para os medidores de potência ativa nas barras
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            Q = Sz[j+m*i+m*m];
            if (i!=j && Q!=0) {
                line = line+1;
                sigma[line] = Sz[j+m*i+m*m*3];
            }
        }
    } //Para os medidores de fluxo de potência reativa nas linhas
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            Q = Sz[j+m*i+m*m];
            if (i==j && Q!=0) {
                line = line+1;
                sigma[line] = Sz[j+m*i+m*m*3];
            }
        }
    } //Para os medidores de potência reativa nas barras
    for (int i=1; i<m; i++) {
        if (Vz[i]!=0) {
            line = line+1;
            sigma[line] = Vz[i+m];
        }
    } //Para os medidores de módulo de tensão nas barras
    return;
}

long double Max(int m, long double* Sz, long double* S, long double* Vz, long double* V, long double* zh) {
    long double Max=0;
    int line=0, module;
    long double P, Q;
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            P = Sz[j+m*i];
            if (i!=j && P!=0) {
                line = line+1;
                zh[line] = P-S[j+m*i];
                module = sqrt(pow(zh[line], 2));
                if (module>Max) {Max = module;}
            }
        }
    } //Para os medidores de fluxo de potência ativa nas linhas
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            P = Sz[j+m*i];
            if (i==j && P!=0) {
                line = line+1;
                zh[line] = P-S[j+m*i];
                module = sqrt(pow(zh[line], 2));
                if (module>Max) {Max = module;}
            }
        }
    } //Para os medidores de potência ativa nas barras
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            Q = Sz[j+m*i+m*m];
            if (i!=j && Q!=0) {
                line = line+1;
                zh[line] = Q-S[j+m*i+m*m];
                module = sqrt(pow(zh[line], 2));
                if (module>Max) {Max = module;}
            }
        }
    } //Para os medidores de fluxo de potência reativa nas linhas
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            Q = Sz[j+m*i+m*m];
            if (i==j && Q!=0) {
                line = line+1;
                zh[line] = Q-S[j+m*i+m*m];
                module = sqrt(pow(zh[line], 2));
                if (module>Max) {Max = module;}
            }
        }
    } //Para os medidores de potência reativa nas barras
    for (int i=1; i<m; i++) {
        if (Vz[i]!=0) {
            line = line+1;
            zh[line] = Vz[i]-V[i];
            module = sqrt(pow(zh[line], 2));
            if (module>Max) {Max = module;}
        }
    } //Para os medidores de módulo de tensão nas barras
    return Max;
}



int DecoupledJacobianAA(int m, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Haa) {
    int line=0, n=m-1;
    long double P;
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            P = Sz[j+m*i];
            if (i!=j && P!=0) {
                line = line+1;
                for (int row=1; row<m; row++) {
                    if (row==i) {
                        Haa[row-1+n*line] = dPijTi(m, i, j, Ybus, V, Theta);
                        //cout << "Elemento P-T-iji (" << line << "," << row << "): " << H[row+n*line] << " Função: dPijTi\n";
                    }
                    else if (row==j) {
                        Haa[row-1+n*line] = dPijTj(m, i, j, Ybus, V, Theta);
                        //cout << "Elemento P-T-ijj (" << line << "," << row << "): " << H[row+n*line] << " Função: dPijTj\n";
                    }
                    else {
                        Haa[row-1+n*line] = 0;
                        //cout << "Elemento P-T-ijk (" << line << "," << row << "): " << H[row+n*line] << " Função: 0\n";
                    }
                }
            }
        }
    } //Para os medidores de fluxo de potência ativa nas linhas
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            P = Sz[j+m*i];
            if (i==j && P!=0) {
                line = line+1;
                for (int row=1; row<m; row++) {
                    if (row==i) {
                        Haa[row-1+n*line] = dPiTi(m, i, Ybus, V, Theta);
                        //cout << "Elemento P-T-iii (" << line << "," << row << "): " << H[row+n*line] << " Função: dPiTi\n";
                    }
                    else {
                        Haa[row-1+n*line] = dPiTj(m, i, row, Ybus, V, Theta);
                        //cout << "Elemento P-T-iij (" << line << "," << row << "): " << H[row+n*line] << " Função: dPiTj\n";
                    }
                }
            }
        }
    } //Para os medidores de potência ativa nas barras
    return line;
}


int DecoupledJacobianRR(int m, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Vz, long double* Hrr) {
    int line=0;
    long double Q;
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            Q = Sz[j+m*i+m*m];
            if (i!=j && Q!=0) {
                line = line+1;
                for (int row=1; row<m; row++) {
                    if (row==i) {
                        Hrr[row+m*line] = dQijVi(m, i, j, Ybus, V, Theta);
                        //cout << "Elemento Q-V-iji (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dQijVi\n";
                    }
                    else if (row==j) {
                        Hrr[row+m*line] = dQijVj(m, i, j, Ybus, V, Theta);
                        //cout << "Elemento Q-V-ijj (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dQijVj\n";
                    }
                    else {
                        Hrr[row+m*line] = 0;
                        //cout << "Elemento Q-V-ijk (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: 0\n";
                    }
                }
            }
        }
    } //Para os medidores de fluxo de potência reativa nas linhas
    for (int i=1; i<m; i++) {
        for (int j=1; j<m; j++) {
            Q = Sz[j+m*i+m*m];
            if (i==j && Q!=0) {
                line = line+1;
                for (int row=1; row<m; row++) {
                    if (row==i) {
                        Hrr[row+m*line] = dQiVi(m, i, Ybus, V, Theta);
                        //cout << "Elemento Q-V-iii (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dQiVi\n";
                    }
                    else {
                        Hrr[row+m*line] = dQiVj(m, i, row, Ybus, V, Theta);
                        //cout << "Elemento Q-V-iij (" << line << "," << row << "): " << H[row+m-1+n*line] << " Função: dQiVj\n";
                    }
                }
            }
        }
    } //Para os medidores de potência reativa nas barras
    for (int i=1; i<m; i++) {
        if (Vz[i]!=0) {
            line = line+1;
            for (int row=1; row<m; row++) {
                Hrr[row+m*line]=0;
                //cout << "Elemento V-T-ij (" << line << "," << row << "): " << H[row+n*line] << "\n";
                if (row==i) {
                    Hrr[row+m*line] = 1;
                    //cout << "Elemento V-V-ii (" << line << "," << row << "): " << H[row+m-1+n*line] << "\n";
                }
            }
        }
    } //Para os medidores de módulo de tensão nas barras
    return line;
}


int DecoupledJacobianObs(int m, long double* Ybus, long double* V, long double* Theta, long double* Sz, long double* Haa) {
	int line = 0, n = m - 1;
	long double P;
	for (int i = 1; i<m; i++) {
		for (int j = 1; j<m; j++) {
			P = Sz[j + m*i];
			if (i != j && P != 0) {
				line = line + 1;
				for (int row = 1; row<m; row++) {
					if (row == i) {
						Haa[row - 1 + n*line] = 1; // dPijTi(m, i, j, Ybus, V, Theta);
						//std::cout << "Elemento P-T-iji (" << line << "," << row-1 << "): " << Haa[row-1+n*line] << " Funcao: dPijTi = 1\n";
					}
					else if (row == j) {
						Haa[row - 1 + n*line] = -1; // dPijTj(m, i, j, Ybus, V, Theta);
						//std::cout << "Elemento P-T-ijj (" << line << "," << row-1 << "): " << Haa[row-1+n*line] << " Funcao: dPijTj = -1\n";
					}
					else {
						Haa[row - 1 + n*line] = 0;
						//std::cout << "Elemento P-T-ijk (" << line << "," << row-1 << "): " << Haa[row-1+n*line] << " Funcao: 0\n";
					}
				}
			}
		}
	} //Para os medidores de fluxo de potência ativa nas linhas
	for (int i = 1; i<m; i++) {
		for (int j = 1; j<m; j++) {
			P = Sz[j + m*i];
			if (i == j && P != 0) {
				line = line + 1;
				for (int row = 1; row<m; row++) {
					if (row == i) {
						Haa[row - 1 + n*line] = 0;
						for (int loop = 1; loop < m; loop++) { if (Ybus[loop + m*i + 1 * m*m] != 0 && loop != i) { Haa[row - 1 + n*line] += 1; } }
						//Haa[row - 1 + n*line] = dPiTi(m, i, Ybus, V, Theta);
						//std::cout << "Elemento P-T-iii (" << line << "," << row-1 << "): " << Haa[row-1+n*line] << " Funcao: dPiTi = soma\n";
					}
					else {
						//Haa[row - 1 + n*line] = Ybus[j + m*i + 1 * m*m];
						if (Ybus[row + m*i + 1 * m*m] != 0) { Haa[row -1 + n*line] = -1; }
						else { Haa[row - 1 + n*line] = 0; }
						//Haa[row - 1 + n*line] = dPiTj(m, i, row, Ybus, V, Theta);
						//std::cout << "Elemento P-T-iij (" << line << "," << row-1 << "): " << Haa[row-1+n*line] << " Funcao: dPiTj = Ybus\n";
					}
				}
			}
		}
	} //Para os medidores de potência ativa nas barras
	return line;
}
