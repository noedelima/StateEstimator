//
//  Main.cpp
//  Estimador de Estado
//
//  Created by Noé de Lima Bezerra on 31/12/00.
//  Copyright © 2000 Noé de Lima Bezerra. All rights reserved.
//

#include <iostream>
#include <complex>
#include <cmath>
#include <new>
#include "AuxFunctions.hpp"
using namespace std;

int n;
int m;

void present() {cout << "\n\n\nEstimador de Estado\nNoe de Lima Bezerra\n\nEste programa Recebe os parametros de topologia de uma rede eletrica e os parametros de medicao para estimar o estado mais provavel de operacao do sistema.\n\n";} //Exibe mensagem de apresentação do programa
int selectcase(){
    int option=0;
    cout << "Escolha entre simular uma rede de teste ou entrar com os dados de uma rede do usuario.\n\n";
    cout << "Opcoes:\n";
    cout << "1 para exemplo 2.2 do livro Power System State Estimation, Abur & Exposito;\n";
    cout << "2 para rede de teste 14-bus do IEEE\n";
    cout << "0 ou outro valor para entrar com os parametros da rede do usuario\n";
    cout << "Opcao escolhida: ";
    cin >> option;
	if (option == 1) {
		return 1;
	}
	if (option == 2) {
		cout << "Escolha a redundancia na distribuicao dos medidores:\n1 para 1.5%\n2 para 2.5%\n3 para 3.5%\nOpcao escolhida ";
		cin >> option;
		if (option == 1 || option == 2 || option == 3) { return option + 1; }
		else { return 2; }
	}
    return 0;
} //Exibe uma mensagem para selecionar a utilização do programa

int main() {
    present(); //Apresenta o programa
    
Start:
    int option = selectcase(); //Pede ao usuário definição se utiliza exemplo ou nova rede
    
    if (option==1) {n=3;}
    else if (option==2 || option==3 || option==4) {n=14;}
    else {
        cout << "Quantas barras possui o sistema? ";
        cin >> n;
    }
    
    m=n+1;
    long double *Ybus, *Sz, *Vz;
	Ybus = new long double [(2)*(m)*(m)]; // Matriz de 3 dim
	Sz = new long double [(4)*(m)*(m)]; // Matriz de 3 dim
	Vz = new long double [(2)*(m)]; // Matriz de 2 dim

    if (option==1) {
        GetYbus1(m, Ybus); //Cria a matriz Ybus
        GetMeansPower1(m, Sz); //Recebe os valores de medição de potência
        GetMeansVoltageMag1(m, Vz); //Recebe os valores de magnitude de tensão medidos
    }
    else if (option==2) {
        GetYbus14bus150(m, Ybus); //Cria a matriz Ybus
        GetMeansPower14bus150(m, Sz); //Recebe os valores de medição de potência
        GetMeansVoltageMag14bus150(m, Vz); //Recebe os valores de magnitude de tensão medidos
    }
	else if (option == 3) {
		GetYbus14bus250(m, Ybus); //Cria a matriz Ybus
		GetMeansPower14bus250(m, Sz); //Recebe os valores de medição de potência
		GetMeansVoltageMag14bus250(m, Vz); //Recebe os valores de magnitude de tensão medidos
	}
	else if (option == 4) {
		GetYbus14bus350(m, Ybus); //Cria a matriz Ybus
		GetMeansPower14bus350(m, Sz); //Recebe os valores de medição de potência
		GetMeansVoltageMag14bus350(m, Vz); //Recebe os valores de magnitude de tensão medidos
	}
	else {
        GetYbus(m, Ybus); //Cria a matriz Ybus
        GetMeansPower(m, Sz); //Recebe os valores de medição de potência
        GetMeansVoltageMag(m, Vz); //Recebe os valores de magnitude de tensão medidos
    }

    //Apresenta os valores de entrada
    cout << "\n\n";
    cout <<"Matriz de Admitancia:\n";
    for (int i=1; i<=n; i++) {
        cout << "\nLinha " << i;
        for (int j=1; j<=n; j++) {
			cout << "           " << Ybus[j + m*i + m*m * 0] << "+j(" << Ybus[j + m*i + m*m * 1] << ")";
        }
    } //Mostra a matriz Ybus montada
    cout << "\n\n";
    cout << "\n\nMatriz de Potencias Medidas:\n";
    for (int i=1; i<m; i++) {
        cout << "\nLinha " << i << ":";
        for (int j=1; j<m; j++) {
            cout << "               " << Sz[j + m*i + m*m * 0] << "+j(" << Sz[j + m*i + m*m * 1] << ") +- " << Sz[j + m*i + m*m * 2] << "+-j" << Sz[j + m*i + m*m * 3];
        }
    } //Mostra a matriz Sz montada
    cout << "\n\n";
    cout << "\n\nVetor de Tensoes Medidas:\n\n";
    for (int i=1; i<m; i++) {
		cout << "Vz" << i << " = " << Vz[i + m * 0] << ", com incerteza = +-" << Vz[i + m * 1] << "\n";
    } //Mostra o vetor Vz montado
    cout << "\n\n";

//Analisys:
    int temp=0;
    long double *V, *Theta, *S, e=4;
	V = new long double [m]; // Vetor
	Theta = new long double [m]; // Vetor
	S = new long double [(2)*(m)*(m)]; //Matriz 3 dim

    for (int i=1; i<m; i++) {
        V[i] = 1;
        Theta[i] = 0;
    } //Inicializa o vetor de estado
    
	//Análise de observabilidade
	option = ResidualAnalisys(n, Ybus, V, Theta, Sz, Vz, S);
	cout << "\n\n\n";
	if (option == 1) { cout << "\n\nO sistema nao e observavel nesta configuracao\n\n"; goto End; }

Solution:
	for (int i = 1; i<m; i++) {
		V[i] = 1;
		Theta[i] = 0;
	} //Reinicializa o vetor de estado
    //Definição da tolerância de convergência
    cout << "Especifique a tolerancia de convergencia e na forma 10^(-e): ";
    cin >> e;
    e = pow(10, -e);
	cout << "\nA tolerancia e de " << e << "\n";
    
    //Escolha do método de solução
    cout << "\nEspecifique Metodo de solucao conforme lista abaixo:\n1 para Desacoplado Rapido;\n2 para Desacoplado;\n3 para Newton-Raphson Rapido; ou\noutro valor para Newton-Raphson completo.\n";
    cin >> temp;
    
    if (temp==1) {option = FastDecoupled(n, Ybus, V, Theta, Sz, Vz, S, e);}
    else if (temp==2) {option = Decoupled(n, Ybus, V, Theta, Sz, Vz, S, e);}
    else if (temp==3) {option = FastNewtonRaphson(n, Ybus, V, Theta, Sz, Vz, S, e);}
    else {option = NewtonRaphson(n, Ybus, V, Theta, Sz, Vz, S, e);}
    
//Results:
    if (option==1) {goto End;}
    else if (option>1) {cout << "\n\nO sistema nao e observavel nesta configuracao\n\n"; goto End;}
    //Apresentação da Solução
    cout << "\n\nO estado estimado final e:\n\n";
    for (int i=1; i<m; i++) {
        cout << "\nLinha: " << i;
        for (int j=1; j<m; j++) {
            cout << "           " << S[j + m*i + m*m * 0] << "+j(" << S[j + m*i + m*m * 1] << ")";
        }
    } //Apresenta o estado do sistema
    for (int i=2; i<m; i++) {Theta[i] = Theta[i]*180/3.14156;} //Converte Theta para graus
    cout << "\n\nOs seguintes valores foram obtidos, com V em PU e Theta em graus:\n\n";
    for (int i=1; i<m; i++) {
        cout << "Barra " << i << ": " << V[i] << " V <" << Theta[i] << " deg V\n";
    } //Apresenta as variáveis de estado calculadas

    
End:
    cout << "\n\n\nFim! Por enquanto...\n\n";
    cout << "Digite 1 para nova simulacao, 2 para nova solucao deste sistema ou outro valor para sair: ";
    cin >> temp;
    if (temp==1) {goto Start;}
    else if (temp==2) {goto Solution;}
    cout << "\n\n\nObrigado e volte sempre!!!\n\n\n";

	delete[] Ybus; // deleta Matriz de 3 dim da memória
	delete[] Sz; // deleta Matriz de 3 dim da memória
	delete[] Vz; // deleta Matriz de 2 dim da memória

    return 0;
}
