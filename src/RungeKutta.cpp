#include "RungeKutta.h"
#include <iostream>
#include<vector>
#include"../include/Ponto.h"
using namespace std;

RungeKutta::RungeKutta()
{

}


vector<Ponto> RungeKutta::rungeKutta1(double x0, double y0, double z0, double x, double h)
{
    double variavelConversao = ((x - x0) / h);
    int n= variavelConversao;

    vector<Ponto> pontos;
    double k1;
    double l1;

    double y = y0;
    double z = z0;
    for (int i = 1; i <= n; i++)
    {

        k1 = RungeKutta::funcaoF(x0, y, z);
        l1 = RungeKutta::funcaoG(x0, y, z);

        y = y + h * (k1);
        z = z + h * (l1);
        x0 = x0 + h;

        Ponto ponto ;
        ponto.setValorX(x0);
        ponto.setValorY(y);
        ponto.setValorZ(z);
        pontos.push_back(ponto);
    }
    return pontos;
}

vector<Ponto> RungeKutta::rungeKutta2(double x0, double y0, double z0, double x, double h)
{
    double variavelConversao = ((x - x0) / h);
    int n= variavelConversao;
    vector<Ponto> pontos ;
    double k1, k2;
    double l1, l2;

    double y = y0;
    double z = z0;
    for (int i = 1; i <= n; i++)
    {

        k1 = RungeKutta::funcaoF(x0, y, z);
        l1 = RungeKutta::funcaoG(x0, y, z);

        k2 = RungeKutta::funcaoF(x0 + h, y + h * k1, z + h * k1);
        l2 = RungeKutta::funcaoG(x0 + h, y + h * l1, z + h * l1);

        y = y + h * (k1 + k2) / 2;
        z = z + h * (l1 + l2) / 2;
        x0 = x0 + h;

        Ponto ponto;
        ponto.setValorX(x0);
        ponto.setValorY(y);
        ponto.setValorZ(z);
        pontos.push_back(ponto);
    }
    return pontos;
}

vector<Ponto> RungeKutta::rungeKutta3(double x0, double y0, double z0, double x, double h)
{
    double variavelConversao = ((x - x0) / h);
    int n= variavelConversao;

    vector<Ponto> pontos;
    double k1, k2, k3;
    double l1, l2, l3;

    double y = y0;
    double z = z0;
    for (int i = 1; i <= n; i++)
    {

        k1 = RungeKutta::funcaoF(x0, y, z);
        l1 = RungeKutta::funcaoG(x0, y, z);

        k2 = RungeKutta::funcaoF(x0 + h / 3, y + h * k1 / 3, z + h * l1 / 3);
        l2 = RungeKutta::funcaoG(x0 + h / 3, y + h * k1 / 3, z + h * l1 / 3);

        k3 = RungeKutta::funcaoF(x0 + 2 * h / 3, y + h * 2 * k2 / 3, z + h * 2 * l2 / 3);
        l3 = RungeKutta::funcaoG(x0 + 2 * h / 3, y + h * 2 * k2 / 3, z + h * 2 * l2 / 3);

        y = y + h * (k1 + 3 * k3) / 4;
        z = z + h * (l1 + 3 * l3) / 4;
        x0 = x0 + h;

        Ponto ponto;
        ponto.setValorX(x0);
        ponto.setValorY(y);
        ponto.setValorZ(z);
        pontos.push_back(ponto);
    }
    return pontos;
}

vector<Ponto> RungeKutta::rungeKutta4(double x0, double y0, double z0, double x, double h)
{
    double variavelConversao = ((x - x0) / h);
    int n= variavelConversao;
    vector<Ponto> pontos;
    double k1, k2, k3, k4;
    double l1, l2, l3, l4;

    double y = y0;
    double z = z0;
    Ponto ponto;
    ponto.setValorX(x0);
    ponto.setValorY(y);
    ponto.setValorZ(z);
    pontos.push_back(ponto);
    for (int i = 1; i <= n; i++)
    {

        k1 = RungeKutta::funcaoF(x0, y, z);
        l1 = RungeKutta::funcaoG(x0, y, z);

        k2 = RungeKutta::funcaoF(x0 + h / 2, y + h * k1 / 2, z + h * l1 / 2);
        l2 = RungeKutta::funcaoG(x0 + h / 2, y + h * k1 / 2, z + h * l1 / 2);

        k3 = RungeKutta::funcaoF(x0 + h / 2, y + h * k2 / 2, z + h * l2 / 2);
        l3 = RungeKutta::funcaoG(x0 + h / 2, y + h * k2 / 2, z + h * l2 / 2);

        k4 = RungeKutta::funcaoF(x0 + h, y + h * k3, z + h * l3);
        l4 = RungeKutta::funcaoG(x0 + h, y + h * k3, z + h * l3);

        y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        z = z + h * (l1 + 2 * l2 + 2 * l3 + l4) / 6;
        x0 = x0 + h;


        ponto.setValorX(x0);
        ponto.setValorY(y);
        ponto.setValorZ(z);
        pontos.push_back(ponto);

    }

    return pontos;
}
double RungeKutta::funcaoF(double x, double y, double z)
{
    if(tipo=="SI")
    {
        return mu*N-alpha*y*z-mu*y;
    }
    else
    {
        if(tipo=="SIS")
        {
            return mu*N-alpha*y*z-mu*y+beta*z;
        }

    }

}

double RungeKutta::funcaoG(double x, double y, double z)
{
    if(tipo=="SI")
    {
        return alpha*y*z-mu*z;
    }
    else
    {
        if(tipo=="SIS")
        {
            return alpha*y*z-beta*z-mu*z;
        }
    }

}

double RungeKutta::setMu(double valorMu)
{
    mu=valorMu;
}

double RungeKutta::setAlpha(double valorAlpha)
{
    alpha=valorAlpha;

}

double RungeKutta::setN(double valorN)
{
    N=valorN;
}

double RungeKutta::setBeta(double valorBeta)
{
    beta=valorBeta;
}
string RungeKutta::setTipo(string valorTipo)
{
    tipo=valorTipo;
}
double RungeKutta::getR0()
{
     if(tipo=="SI")
    {
        return alpha*N/mu;
    }
    else
    {
        if(tipo=="SIS")
        {
            return alpha*N/(mu+beta);
        }
    }
}

