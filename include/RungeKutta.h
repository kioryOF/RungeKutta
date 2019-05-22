#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H
#include"../include/Ponto.h"
#include<vector>
#include <iostream>
#include <fstream>

using namespace std;
class RungeKutta
{
public:
    RungeKutta();

    static vector<Ponto> rungeKutta1(double, double, double, double, double);
    static vector<Ponto> rungeKutta2(double,double,double,double,double);
    static vector<Ponto> rungeKutta3(double,double,double,double,double);
    static vector<Ponto> rungeKutta4(double,double,double,double,double);

    static double getR0();
    static double setMu(double);
    static double setAlpha(double);
    static double setN(double);
    static double setBeta(double);
    static string setTipo(string);
    static string tipo;

protected:

private:
    static double mu;
    static double alpha;
    static double N;
    static double beta;

    static double funcaoF(double,double,double);
    static double funcaoG(double,double,double);
};

#endif // RUNGEKUTTA_H

