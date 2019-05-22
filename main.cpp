#include <iostream>
#include<vector>
#include"../include/Ponto.h"
#include"../include/RungeKutta.h"
#include <fstream>

using namespace std;
double RungeKutta::mu=0.0;
double RungeKutta::alpha=0.0;
double RungeKutta::N=0.0;
double RungeKutta::beta=0.0;
string RungeKutta::tipo="";

int main()
{
    std::cout.precision(15);
    double x0 = 0, y0 = 950,z0 = 50, x = 1400, h = 1;
    RungeKutta::setMu(0.00027472527); // Taxa de mortalidade e natalidade - 1/70*52 - Expectativa de vida brasil em semanas
    RungeKutta::setAlpha(0.0001);//Taxa de contato
    RungeKutta::setN(y0+z0);//População total
    RungeKutta::setBeta(0.125);//Taxa de recuperação - para a tuberculose 1/8 semanas
    string tipo ="SIS";
    RungeKutta::tipo=tipo;//SI ou SIS

    vector<Ponto> sol =  RungeKutta::rungeKutta4(x0, y0, z0, x, h);
    cout<<RungeKutta::tipo;
    cout<<"\nRUNGE KUTTA 4:";
    cout<<"\nR0 :"<<RungeKutta::getR0();
    cout<<"\nComeco:";
    string xSTR="x:[";
    string ySTR="y:[";
    string zSTR="z:[";

    for (int i=0; i<sol.size(); i++)
    {

        Ponto ponto = sol[i];
        if(i+1==sol.size())
        {
            xSTR+=to_string(ponto.getValorX());
            ySTR+=to_string(ponto.getValorY());
            zSTR+=to_string(ponto.getValorZ());
        }
        else
        {
            xSTR+=to_string(ponto.getValorX())+", ";
            ySTR+=to_string(ponto.getValorY())+", ";
            zSTR+=to_string(ponto.getValorZ())+", ";
        }
        cout<<"\nx: " << std::fixed << ponto.getValorX();
        cout<<"\nyn: " << std::fixed << ponto.getValorY();
        cout<<"\nzn: " << std::fixed <<ponto.getValorZ();
        cout<<"\n";

    }

    cout<<"\nfim";
    cout<<"\n--------------------";
    xSTR+="];";
    ySTR+="];";
    zSTR+="];";
    ofstream myfile;
    myfile.open ("example.txt");
    myfile << xSTR+"\n";
    myfile << ySTR+"\n";
    myfile << zSTR+"\n";
    myfile << "plot2d([[discrete, x,y],[discrete, x,z]],[legend,\"S\",\"I\"],[style,[lines,5]]);";
    myfile.close();
}
/*

x' =4x-3y
y' =6x-7y
x(0)=2,y(0)=-1

x(t)=3e^(2t)-e^(-5t)
y(t)=2e^(2t)-3e^(-5t)
EM t=1
x=22.1604303498
y=14.7578983569
 */
