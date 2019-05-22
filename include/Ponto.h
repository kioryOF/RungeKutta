#ifndef PONTO_H
#define PONTO_H


class Ponto
{
private:
    double valorX;
    double valorY;
    double valorZ;
public:
    Ponto();

    double getValorZ();
    double getValorY();
    double getValorX();
    void setValorZ(double);
    void setValorY(double);
    void setValorX(double);

protected:




};

#endif // PONTO_H
