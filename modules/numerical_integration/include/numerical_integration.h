// Copyright 2018 Lvova Alina

#ifndef MODULES_NUMERICAL_INTEGRATION_INCLUDE_NUMERICAL_INTEGRATION_H_
#define MODULES_NUMERICAL_INTEGRATION_INCLUDE_NUMERICAL_INTEGRATION_H_

class Integral {
 public:
    Integral();
    Integral(double lower_limit, double upper_limit, int divisions);
    double RiemannSumLeft();
    double TrapezoidalRule();
    double SimpsonRule();
    double Simpson3_8Rule();
    double BooleRule();
    double NewtonCotes5();
    double GaussianQuadrature();
    void setLower(double _low);
    void setUpper(double _up);
    void setDivisions(int _div);

 private:
    double res;
    double step;
    int div;
    double low, up;
    double function(double x);
};
#endif  // MODULES_NUMERICAL_INTEGRATION_INCLUDE_NUMERICAL_INTEGRATION_H_
