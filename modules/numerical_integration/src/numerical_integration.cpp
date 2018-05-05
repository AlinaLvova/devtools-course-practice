// Copyright 2018 Lvova Alina

#include "include/numerical_integration.h"
#include <math.h>
#include <stdexcept>

double Integral::function(double x) {
    return (x * x + x + 2);
}

Integral::Integral() {
    low = 0;
    up = 0;
    div = 0;
    res = 0;
    step = 0;
}

Integral::Integral(double lower_limit, double upper_limit, int divisions) {
    if (lower_limit < upper_limit) {
        up = upper_limit;
        low = lower_limit;
    } else {
        throw std::out_of_range("limits out of range");
    }
    setDivisions(divisions);
    res = 0;
    step = (upper_limit - lower_limit) / divisions;
}

double Integral::RiemannSumLeft() {
    res = 0;
    double A = low - step;
    for (int i = 0; i < div; i++) {
        A += step;
        res += step * function(A);
    }
    return res;
}

double Integral::TrapezoidalRule() {
    res = 0;
    double A = low - step;
    for (int i = 0; i < div; i++) {
        A += step;
        res += step / 2 * (function(A) + function(A + step));
    }
    return res;
}

double Integral::SimpsonRule() {
    res = 0;
    double A = low - step;
    for (int i = 0; i < div; i++) {
        A += step;
        res += step / 6 *
            (function(A) + 4 * function((2 * A + step) / 2)
             + function(A + step));
    }
    return res;
}

double Integral::Simpson3_8Rule() {
    res = 0;
    double A = low - step;
    for (int i = 0; i < div; i++) {
        A += step;
        double B = A + step;
        res += (B - A) / 8 *
            (function(A) + 3 * function((2 * A + B) / 3)
             + 3 * function((A + 2 * B) / 3) + function(B));
    }
    return res;
}

double Integral::BooleRule() {
    res = 0;
    double A = low - step;
    for (int i = 0; i < div; i++) {
        A += step;
        double B = A + step;
        res += (B - A) / 90 * (7 * function(A) + 32 * function((3 * A + B) / 4)
             + 12 * function((A + B) / 2) +
             32 * function((A + 3 * B) / 4) + 7 * function(B));
    }
    return res;
}

double Integral::NewtonCotes5() {
    res = 0;
    double A = low - step;
    for (int i = 0; i < div; i++) {
        A += step;
        double B = A + step;
        res += (B - A) / 288 * (19 * function(A)
               + 75 * function((4 * A + B) / 5)
               + 50 * function((3 * A + 2 * B) / 5)
               + 50 * function((2 * A + 3 * B) / 5)
               + 75 * function((A + 4 * B) / 5) + 19 * function(B));
    }
    return res;
}

double Integral::GaussianQuadrature() {
    res = 0;
    double A = low - step;
    for (int i = 0; i < div; i++) {
        A += step;
        double B = A + step;
        double f1 = (A + B) / 2;
        double f2 = (B - A) / (2 * sqrt(3));
        res += (B - A) / 2 * (function(f1 - f2) + function(f1 + f2));
    }
    return res;
}

void Integral::setLower(double _low) {
    if (_low < up)
        low = _low;
    else
        throw std::out_of_range("lower limit out of range");
}

void Integral::setUpper(double _up) {
    if (_up > low)
        up = _up;
    else
        throw std::out_of_range("upper limit out of range");
    step = (up - low) / div;
    res = 0;
}

void Integral::setDivisions(int _div) {
    if (_div > 1)
        div = _div;
    else
        throw std::out_of_range("Divisions must be greater than 1");
    step = (up - low) / div;
    res = 0;
}
