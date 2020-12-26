//
// Created by Admin on 25.12.2020.
//

#ifndef GINAC_RELAX_PADIC_HPP
#define GINAC_RELAX_PADIC_HPP

#include <vector>
#include <iostream>
#include <cstdlib>

typedef unsigned long long uslong;

class padicRepresentation {
public:
    virtual uslong next(long long i) = 0;
    long long int prec; // long long because arithmetic operations leads to overflow
    std::vector<uslong> coef;
    long long int val;
    uslong prime_base;
    friend std::ostream &operator<<(std::ostream &os, const padicRepresentation &number);
    padicRepresentation() {};

private:
};

class padicNumber : public padicRepresentation {
public:
    padicNumber(long long int base10, uslong prime_base, uslong init_prec = 20);

    uslong next(long long i) override;
private:
    bool is_negative;
    uslong Ox;
    uslong computeOX();

};


class padicSum : public padicRepresentation {
private:
    padicRepresentation &op1;
    padicRepresentation &op2;
    uslong excess = 0;
public:
    uslong next(long long i) override;
    padicSum(padicRepresentation &op1, padicRepresentation &op2);
    void compute_to_N(long long int N);
    void compute_to_max();
    uslong computeSum();
};

#endif //GINAC_RELAX_PADIC_HPP
