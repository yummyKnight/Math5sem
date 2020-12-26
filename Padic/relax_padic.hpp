//
// Created by Admin on 25.12.2020.
//

#ifndef GINAC_RELAX_PADIC_HPP
#define GINAC_RELAX_PADIC_HPP

#include <vector>
#include <iostream>
#include <cstdlib>

typedef unsigned long long uslong;
typedef long long int slong;

class padicRepresentation {
public:
    virtual long long next(long long i) = 0;
    slong prec; // long long because arithmetic operations leads to overflow
    std::vector<long long> coef;
    slong val;
    uslong prime_base;
    friend std::ostream &operator<<(std::ostream &os, const padicRepresentation &number);
    padicRepresentation() {};

private:
};

class padicNumber : public padicRepresentation {
public:
    padicNumber(slong base10, uslong prime_base, uslong init_prec = 20);

    slong next(long long i) override;
private:
    bool is_negative;
    slong Ox;
    slong computeOX();

};


class padicSum : public padicRepresentation {
private:
    padicRepresentation &op1;
    padicRepresentation &op2;
    slong excess = 0;
public:
    slong next(long long i) override;
    padicSum(padicRepresentation &op1, padicRepresentation &op2);
    void compute_to_N(slong N);
    void compute_to_max();
    slong computeSum();
};

class padicSub : public padicRepresentation {
private:
    padicRepresentation &op1;
    padicRepresentation &op2;
    slong excess = 0;
public:
    slong next(long long i) override;
    padicSub(padicRepresentation &op1, padicRepresentation &op2);
    void compute_to_N(slong N);
    void compute_to_max();
    slong computeSub();
};

#endif //GINAC_RELAX_PADIC_HPP
