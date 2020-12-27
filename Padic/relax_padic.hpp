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

long getLn(long long l);

class padicRepresentation {
public:
    virtual long long next(long long i) = 0;

    slong prec; // long long because arithmetic operations leads to overflow
    std::vector<long long> coef;
    slong val;
    uslong prime_base;

    friend std::ostream &operator<<(std::ostream &os, const padicRepresentation &number);

    padicRepresentation() {};

    slong to10base();

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

class padicOperator : public padicRepresentation {
public:
    slong excess = 0;
    padicRepresentation &op1;
    padicRepresentation &op2;
    padicOperator(slong excess, padicRepresentation &op1, padicRepresentation &op2);
    void compute_to_max();
    void compute_to_N(slong N);
};

class padicSum : public padicOperator {
private:
public:
    slong next(long long i) override;
    padicSum(padicRepresentation &op1, padicRepresentation &op2);
    slong computeSum();
};

class padicSub : public padicOperator {
private:
public:
    slong next(long long i) override;
    padicSub(padicRepresentation &op1, padicRepresentation &op2);
    slong computeSub();
};

class padicMul : public padicOperator {
private:
    std::vector<std::vector<long long>> ya;
    std::vector<std::vector<long long>> yb;
public:
    slong next(long long j) override;
    padicMul(padicRepresentation &op1, padicRepresentation &op2);
    slong computeMul();
};


#endif //GINAC_RELAX_PADIC_HPP
