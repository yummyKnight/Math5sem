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

slong getLn(long long l);

class padicRepresentation {
public:
    long long get(long long i);
    padicRepresentation() = default;
public:
    slong prec; // long long because arithmetic operations leads to overflow
    std::vector<long long> coef;
    slong val;
    uslong prime_base;
    bool is_negative;

    friend std::ostream &operator<<(std::ostream &os, const padicRepresentation &number);

    slong to10base();

    slong convertToNegative(long long i);

private:
    virtual long long next() = 0;
};

class padicNumber : public padicRepresentation {
public:
    padicNumber(slong base10, uslong prime_base, uslong init_prec = 20);

    slong next() override;

private:
    slong Ox;
};

class padicOperator : public padicRepresentation {
public:
    slong excess = 0;
    padicRepresentation *op1;
    padicRepresentation *op2;
    void compute_to_max();

    padicOperator(padicRepresentation *op1, padicRepresentation *op2);

    padicOperator();

    void compute_to_N(slong N);

protected:
    slong max_prec = 0;
};

class padicSumSubOperator : public padicOperator {
protected:
    bool mode; // True for sub, false for sum
    slong nextSum();
    slong nextSub();
public:
    padicSumSubOperator(padicRepresentation &op1, padicRepresentation &op2);
private:
    long long int next() override;
};

class padicSum : public padicSumSubOperator {
private:
public:
    padicSum(padicRepresentation &op1, padicRepresentation &op2);
};

class padicSub : public padicSumSubOperator {
private:
public:
    padicSub(padicRepresentation &op1, padicRepresentation &op2);
};

class padicMul : public padicOperator {
private:
    std::vector<std::vector<long long>> ya;
    std::vector<std::vector<long long>> yb;
public:
    slong next() override;
    padicMul(padicRepresentation &op1, padicRepresentation &op2);

};


#endif //GINAC_RELAX_PADIC_HPP
