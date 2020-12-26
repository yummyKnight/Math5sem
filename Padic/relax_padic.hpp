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
    virtual uslong next(uslong i) = 0;

    uslong prec;
    std::vector<uslong> coef;
    long long val;
    uslong prime_base;
private:
};

class padicNumber : private padicRepresentation {
public:
    padicNumber(long long int base10, uslong prime_base, uslong init_prec = 20);

    uslong next(uslong i) override;

    friend std::ostream &operator<<(std::ostream &os, const padicNumber &number);

private:
    bool is_negative;
    uslong Ox;

    friend class padicSum;

    uslong computeOX();

};


class padicSum : public padicRepresentation {
private:
    padicRepresentation &op1;
    padicRepresentation &op2;
    uslong val_diff;
    uslong excess = 0;
public:
    uslong next(uslong i) override;

    padicSum(padicRepresentation &op1, padicRepresentation &op2);

    uslong computeSum();
};

#endif //GINAC_RELAX_PADIC_HPP
