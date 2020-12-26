//
// Created by Admin on 25.12.2020.
//

#include "relax_padic.hpp"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>

using namespace std;

bool is_prime(uslong n) {
    for (int i = 2; i * i <= n; ++i) //no sqrt, please
    {
        if (n % i == 0) return false;
    }
    return true;
}

uslong reduce(uslong &x, uslong prime_base) {
    uslong d;
    d = x % prime_base;
    x -= d;
    x = x / prime_base;
    return d;
}


padicNumber::padicNumber(long long base10, uslong prime_base, uslong init_prec) {
    uslong x;
    uslong d = 0;
    uslong j;
    long long val = 0;
    std::vector<uslong> coef;
    is_negative = base10 < 0;
    x = base10;
    if (!is_prime(prime_base))
        throw std::invalid_argument("P should be prime");

    d = reduce(x, prime_base);
    while (d == 0) {
        val++;
        d = reduce(x, prime_base);
    }
    if (is_negative) {
        coef.push_back(prime_base - d);
    } else {
        coef.push_back(d);
    }
    for (j = 0; x != 0 || j > init_prec; j++) {
        d = reduce(x, prime_base);
        if (is_negative) {
            coef.push_back(prime_base - d - 1);
        } else {
            coef.push_back(d);
        }
    }

    this->prec = coef.size();
    this->Ox = x;
    this->prime_base = prime_base;
    this->coef = coef;
    this->val = val;
}

ostream &operator<<(ostream &os, const padicNumber &number) {
    uslong i = 0;
    auto coef = number.coef.at(i);
    if (!coef) {
        os << 0;
        return os;
    }
    cout << number.coef.at(i) << "*" << std::to_string(number.prime_base) + "^" + std::to_string(i + number.val);

    for (i = i + 1; i < number.coef.size(); i++) {
        coef = number.coef.at(i);
        if (coef)
            cout << " + " << number.coef.at(i) << "*"
                 << std::to_string(number.prime_base) + "^" + std::to_string(i + number.val);
    }
    return os;
}

uslong padicNumber::next(uslong i) {
    if (i - this->prec > 1)
        throw invalid_argument("Padic number next should get i - this->prec >= 1");
    if (i > this->prec) {
        this->prec++;
        return computeOX();
    } else return coef.at(i);
}

uslong padicNumber::computeOX() {
    uslong d;
    d = reduce(this->Ox, prime_base);
    if (is_negative) {
        coef.push_back(prime_base - d - 1);
    } else {
        coef.push_back(d);
    }
    return coef.back();
}

padicSum::padicSum(padicRepresentation &op1, padicRepresentation &op2) : op1(op1), op2(op2) {
    if (op1.prime_base != op2.prime_base)
        throw invalid_argument("Bases in sum should be identical");
    val_diff = op1.val - op2.val;
    prec = 0;
    this->val = min(op1.val, op2.val);
}

uslong padicSum::computeSum() {
    if (val_diff < 0) {
        val_diff++;
        return op1.next(prec);
    }
    if (val_diff > 0) {
        val_diff--;
        return op2.next(prec);
    }
    if (val_diff == 0) {
        uslong base = op2.prime_base;
        uslong t = op1.next(prec) + op2.next(prec) + excess;
        if (t < base) {
            excess = 0;
            return t;
        } else {
            excess = 1;
            return t - base;
        }
    }
    throw logic_error("computeSum val not fir in any if");
}

uslong padicSum::next(uslong i) {
    if (i - this->prec > 1)
        throw invalid_argument("Padic number next should get i - this->prec >= 1");
    if (i > this->prec) {
        this->prec++;
        return computeSum();
    } else return coef.at(i);
}




