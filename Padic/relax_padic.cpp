//
// Created by Admin on 25.12.2020.
//

#include "relax_padic.hpp"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>
#include <cassert>

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
    for (j = 0; x != 0 || j > init_prec - val; j++) {
        d = reduce(x, prime_base);
        if (is_negative) {
            coef.push_back(prime_base - d - 1);
        } else {
            coef.push_back(d);
        }
    }
    // prec тот макимум в масиве который уже посчитан
    this->prec = coef.size() + val;
    this->Ox = x;
    this->prime_base = prime_base;
    this->coef = coef;
    this->val = val;
}

ostream &operator<<(ostream &os, const padicRepresentation &number) {
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


uslong padicNumber::next(long long i) {
    assert(i > 0);
    if (i - this->prec > 1)
        throw invalid_argument("Padic number next should get i - this->prec >= 1");
    if (i > this->prec || this->prec == 0) {
        auto res = computeOX();
        this->prec++;
        return res;
    } else {
        if (i < val + 1)
            return 0;
        else
            return coef.at(i - 1 - val);
    }

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
//  because start from 0
    this->val = min(op1.val, op2.val);
    this->prec = val;
    this->prime_base = op2.prime_base;
}

uslong padicSum::computeSum() {
    uslong base = op2.prime_base;
    uslong t = op1.next(prec + 1) + op2.next(prec + 1) + excess;
    if (t < base) {
        excess = 0;
        coef.push_back(t);
        return t;
    } else {
        excess = 1;
        uslong res = t - base;
        coef.push_back(res);
        return res;
    }
}

uslong padicSum::next(long long i) {
    assert(i > 0);
    if (i - this->prec > 1)
        throw invalid_argument("Padic number next should get i - this->prec >= 1");
    if (i > this->prec || this->prec == 0) {
        auto res = computeSum();
        this->prec++;
        return res;
    } else {
        if (i < val + 1)
            return 0;
        else
            return coef.at(i - 1 - val);
    }
}

void padicSum::compute_to_N(long long int N) {
    for (long long int i = 0; i < N; ++i) {
        next(i + 1);
    }
}

void padicSum::compute_to_max() {
    uslong max_N = max(op1.prec, op2.prec);
    compute_to_N(max_N);
}





