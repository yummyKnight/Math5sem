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

bool is_prime(slong n) {
    for (int i = 2; i * i <= n; ++i) //no sqrt, please
    {
        if (n % i == 0) return false;
    }
    return true;
}

slong reduce(slong &x, uslong prime_base) {
    slong d;
    d = x % prime_base;
    x -= d;
    x = x / prime_base;
    return d;
}

padicNumber::padicNumber(long long base10, uslong prime_base, uslong init_prec) {
    slong x;
    slong d = 0;
    slong j;
    long long val = 0;
    std::vector<slong> coef;
    is_negative = base10 < 0;
    x = abs(base10);
    if (!is_prime(prime_base) and prime_base == 1)
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

    if (is_negative) {
        for (j = 0; j < init_prec - val; j++) {
            d = reduce(x, prime_base);
            coef.push_back(prime_base - d - 1);
        }
    } else {
        for (j = 0; x != 0 && j < init_prec - val; j++) {
            d = reduce(x, prime_base);
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

    if (number.coef.empty())
        os << 0;
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


slong padicNumber::next(long long i) {
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

slong padicNumber::computeOX() {
    slong d;
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

slong padicSum::computeSum() {
    uslong base = op2.prime_base;
    slong t = op1.next(prec) + op2.next(prec) + excess;
    if (t < base) {
        excess = 0;
        coef.push_back(t);
        return t;
    } else {
        excess = 1;
        slong res = t - base;
        // check if leading zero
        if (res == 0 && prec - 1 - val == 0)
            this->val++;
        else
            coef.push_back(res);
        return res;
    }
}

slong padicSum::next(long long i) {
    assert(i > 0);
    if (i - this->prec > 1)
        throw invalid_argument("Padic number next should get i - this->prec >= 1");
    if (i > this->prec || this->prec == 0) {
        this->prec++;
        auto res = computeSum();
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
    slong max_N = max(op1.prec, op2.prec) + 1; // +1 to overflow
    compute_to_N(max_N);
}


slong padicSub::next(long long i) {
    assert(i > 0);
    if (i - this->prec > 1)
        throw invalid_argument("Padic number next should get i - this->prec >= 1");
    if (i > this->prec || this->prec == 0) {
        this->prec++;
        auto res = computeSub();
        return res;
    } else {
        if (i < val + 1)
            return 0;
        else
            return coef.at(i - 1 - val);
    }
}

padicSub::padicSub(padicRepresentation &op1, padicRepresentation &op2) : op1(op1), op2(op2) {
    if (op1.prime_base != op2.prime_base)
        throw invalid_argument("Bases in sum should be identical");
//  because start from 0
    this->val = min(op1.val, op2.val);
    this->prec = val;
    this->prime_base = op2.prime_base;
}

void padicSub::compute_to_N(long long int N) {
    for (long long int i = 0; i < N; ++i) {
        next(i + 1);
    }
}

void padicSub::compute_to_max() {
    slong max_N = max(op1.prec, op2.prec) + 1; // +1 to overflow
    compute_to_N(max_N);
}

slong padicSub::computeSub() {
    uslong base = op2.prime_base;
    long long t = op1.next(prec) - op2.next(prec) - excess;
    if (t >= 0) {
        excess = 0;
        // check if leading zero
        if (t == 0 && prec - 1 - val == 0)
            this->val++;
        else
            coef.push_back(t);
        return t;
    } else {
        excess = 1;
        slong res = t + base;
        coef.push_back(res);
        return res;
    }
    return 0;
}
