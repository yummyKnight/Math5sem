//
// Created by Admin on 25.12.2020.
//

#include "relax_padic.hpp"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>
#include <cmath>
#include <cassert>

using namespace std;


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

uslong get_invertible(uslong scalar, uslong prime_base) {
    for (uslong i = 0; i < prime_base; ++i) {
        if ((scalar * i) % prime_base == 1)
            return i;
    }
    return 0;
}

uslong converToRing(slong scalar, uslong prime_base) {
    uslong res;
    if (scalar < 0) {
        slong base = prime_base;
        scalar = scalar % base;
        if (scalar < 0) res = scalar + prime_base;
        else res = 0;
    } else if (scalar >= prime_base)
        res = scalar % prime_base;
    else
        res = scalar;
    return res;
}

slong padicRepresentation::to10base() {
    slong base10 = 0;
//    if negative
//    slong prime = prime_base;
//    base10 += abs(coef.at(0) - prime) * pow(prime_base, val);
//    for (uslong i = 1; i < coef.size(); ++i) {
//        base10 += abs(coef.at(i) - prime + 1) * pow(prime_base, i + val);
//    }
    for (uslong i = 0; i < coef.size(); ++i) {
        base10 += coef.at(i) * pow(prime_base, i + val);
    }
    return base10;
}

long long padicRepresentation::get(long long i) {
    assert(i > 0);
    if (i - this->prec > 1)
        throw invalid_argument("Padic number next should get i - this->prec >= 1");
    if (i > this->prec || this->prec == 0) {
        this->prec++;
        return next();
    } else {
        if (i < val + 1)
            return 0;
        else
            return coef.at(i - 1 - val);
    }
}

void padicRepresentation::increase_val() {
    val++;
}

padicOperator::padicOperator(slong excess, padicRepresentation &op1, padicRepresentation &op2) : excess(excess),
                                                                                                 op1(op1), op2(op2) {}

void padicOperator::compute_to_N(long long int N) {
    for (long long int i = 1; i < N + 1; ++i) {
        get(i);
    }
}

void padicOperator::compute_to_max() { // +1 to overflow
    assert(max_prec > 0);
    compute_to_N(max_prec);
}

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

slong padicNumber::next() {
    slong d;
    d = reduce(this->Ox, prime_base);
    if (is_negative) {
        coef.push_back(prime_base - d - 1);
    } else {
        coef.push_back(d);
    }
    return coef.back();
}

padicSum::padicSum(padicRepresentation &op1, padicRepresentation &op2) : padicOperator(0, op1, op2) {
    if (op1.prime_base != op2.prime_base)
        throw invalid_argument("Bases in sum should be identical");
//  because start from 0
    this->val = min(op1.val, op2.val);
    this->prec = val;
    this->prime_base = op2.prime_base;
    max_prec = max(op1.prec, op2.prec) + 1; // 1 for overflow
}

slong padicSum::next() {
    uslong base = op2.prime_base;
    slong t = op1.get(prec) + op2.get(prec) + excess;
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

padicSub::padicSub(padicRepresentation &op1, padicRepresentation &op2) : padicOperator(0, op1, op2) {
    if (op1.prime_base != op2.prime_base)
        throw invalid_argument("Bases in sum should be identical");
//  because start from 0
    this->val = min(op1.val, op2.val);
    this->prec = val;
    this->prime_base = op2.prime_base;
    max_prec = max(op1.prec, op2.prec);
}

slong padicSub::next() {
    uslong base = op2.prime_base;
    long long t = op1.get(prec) - op2.get(prec) - excess;
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
}

/*
     if (scalar < 0) {
        scalar = scalar / (slong) prime_base;
        if (scalar < 0) u_scalar = scalar + prime_base;
        else u_scalar = 0;
    }
*/
mulRemPadic::mulRemPadic(padicRepresentation &op1, uslong uScalar) : op1(op1), u_scalar(uScalar) {}

slong mulRemPadic::do_job(uslong prec) {
    return (op1.get(prec) * u_scalar) % op1.prime_base;
}

mulQuoPadic::mulQuoPadic(padicRepresentation &op1, uslong uScalar) : op1(op1), u_scalar(uScalar) {}

slong mulQuoPadic::do_job() {
    return (op1.coef.back() * u_scalar) / op1.prime_base;
}

scalarDivPadic::scalarDivPadic(padicRepresentation &op1, slong scalar) : op1(op1) {
    prec = op1.val;
    val = op1.val;
    prime_base = op1.prime_base;
    coef = vector<slong>();
    if (scalar < 0) {
        scalar = scalar / (slong) prime_base;
        if (scalar < 0) u_scalar = scalar + prime_base;
        else u_scalar = 0;
    } else if (scalar >= prime_base)
        u_scalar = scalar % prime_base;
    else u_scalar = scalar;

    scalar_inv = get_invertible(scalar, prime_base);
    if (scalar_inv == 0)
        throw invalid_argument("Scalar should be invertible in Z_p");


//  push 0 element (precision still == prec from op1)


}


long long scalarDivPadic::next() {
    uslong c;
    if (coef.empty()) { // if 1st time
        c = (op1.get(prec) * scalar_inv) % prime_base;
        coef.push_back(c);
        return c;
    } else {
        slong mul_quo = (coef.back() * u_scalar) / prime_base;
        slong part = op1.get(prec) - mul_quo;
        uslong u_part = converToRing(part, prime_base);
        c = (u_part * scalar_inv) % prime_base; // mul_rem
        if (c == 0 && prec - 1 - val == 0)
            val++;
        else
            coef.push_back(c);
        return c;
    }
}

void scalarDivPadic::compute_to_max() {
    for (slong i = 1; i < op1.prec + 1; ++i) {
        get(i);
    }
}

void scalarDivPadic::compute_to_N(slong N) {
    for (slong i = 1; i < N; ++i) {
        get(i);
    }
}

long long int divPadic::next() {
    uslong c;
    uslong scalar_inv;
    if (coef.empty()) {
        c = u * op1.get(op1.val + 1) % prime_base;
        coef.push_back(c);
        return c;
    } else {
        slong part_a = op1.get(prec);
        slong part_b = op2.get(prec) - b0;
        slong part_c = part_b * coef.back();
        slong part_d = converToRing(part_a - part_c, prime_base);
        c = part_d / b0;
        coef.push_back(c);
        return c;
    }
}

divPadic::divPadic(padicRepresentation &op1, padicRepresentation &op2) : padicOperator(0, op1, op2) {

    if (op1.prime_base != op2.prime_base)
        throw invalid_argument("Bases in div should be identical");
//  because start from 0
    this->val = max(op1.val, op2.val);
    this->prec = val;
    this->prime_base = op2.prime_base;
    max_prec = max(op1.prec, op2.prec);
    b0 = op2.get(val + 1);
    u = get_invertible(get_invertible(b0, prime_base), prime_base);
    if (u == 0)
        throw invalid_argument("Scalar should be invertible in Z_p");
}
