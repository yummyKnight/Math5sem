//
// Created by Admin on 25.12.2020.
//

#include "relax_padic.hpp"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>
#include <cassert>
#include <cmath>
#include <algorithm>
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
    //std::cout <<"i  "<< i <<endl;
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

padicOperator::padicOperator(slong excess, padicRepresentation &op1, padicRepresentation &op2) : excess(excess),
                                                                                                 op1(op1), op2(op2) {}

void padicOperator::compute_to_N(long long int N) {
    for (long long int i = 1; i < N + 1; ++i) {
        get(i);
    }
}

void padicOperator::compute_to_max() {// +1 to overflow
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
//    if (!is_prime(prime_base) and prime_base == 1)
//        throw std::invalid_argument("P should be prime");

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

padicMul::padicMul(padicRepresentation &op1, padicRepresentation &op2) : padicOperator(0, op1, op2) {
    if (op1.prime_base != op2.prime_base)
        throw invalid_argument("Bases in Mul should be identical");
//  because start from 0
    // TODO: val == 0 ONLY for tests
    this->val = 0;
    // TODO: prec = val
    this->prec = 0;
    this->prime_base = op2.prime_base;
    max_prec = op1.prec + op2.prec - 1;
}

slong padicMul::computeMul() {
    // НЕ готов
    uslong base = op2.prime_base;
    long long t =0;// op1.next(prec) - op2.next(prec) - excess;

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

slong getLn(slong n) {

    slong l = 0, mul_l = 2, n_2 = n + 2;
    while (n_2 % mul_l == 0 && n_2 >= mul_l) {
        if (mul_l == n_2)
            break;
        l++;
        mul_l *= 2;
    }
    return l;
}

slong padicMul::next() {
    slong n = prec-1;
    slong ta = 0, tb = 0;
    //std::cout << "     " <<n <<"      "<< std::endl;
    // начало алгоритма умножения
    // дополнения векторов ya и yb
    ;

    slong l2n = getLn(n * 2), ln = getLn(n), l2n1 = getLn(n * 2 + 1);
    std::vector<slong> ya_new(l2n + 1, 0), yb_new( l2n + 1, 0);
    std::vector<slong> ya1_new( l2n1 + 1, 0), yb1_new( l2n1 + 1, 0);
    ya.push_back(ya_new);
    ya.push_back(ya1_new);
    yb.push_back(yb_new);
    yb.push_back(yb1_new);

    op1.get(n+1);
    op2.get(n+1);
    // вычисление алгоритмов
    slong coef_a, coef_b, left, right;
    slong pow_prime = 1;

    for (slong q = 0; q <= ln; q++) {
        //std::cout << "q " << q <<" "<< std::endl;
        slong q_2 = pow(2, q);
        slong  k = (n + 2) / q_2;
        ta += ya.at(n).at(q);
        tb += yb.at(n).at(q);
        slong coef_a = 0, coef_b = 0, left = q_2, right = q_2 * 2-1;
        pow_prime = 1;
        for (slong i = left-1, i_k = left*(k-1)-1 ; i < right; i++,i_k++) {
            // compute to N
            //std::cout << "+ " << i <<" "<< i * (k-1) - 1<< std::endl;
            coef_a += op1.get(i+1) * pow_prime;
            coef_b += op2.get(i_k+1) * pow_prime;
            //assert(n != 4 || q !=1 || i != 3) ;
            pow_prime *= this->prime_base;
        }

        ta += coef_a * coef_b;
        if (k == 2) {break;}

        coef_a = 0;
        coef_b = 0;
        pow_prime = 1;

        for (slong i = left-1, i_k = left*(k-1)-1 ; i < right; i++,i_k++) {
            coef_a += op1.get( i_k+1 ) * pow_prime;
            coef_b += op2.get( i +1) * pow_prime;
            pow_prime *= this->prime_base;
        }

        tb += coef_a * coef_b;

    }
    //assert(n != 4) ;
    right = pow(2, ln + 1);
    slong sa = ta;

    pow_prime = 1;
    if (this->coef.size() < n + right+1)
        this->coef.resize(n + right, 0);
    //std::cout << coef.size();
    for (slong i = n; i < n + right; i++) {
        sa += this->coef.at(i- val) * pow_prime;
        pow_prime *= this->prime_base;
    }

    pow_prime = 1;
    for (slong i = 0; i < right; i++) {
        this->coef.at(n + i- val) = sa / pow_prime % this->prime_base;
        pow_prime *= this->prime_base;
    }
    if (n + 2 != right)
        ya.at(n + right).at(ln) = sa / pow_prime;

    slong sb = tb;
    pow_prime = 1;
    right = pow(2, ln + 1);

    for (slong i = n; i < n + right; i++) {
        sb += this->coef.at(i- val) * pow_prime;
        pow_prime *= this->prime_base;
    }
    pow_prime = 1;
    for (slong i = 0; i < right; i++) {
        this->coef.at(n + i - val) = sb / pow_prime % this->prime_base;
        pow_prime *= this->prime_base;
    }

    if (n + 2 != right) {
        yb.at(n + right ).at(ln) = sb / pow_prime;
    }

    //std::cout << "sa sb       " << sa <<" "<< sb << " " << this->coef.at(n-val)<< std::endl;

    if (coef.at(n-val) == 0 && n == this->val) {
        //std::cout << "coef  n  val  " << coef.at(n-val) <<" "<< n << " " << val<< std::endl;
        this->val++;
        coef.erase(coef.begin() );
    }
    return this->coef.at(n);
}