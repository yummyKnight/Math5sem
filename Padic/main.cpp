#include <iostream>
#include "relax_padic.cpp"

int main() {
    padicNumber test = padicNumber(-958, 5);
    cout << test << endl;
//    padicNumber test2 = padicNumber(5, 2);
//    test2.next(3);
//    cout << test2 << endl;
//    padicNumber test3 = padicNumber(-1, 5);
//    test2.next(3);
//    cout << test3 << endl;
//    padicNumber test2_1 = padicNumber(375 ,5);
//    padicNumber test2_2 = padicNumber(125, 5);
//    padicSub tmp_sum = padicSub(test2_1, test2_2);
//    tmp_sum.compute_to_max();
//    cout << tmp_sum << endl;
    padicNumber test2_1 = padicNumber(1000000478781 ,5);
    padicNumber test2_2 = padicNumber(100000478781, 5);
    padicSub tmp_sum = padicSub(test2_1, test2_2);
    // 3*5^0 + 1*5^1 + 3*5^2 + 2*5^3 + 1*5^4
    tmp_sum.compute_to_max();
    cout << tmp_sum.to10base() << endl;
    cout << tmp_sum << endl;
    cout << -6 / (slong)(uslong)5 << endl;
}