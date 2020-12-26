#include <iostream>
#include "relax_padic.cpp"

int main() {
    padicNumber test = padicNumber(375, 5);
    cout << test << endl;
//    padicNumber test2 = padicNumber(5, 2);
//    test2.next(3);
//    cout << test2 << endl;
//    padicNumber test3 = padicNumber(-1, 5);
//    test2.next(3);
//    cout << test3 << endl;
    padicNumber test2_1 = padicNumber(375 ,5);
    padicNumber test2_2 = padicNumber(125, 5);
    padicSub tmp_sum = padicSub(test2_1, test2_2);
    tmp_sum.compute_to_max();
    cout << tmp_sum << endl;
}