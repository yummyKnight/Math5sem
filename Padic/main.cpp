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
//    padicNumber test2_1 = padicNumber(1000000478781 ,5);
//    padicNumber test2_2 = padicNumber(100000478781, 5);
//    padicSub tmp_sum = padicSub(test2_1, test2_2);
    // 3*5^0 + 1*5^1 + 3*5^2 + 2*5^3 + 1*5^4

//    padicNumber test2_1 = padicNumber(5 ,7);
//    padicNumber test2_2 = padicNumber(13, 7);
//    padicNumber test2_3 = padicNumber(65, 7);

//    padicNumber test2_1 = padicNumber(425 ,7);
//    padicNumber test2_2 = padicNumber(739, 7);
//    padicNumber test2_3 = padicNumber(314075, 7);

//    padicNumber test2_1 = padicNumber(7 ,7);
//    padicNumber test2_2 = padicNumber(49, 7);
//    padicNumber test2_3 = padicNumber(343, 7);

    padicNumber test2_1 = padicNumber(-100500 ,7);
    padicNumber test2_2 = padicNumber(492400, 7);
    padicNumber test2_3 = padicNumber(494912250000, 7);

    cout << test2_1 << endl;
    cout << test2_2 << endl;
    cout << test2_3 << endl;

    //padicMul tmp_sub = padicMul(test2_1, test2_2);
    //padicSub tmp_sub = padicSub(test2_1, test2_2);
    padicSum tmp_sum = padicSum(test2_1, test2_2);
    //tmp_sub.compute_to_max();
    //cout <<"result\n "<< tmp_su << endl;
    cout <<"result\n "<< tmp_sum.is_negative <<"   "<< test2_2.is_negative << endl;
    cout << endl;

//    for (const auto& obj : test2_1.coef)
//        std::cout << ' ' << obj;
//    cout << endl;
//    for (const auto& obj : test2_2.coef)
//        std::cout << ' ' << obj;
//    cout << endl;
//    for (const auto& obj : test2_3.coef)
//        std::cout << ' ' << obj;
//    cout << endl;
//    for (const auto& obj : tmp_sum.coef)
//        std::cout << ' ' << obj;
//    cout << endl << tmp_sum.val;
//    cout << endl << test2_1.val;
    //cout << tmp_sum << endl;
    //cout << -6 / (slong)(uslong)5 << endl;
}