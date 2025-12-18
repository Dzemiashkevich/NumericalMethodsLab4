//
// Created by Dzemiashkevich Vlad on 15.11.25.
//

#ifndef EILER_H
#define EILER_H
#include <fstream>
#include <functional>
#include <vector>
#include "../../NumericalMethodsLab2.1/Methods/Newton.h"
template <typename T>
class EilerEq {


public:
    std::function<T(T)> f;
    T tau, t0, tn;

    EilerEq(T tau,T t0, T tn,  std::function<T(T)> f) : tau(tau), f(f), t0(t0), tn(tn) {}

};

template <typename T>
class EilerSys:public EilerEq<T> {
    std::function<std::vector<T>(T t, std::vector<T> var)> system;
    public:

    EilerSys(T tau, T t0, T tn, std::function<std::vector<T>(T,std::vector<T>)> sys): EilerEq<T>(tau, t0, tn, nullptr), system(sys) {}

    void solve(const std::pair<T,std::vector<T>>& init) {
            std::vector<T> res = init.second;
            std::vector<std::pair<T, std::vector<T> > > answers;
        std::fstream file;
        file.open("../EilerSys.txt", std::ios::out);
        file<<0;
        for (int i=0; i<res.size(); ++i) {
            file<<"\t";
            file<<init.second[i];
        }
        file<<"\n";
            for (T t = this->t0; t < this->tn; t += this->tau) {
                file<<t+this->tau;
                auto sys = system(t, res);
                for (int i=0;i<res.size();i++) {
                    file<<"\t";
                    res[i]+=this->tau*sys[i];
                    file<<res[i];
                }
                file<<"\n";
            }
        file.close();
    }
    static void answer(std::function<std::vector<T>(T, std::vector<T>)> sys, std::pair<T,std::vector<T>> init,T tau0, T t0=0, T tn=100) {
        EilerSys<T> Ei(tau0, t0, tn, sys);
        Ei.solve(init);
    }
};



#endif //EILER_H
