//
// Created by Dzemiashkevich Vlad on 4.12.25.
//

#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H
#include <vector>
#include <functional>
#include <fstream>

template <typename T>
class RK4 {


public:
    std::function<T(T)> f;

    T tau, t0, tn, tau2;

    RK4(T tau,T t0, T tn,  std::function<T(T)> f) : tau(tau), f(f), t0(t0), tn(tn) {}

};

template <typename T>
class RK4Sys:public RK4<T> {
    std::function<std::vector<T>(T t, std::vector<T> var)> system;
public:

    RK4Sys(T tau, T t0, T tn, std::function<std::vector<T>(T,std::vector<T>)> sys): RK4<T>(tau, t0, tn, nullptr), system(sys) {}

    std::pair<T, std::vector<T>> stepRK4(const std::pair<T, std::vector<T>>& init, T tau){
        int n=init.second.size();
        std::vector<T> k1(n),k2(n),k3(n),k4(n);
         std::vector<T> res=init.second;
         std::vector<T> tmp(res.size());
        T t=init.first;
             k1 = system(t, res);
            for (int i=0; i<res.size(); ++i) {
                tmp[i]=res[i]+tau*k1[i]/2;
            }
             k2=system(t+tau/2, tmp);
            for (int i=0; i<res.size(); ++i) {
                tmp[i]=res[i]+tau*k2[i]/2;
            }
             k3=system(t+tau/2, tmp);
            for (int i=0; i<res.size(); ++i) {
                tmp[i]=res[i]+tau*k3[i];
            }
            k4=system(t+tau, tmp);
            for (int i=0;i<res.size();i++) {
                res[i]+=tau*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
            }
            return std::make_pair(t+tau, res);
    }

    std::vector<std::pair<double, std::vector<double>>> solve(const std::pair<T, std::vector<T>>& init) {
        std::pair<T, std::vector<T>> res = init;
        std::vector<std::pair<T, std::vector<T>>> answers;
        answers.push_back(res);
            for (T t = this->t0; t < this->tn; t += this->tau) {
                res=stepRK4( res, this->tau);
                answers.pushback(res);
            }
            return answers;
    }

    std::vector<std::pair<T, std::vector<T>>> err_solve(const std::pair<T, std::vector<T>>& init, T err) {
        for (long int i=1;i<100000 ; i++) {

            bool exit=false;
            std::pair<T, std::vector<T>> res0 = init;
            std::pair<T, std::vector<T>> res1 = init;

            std::vector<std::pair<T, std::vector<T>>> answers0 = {res0};
            std::vector<std::pair<T, std::vector<T>>> answers1 = {res1};

            T tau0=static_cast<T>((this->tn-this->t0)/(i));
            T tau1=tau0/2.;
            for (int step = 0; step < i; ++step) {

                res0=stepRK4(res0, tau0);
                answers0.push_back(res0);

                for (int j=step; j<step+2; ++j) {
                    res1=stepRK4( res1, tau1);
                    answers1.push_back(res1); // 0 0.5 1   0 1
                }

                /*                                    //Проверка работы данной части, не знаю почему, но при отладке нельзя посмотреть что зранится в векторе
                 std::cout<<"\n"<<i<<"\n";
                 std::cout<<res1.first<<"\t";
                for (int i=0; i<res1.second.size(); ++i) {
                    std::cout<<res1.second[i]<<"\t";
                }
                std::cout <<"\n"<<"---------"<<"\n";

                std::cout<<res0.first<<"\t";
                for (int i=0; i<res1.second.size(); ++i) {
                    std::cout<<res0.second[i]<<"\t";
                }
                std::cout <<"\n"<<"---------"<<"\n"<<"\n";*/

                for (int j=0; j<init.second.size(); ++j) {
                    if (abs(res1.second[j]-res0.second[j])>err) {
                        exit=true;
                        break;
                    }
                }
                if (exit) break;
            }
            if (!exit) {
                this->tau=tau0;
                return answers0;
            }
            exit=false;
        }
        throw std::invalid_argument("Can't solve RK4Sys with this tolerance");
    }



    static std::vector<std::pair<double, std::vector<double>>> answer(std::function<std::vector<T>(T, std::vector<T>)> sys, std::pair<T,std::vector<T>> init,T tau0, T t0=0, T tn=100) {
        RK4Sys<T> RK4(tau0, t0, tn, sys);
        std::vector<std::pair<double, std::vector<double>>> answer = RK4.solve(init);
        std::fstream FileRK4;
        FileRK4.open("../RK4Sys.txt",std::ios::out);
        for (auto & i : answer) {
            FileRK4<<i.first<<"\t";
            for (double j : i.second) {
                FileRK4<<j<<"\t";
            }
            FileRK4<<"\n";
        }
        FileRK4.close();
         return answer;
    }
    static std::vector<std::pair<double, std::vector<double>>> answer_with_err(std::function<std::vector<T>(T, std::vector<T>)> sys, std::pair<T,std::vector<T>> init, T t0=0, T tn=100, T err=1e-8) {
        RK4Sys<T> RK4(1, t0, tn, sys);
        std::vector<std::pair<double, std::vector<double>>> answer=  RK4.err_solve(init, err);
        std::fstream FileRK4;
        FileRK4.open("../RK4Sys.txt",std::ios::out);
        for (auto & i : answer) {
            FileRK4<<i.first<<"\t";
            for (double j : i.second) {
                FileRK4<<j<<"\t";
            }
            FileRK4<<"\n";
        }
        FileRK4.close();
        return answer;
    }
};



#endif //RUNGEKUTTA4_H
