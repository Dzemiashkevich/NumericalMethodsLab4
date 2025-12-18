#include <iostream>
#include "Methods/Eiler.h"
#include <functional>
#include <ctime>
#include "Methods/RungeKutta4.h"
// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
int main() {
    double omega=1;
    const std::vector<std::string> var{"dx", "x"};
    std::vector<std::string> name{"EilerSys", "RK4Sys"};
    std::function<std::vector<double>(double, std::vector<double>)> func = [omega
            ](double t, std::vector<double> var)-> std::vector<double> {
        return { -omega * omega * var[1], var[0]};
    };
    std::pair<double, std::vector<double>> init {std::make_pair(0, std::vector<double>{0,1})};
    double tau0=1e-4;
    double taurk4=1e-3;
    double t0=0,tn=20;
    std::clock_t start_time = std::clock();
    EilerSys<double>::answer(func, init,tau0, t0, tn);
    std::clock_t end_time = std::clock();
    std::cout<<(end_time-start_time)<<"\n";
    start_time=std::clock();
    RK4Sys<double>::answer_with_err(func, init, t0, tn, 1e-8);
    end_time=std::clock();
    std::cout<<(end_time-start_time)<<"\n";
    std::ofstream AnalyticRK4;
    AnalyticRK4.open("../Analytic_X(t)_RK4.txt",std::ios::out);
    for (int i=0; i<static_cast<int>((-t0 +tn) / taurk4);++i) {
        AnalyticRK4<<t0+taurk4*i<<"\t"<<-omega*std::sin(omega*(t0+i*taurk4))<<"\t"<<cos(omega*(t0+i*taurk4))<<"\n";
    }
    AnalyticRK4.close();
    for (int i=0; i<name.size(); i++) {
        for (int j=2; j<=3;++j) {
            std::ofstream script("temp_plot.gp", std::ios::out);
            script << "set terminal pngcairo size 1200,1200 enhanced font 'Arial,12'\n";
            script << "set output '../Graphics/"+name[i]+'_'+var[j-2]+"(t).png'\n";
            script << "set title 'Решение системы ОДУ методом Эйлера'\n";
            script << "set xlabel 'Время t'\n";
            script << "set ylabel \'"+ var[j-2]<<"\'\n";
            script << "set grid\n";
            script<<"set xrange ["<<t0<<":"<<tn<<"]\n";
            script << "plot \'../"+name[i]+".txt\' using 1:"+std::to_string(j)+" with lines title '"+ var[j-2]+"(t)', \\\n";
            script << " \'../Analytic_X(t)_RK4.txt\' using 1:"+std::to_string(j)+" with lines title '"+ var[j-2]+"(t)', \\\n";
            //script << " \'../Analytics.txt\' using 1:"+std::to_string(j)+" with lines title 'A"+ var[j-2]+"(t)', \\\n";
            script.close();
            std::system("gnuplot -persist temp_plot.gp");
        }
    }
    return 0;
}