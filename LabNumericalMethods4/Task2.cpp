#include <iostream>
#include "Methods/Eiler.h"
#include <functional>
#include <ctime>
#include "Methods/RungeKutta4.h"
int main() {
    double sigma=10, r=13,b=8/3.;
    const std::vector<std::string> var{"x", "y","z"};
    std::vector<std::string> name{"EilerSys", "RK4Sys"};
    std::function<std::vector<double>(double, std::vector<double>)> func = [sigma, r, b
            ](double t, std::vector<double> var)-> std::vector<double> {
        return { sigma*(var[1]-var[0]), var[0]*(-var[2]+r)-var[1], var[0]*var[1]-b*var[2]};
    };
    std::pair<double, std::vector<double>> init=std::make_pair(0, std::vector<double>{0,1,0});
    double tau0=1e-2;
    double taurk4=1e-2;
    double t0=0,tn=10;
    std::vector<std::pair<double, std::vector<double>>> answer;
    std::clock_t start_time = std::clock();
    EilerSys<double>::answer(func, init,tau0, t0, tn);
    std::clock_t end_time = std::clock();
    std::cout<<(end_time-start_time)<<"\n";
    start_time=std::clock();
    answer = RK4Sys<double>::answer_with_err(func, init,taurk4, t0, tn);
    end_time=std::clock();
    std::cout<<(end_time-start_time)<<"\n";
    for (int i=0; i<name.size(); i++) {
        for (int j=2; j<=4;++j) {
            std::ofstream script("temp_plot.gp", std::ios::out);
            script << "set terminal pngcairo size 1200,1200 enhanced font 'Arial,12'\n";
            script << "set output '../Graphics/"+name[i]+'_'+var[j-2]+"(t).png'\n";
            script << "set title 'Решение системы ОДУ методом Эйлера'\n";
            script << "set xlabel 'Время t'\n";
            script << "set ylabel \'"+ var[j-2]<<"\'\n";
            script << "set grid\n";
            script<<"set xrange ["<<t0<<":"<<tn<<"]\n";
            script << "plot \'../"+name[i]+".txt\' using 1:"+std::to_string(j)+" with lines title '"+ var[j-2]+"(t)', \\\n";
            script.close();
            std::system("gnuplot -persist temp_plot.gp");
        }
    }
    for (int i=1; i<name.size(); i++) {
        for (int phi=0; phi<360; phi+=1){
        std::ofstream script("temp_plot_3d.gp", std::ios::out);

        script << "set terminal pngcairo size 1200,1200 enhanced font 'Arial,12'\n";
        script << "set output '../Graphics/GIF/Phase3D_" << name[i] <<"phi="<<phi<< ".png'\n";
        script << "set title '3D фазовая траектория системы ОДУ (RK4)'\n";
        script << "set xlabel 'x'\n";
        script << "set ylabel 'y'\n";
        script << "set zlabel 'z'\n";
        script << "set grid\n";
        script << "set view 70,"<<phi<<"\n";
        script << "splot '../" << name[i]<<".txt' using 2:3:4 with points lc rgb 'red' title 'траектория'\n";
        script.close();
        std::system("gnuplot temp_plot_3d.gp");
    }
}
    return 0;
}