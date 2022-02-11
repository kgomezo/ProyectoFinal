#ifndef PROYECTO_H_
#define PROYECTO_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>

void initial_conditions(std::vector<double> & y, double inithe, double initomega);
void fderiv(const std::vector<double> & y, std::vector<double> & dydt, double t,double l, double q, double F, double omegaf);

void Simulacion(std::vector<double>  y, double tinit, double tend, double dt,  double l, double q, double F, double omegaf, int id);

void Simulacion1(std::vector<double>  y, double tinit, double tend, double dt, double l, double q, double F, double omegaf,int id);

void Simulacion2(std::vector<double>  y, double tinit, double tend, double dt, double l, double q, double F, double omegaf,int id);

std::vector<double> Simulacion_test(std::vector<double>  y, double tinit, double tend, double dt, double l, double q, double F, double omegaf);

#endif
