#include "ProyectoF.h"
#include <vector>
#include <cstdlib>
#include <fstream>
#include <omp.h>


int main(void)
{
  std::vector<double> F(4);
  std::vector<double> INt(7);
  std::ifstream fin("input.txt");
  fin >> F[0];
  fin >> F[1];
  fin >> F[2];
  fin >> F[3];
  fin >> F[4];
  fin >> INt[0];
  fin >> INt[1];
  fin >> INt[2];
  fin >> INt[3];
  fin >> INt[4];
  fin >> INt[5];
  fin >> INt[6];

  const double DT = 0.04; //0.01*2*M_PI*std::sqrt(l/G);
  const double TF = 100;
  const int N = 2;
  double q=0.1, q1=0.5;           // lectura de archivo.txt
  double omegaf=0.66;





  //int nth = omp_get_num_threads();
  int thid=omp_get_thread_num();


#pragma omp parallel private(thid)  num_threads(10)
  {

    thid=omp_get_thread_num();
    std::vector<double> y(N);
    //nth = omp_get_num_threads();
    if(thid==0){
      initial_conditions(y, 0.5, INt[6]);
      Simulacion(y, 0.0, TF, DT, 9.8, q, F[0], omegaf, 1);//punto 1
    }
    if(thid==1){
  initial_conditions(y, INt[1], INt[6]); // 0.2 theta
  Simulacion(y, 0.0, TF, DT, 9.8, q1, F[1], omegaf, 21);//punto 2
    }
    if(thid==2){
      initial_conditions(y, INt[1], INt[6]);
  Simulacion(y, 0.0, TF, DT, 9.8, q1, F[2], omegaf, 22);// punto 2
      }
    if(thid==3){
      initial_conditions(y, INt[1], INt[6]);
      Simulacion(y, 0.0, TF, DT, 9.8, q1, F[3], omegaf, 23);} //punto2
   if(thid==4){
     initial_conditions(y, INt[1], INt[6]);
     Simulacion(y, 0.0, 500000, DT, 9.8, q1, F[4], omegaf, 3);}//punto3

     if(thid==5){
    initial_conditions(y, 0.2, INt[6]);
    Simulacion1(y, 0.0, 60, DT, 9.8, q1, F[4], 0.666667 ,4); //punto 4
    }
    if(thid==6){
  initial_conditions(y, INt[1], INt[6]);
  Simulacion2(y, 0.0, 60, DT, 9.8, q1, F[4], 0.666667, 5); //punto 5
    }
     if(thid==7){
  initial_conditions(y, INt[2], INt[6]);
  Simulacion(y, 0.0, 10000, DT, 9.8, q1, F[4], omegaf,611); //punto 6
  Simulacion(y, 0.0, 18000, 0.043, 9.8, q1, F[4], omegaf,612); //punto 6 (distinto dt)
  Simulacion(y, 0.0, 10000, 0.045, 9.8, q1, F[4], omegaf,613); //punto 6 (distinto dt)
    }
    if(thid==8){
  initial_conditions(y, INt[3], INt[6]);
  Simulacion(y, 0.0, 12000, DT, 9.8, q1, F[4], omegaf,621); //punto 6
  Simulacion(y, 0.0, 20000, 0.043, 9.8, q1, F[4], omegaf,622); //punto 6
  Simulacion(y, 0.0, 18000, 0.045, 9.8, q1, F[4], omegaf,623); //punto 6
    }
    if(thid==9){
  initial_conditions(y, INt[4], INt[6]);
  Simulacion(y, 0.0, 3000, DT, 9.8, q1, F[4], omegaf,631); //punto 6
  Simulacion(y, 0.0, 15000, 0.043, 9.8, q1, F[4], omegaf,632); //punto 6
  Simulacion(y, 0.0, 18000, 0.045, 9.8, q1, F[4], omegaf,633);} //punto 6

  }


return 0; }
