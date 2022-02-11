#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>

const double G = 9.8;
void initial_conditions(std::vector<double> & y, double inithe, double initomega);
void fderiv(const std::vector<double> & y, std::vector<double> & dydt, double t,double l, double q, double F, double omegaf);

void Simulacion(std::vector<double>  y, double tinit, double tend, double dt,  double l, double q, double F, double omegaf, int id);

void Simulacion1(std::vector<double>  y, double tinit, double tend, double dt, double l, double q, double F, double omegaf,int id);

void Simulacion2(std::vector<double>  y, double tinit, double tend, double dt, double l, double q, double F, double omegaf,int id);

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
  double l=9.8;
  const double DT = 0.04; //0.01*2*M_PI*std::sqrt(l/G);
  const double TF = 100;
  const int N = 2;
  double q=0.1, q1=0.5;           // lectura de archivo.txt
  double omegaf=0.66;
  std::vector<double> y(N);
  initial_conditions(y, INt[0], INt[6]);
  Simulacion(y, 0.0, TF, DT, 9.8, q, F[0], omegaf, 1); //punto 1

  initial_conditions(y, INt[1], INt[6]); // 0.2 theta
  Simulacion(y, 0.0, TF, DT, 9.8, q1, F[1], omegaf, 21);//punto 2
  Simulacion(y, 0.0, TF, DT, 9.8, q1, F[2], omegaf, 22);// punto 2
  Simulacion(y, 0.0, TF, DT, 9.8, q1, F[3], omegaf, 23);//punto2

  Simulacion(y, 0.0, 500000, DT, 9.8, q1, F[4], omegaf, 3);//punto3

  Simulacion1(y, 0.0, 150, DT, 9.8, q1, F[4], omegaf,4); //punto 4

  Simulacion2(y, 0.0, 60, DT, 9.8, q1, F[4], omegaf,5); //punto 5


  initial_conditions(y, INt[2], INt[6]);
  Simulacion(y, 0.0, 1200, DT, 9.8, q1, F[4], omegaf,611); //punto 6
  Simulacion(y, 0.0, 12000, 0.043, 9.8, q1, F[4], omegaf,612); //punto 6 (distinto dt)
  Simulacion(y, 0.0, 1200, 0.045, 9.8, q1, F[4], omegaf,613); //punto 6 (distinto dt)
  initial_conditions(y, INt[3], INt[6]);
  Simulacion(y, 0.0, 1200, DT, 9.8, q1, F[4], omegaf,621); //punto 6
  Simulacion(y, 0.0, 12000, 0.043, 9.8, q1, F[4], omegaf,622); //punto 6
  Simulacion(y, 0.0, 12000, 0.045, 9.8, q1, F[4], omegaf,623); //punto 6
  initial_conditions(y, INt[4], INt[6]);
  Simulacion(y, 0.0, 3000, DT, 9.8, q1, F[4], omegaf,631); //punto 6
  Simulacion(y, 0.0, 12000, 0.043, 9.8, q1, F[4], omegaf,632); //punto 6
  Simulacion(y, 0.0, 12000, 0.045, 9.8, q1, F[4], omegaf,633); //punto 6

  return 0;
}

void initial_conditions(std::vector<double> & y, double inithe, double initomega)
{
  y[0] = inithe;
  y[1] = initomega;
}



void fderiv(const std::vector<double> & y, std::vector<double> & dydt, double t, double l, double q, double F, double omegaf)
{
  dydt[0] = y[1];//omega
  dydt[1] = -(9.8/l)*std::sin(y[0])-q*y[1]+F*std::sin(omegaf*t); //-(g/l)sin(theta)-q*dtheta/dt + F*sin(omega t)=dW/dt
}



void Simulacion(std::vector<double>  y, double tinit, double tend, double dt, double l, double q, double F, double omegaf, int id)

{

    std::string fname = "Simulacion" + std::to_string(id) + ".txt";
    std::ofstream fout(fname);
    int N = y.size();
    double teo = 0.0;
    std::vector<double> dydt(N, 0.0);
    std::vector<double> k1(N), k2(N), k3(N), k4(N), aux(N);
    double time = 0.0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
    time = 0.0 + ii*dt;
    // k1
    fderiv(y, dydt, time, l, q, F,omegaf );
    for(int ii = 0; ii < N; ++ii) {
      k1[ii] = dt*dydt[ii];
    }
    // k2 aux
    for(int ii = 0; ii < N; ++ii) {
      aux[ii] = y[ii] + k1[ii]/2;
    }
    //k2
    fderiv(aux, dydt, time + dt/2, l, q, F, omegaf);
    for(int ii = 0; ii < N; ++ii) {
      k2[ii] = dt*dydt[ii];
    }
    // k3 aux
    for(int ii = 0; ii < N; ++ii) {
      aux[ii] = y[ii] + k2[ii]/2;
    }
    //k3
    fderiv(aux, dydt, time + dt/2, l, q, F,omegaf);
    for(int ii = 0; ii < N; ++ii) {
      k3[ii] = dt*dydt[ii];
    }
    // k4 aux
    for(int ii = 0; ii < N; ++ii) {
      aux[ii] = y[ii] + k3[ii];
    }
    //k4
    fderiv(aux, dydt, time + dt, l, q, F, omegaf);
    for(int ii = 0; ii < N; ++ii) {
      k4[ii] = dt*dydt[ii];
    }

    // write new data
    for(int ii = 0; ii < N; ++ii) {
      y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
    }


    fout << time << "\t" << y[0] << "\t" << y[1] <<std::endl;

  }
  fout.close();
}

void Simulacion1(std::vector<double>  y, double tinit, double tend, double dt, double l, double q, double F, double omegaf,int id)
{
  std::string fname = "Simulacion" + std::to_string(id) + ".txt";
  std::ofstream fout(fname);
  int N = y.size();
 std::vector<double> w(N);
  initial_conditions(w, 0.201, 0);
 std::vector<double> dydt(N, 0.0);
 std::vector<double> dydt1(N, 0.0);
 std::vector<double> k1(N), k2(N), k3(N), k4(N), aux(N), h1(N), h2(N), h3(N), h4(N), auxi(N);

  double time = 0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
    time = 0.0 + ii*dt;
    // k1
    fderiv(y, dydt, time, l, q, F,omegaf );
    fderiv(w, dydt1, time, l, q, F,omegaf );
    for(int ii = 0; ii < N; ++ii) {
      k1[ii] = dt*dydt[ii];
      h1[ii] = dt*dydt1[ii];
    }
    // k2 aux
    for(int ii = 0; ii < N; ++ii) {
      aux[ii] = y[ii] + k1[ii]/2;
      auxi[ii] = w[ii] + h1[ii]/2;
    }
    //k2
    fderiv(aux, dydt, time + dt/2, l, q, F, omegaf);
    fderiv(auxi, dydt1, time + dt/2, l, q, F, omegaf);
    for(int ii = 0; ii < N; ++ii) {
      k2[ii] = dt*dydt[ii];
      h2[ii] = dt*dydt1[ii];
    }
    // k3 aux
    for(int ii = 0; ii < N; ++ii) {
      aux[ii] = y[ii] + k2[ii]/2;
      auxi[ii] = w[ii] + h2[ii]/2;
    }
    //k3
    fderiv(aux, dydt, time + dt/2, l, q, F,omegaf);
    fderiv(auxi, dydt1, time + dt/2, l, q, F, omegaf);
    for(int ii = 0; ii < N; ++ii) {
      k3[ii] = dt*dydt[ii];
      h3[ii] = dt*dydt1[ii];
    }
    // k4 aux
    for(int ii = 0; ii < N; ++ii) {
      aux[ii] = y[ii] + k3[ii];
      auxi[ii] = w[ii] + h3[ii];
    }
    //k4
    fderiv(aux, dydt, time + dt, l, q, F, omegaf);
    fderiv(auxi, dydt1, time + dt/2, l, q, F, omegaf);
    for(int ii = 0; ii < N; ++ii) {
      k4[ii] = dt*dydt[ii];
      h4[ii] = dt*dydt1[ii];
    }

    // write new data
    for(int ii = 0; ii < N; ++ii) {
      y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
      w[ii] = w[ii] + (h1[ii] + 2*h2[ii] + 2*h3[ii] + h4[ii])/6.0;
    }

     fout << time << "\t" <<std::abs( w[0]-y[0]) << std::endl;

  }
  fout.close();
}

void Simulacion2(std::vector<double>  y, double tinit, double tend, double dt, double l, double q, double F, double omegaf,int id)
{
  std::string fname = "Simulacion" + std::to_string(id) + ".txt";
  std::ofstream fout(fname);
  int N = y.size();
  std::vector<double> w(N);
  initial_conditions(w, 0.201, 0);
  double q1=0.101;
  std::vector<double> dydt(N, 0.0);
  std::vector<double> dydt1(N, 0.0);
  std::vector<double> k1(N), k2(N), k3(N), k4(N), aux(N), h1(N), h2(N), h3(N), h4(N), auxi(N);

  double time = 0;
  int nsteps = (tend - tinit)/dt;
  for(int ii = 0; ii < nsteps; ++ii) {
    time = 0.0 + ii*dt;
    // k1
    fderiv(y, dydt, time, l, q, F,omegaf );
    fderiv(w, dydt1, time, l, q1, F,omegaf );
    for(int ii = 0; ii < N; ++ii) {
      k1[ii] = dt*dydt[ii];
      h1[ii] = dt*dydt1[ii];
    }
    // k2 aux
    for(int ii = 0; ii < N; ++ii) {
      aux[ii] = y[ii] + k1[ii]/2;
      auxi[ii] = w[ii] + h1[ii]/2;
    }
    //k2
    fderiv(aux, dydt, time + dt/2, l, q, F, omegaf);
    fderiv(auxi, dydt1, time + dt/2, l, q1, F, omegaf);
    for(int ii = 0; ii < N; ++ii) {
      k2[ii] = dt*dydt[ii];
      h2[ii] = dt*dydt1[ii];
    }
    // k3 aux
    for(int ii = 0; ii < N; ++ii) {
      aux[ii] = y[ii] + k2[ii]/2;
      auxi[ii] = w[ii] + h2[ii]/2;
    }
    //k3
    fderiv(aux, dydt, time + dt/2, l, q, F,omegaf);
    fderiv(auxi, dydt1, time + dt/2, l, q1, F, omegaf);
    for(int ii = 0; ii < N; ++ii) {
      k3[ii] = dt*dydt[ii];
      h3[ii] = dt*dydt1[ii];
    }
    // k4 aux
    for(int ii = 0; ii < N; ++ii) {
      aux[ii] = y[ii] + k3[ii];
      auxi[ii] = w[ii] + h3[ii];
    }
    //k4
    fderiv(aux, dydt, time + dt, l, q, F, omegaf);
    fderiv(auxi, dydt1, time + dt/2, l, q1, F, omegaf);
    for(int ii = 0; ii < N; ++ii) {
      k4[ii] = dt*dydt[ii];
      h4[ii] = dt*dydt1[ii];
    }

    // write new data
    for(int ii = 0; ii < N; ++ii) {
      y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
      w[ii] = w[ii] + (h1[ii] + 2*h2[ii] + 2*h3[ii] + h4[ii])/6.0;
    }

     fout << time << "\t" <<std::abs( w[0]-y[0]) << std::endl;

  }
  fout.close();
}
