#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"
#include "ProyectoF.h"



TEST_CASE( "Verificando Simulaciones", "[pendulo]" ) {

  std::vector<double> F(4);
  std::ifstream fin("input.txt");
  fin >> F[0];
  fin >> F[1];
  fin >> F[2];
  fin >> F[3];

  std::cout << "Verificando que la fuerza de la simulación del primer punto sea nula" << "\n";
  REQUIRE (F[0] == 0);
  std::cout << "Verificando que las fuerzas de las simulaciones del segundo punto sean diferentes" << "\n";
  REQUIRE (F[1] != F[2]);
  REQUIRE (F[2] != F[3]);
  REQUIRE (F[1] != F[3]);
  fin.close();

  double tinit,thetainit=0.2;
  const int N = 2;

    std::vector<double> y(N);
    std::vector<double> m(150);
    std::vector<double> n(150);
    initial_conditions(y,0.2, 0.0 );

  
    m = Simulacion_test(y, 0.0, 60, 0.4, 9.8, 0.5, 0.1, 0.66);
    std::cout << "Verificando para amplitud inicial del pendulo: 0.2 y una fuerza menor al amortiguamiento" <<"\n";
    for(int ii = 0; ii < 150; ++ii) {
     // std::cout << m[ii] << "\n";

        double a = std::abs(m[ii]);
        REQUIRE(a <= thetainit);} //se verifica que la amplitud con una fuerza muy pequeña sea menor a la amplitud inicial

    std::cout << "Verificando para amplitud inicial del péndulo: 1.0 y una fuerza menor al amortiguamiento" << "\n";
    initial_conditions(y, 1.0, 0.0);
    thetainit = 1.0;
    n = Simulacion_test(y, 0.0, 60, 0.4, 9.8, 0.5, 0.1, 0.66);
     for(int ii = 0; ii < 150; ++ii) {

        double a = std::abs(n[ii]);
        REQUIRE(a <=thetainit );}



}

