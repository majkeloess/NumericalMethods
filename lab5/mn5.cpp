#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <fstream>

// Definicja funkcji
double function(double D, void *params)
{
  double c0 = *(double *)params;
  double K = 0.25;

  return (D / ((c0 - D) * (c0 - D))) - K;
}

// Definicja pochodnej funkcji
double derivative(double D, void *params)
{
  double h = 1e-7;
  double c0 = *(double *)params;
  double K = 0.25;

  return ((function(D + h, params) - function(D, params)) / h);
}

// Implementacja metody Newtona
double newtonMethod(double c0)
{
  int iter = 0, max_iter = 100;
  double D = 0.5, D0, abs, tol = 1e-6;

  do
  {
    iter++;
    D0 = D;
    D = D0 - function(D0, &c0) / derivative(D0, &c0);
    abs = fabs(D - D0);
  } while (abs > tol && iter < max_iter);

  return D;
}

int main()
{
  std::ofstream file;
  file.open("point2.csv");
  file << "\t"
       << "c0"
       << "\t"
       << "D"
       << "\t"
       << "B" << std::endl;
  for (double c0 = 2.0; c0 <= 10.0; c0 += 0.01)
  {
    double D = newtonMethod(c0);
    double B = c0 - D;
    file << "\t" << c0 << "\t" << D << "\t" << B << std::endl;
  }

  file.close();

  return 0;
}
