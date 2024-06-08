#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

// Funkcje do całkowania
double f1(double x) { return sin(x); }
double f2(double x) { return x * sin(5 * x); }

// Metoda trapezów
double trapezoidal(double a, double b, int n, double (*f)(double))
{
  double h = (b - a) / (n - 1);
  double sum = 0.0;
  for (int i = 0; i < n; i++)
  {
    double x = a + i * h;
    if (i == 0 || i == n - 1)
    {
      sum += f(x) / 2;
    }
    else
    {
      sum += f(x);
    }
  }
  return sum * h;
}

// Metoda Simpsona
double simpson(double a, double b, int n, double (*f)(double))
{

  double h = (b - a) / (n - 1);
  double sum = 0.0;
  for (int i = 0; i < n; i++)
  {
    double x = a + i * h;
    if (i == 0 || i == n - 1)
    {
      sum += f(x);
    }
    else if (i % 2 == 0)
    {
      sum += 2 * f(x);
    }
    else
    {
      sum += 4 * f(x);
    }
  }
  return sum * h / 3.0;
}

int main()
{
  const double PI = acos(-1);
  const double I1_exact = 2.0;
  const double I2_exact = PI / 5;
  vector<int> numNodes = {11, 21, 41, 61, 81, 101, 121, 141, 161, 181, 201};

  ofstream file1("wyniki_I1.csv");
  ofstream file2("wyniki_I2.csv");

  file1 << "Liczba wezlow,Blad trapezow,Blad Simpsona" << endl;
  file2 << "Liczba wezlow,Blad trapezow,Blad Simpsona" << endl;

  for (int n : numNodes)
  {
    double I1_trapez = trapezoidal(0, PI, n, f1);
    double I1_simpson = simpson(0, PI, n, f1);
    double I2_trapez = trapezoidal(0, PI, n, f2);
    double I2_simpson = simpson(0, PI, n, f2);

    double error_I1_trapez = abs(I1_exact - I1_trapez);
    double error_I1_simpson = abs(I1_exact - I1_simpson);
    double error_I2_trapez = abs(I2_exact - I2_trapez);
    double error_I2_simpson = abs(I2_exact - I2_simpson);

    // Logarytm z błędu (z zabezpieczeniem przed log(0))
    double log_error_I1_trapez = (error_I1_trapez > 0) ? log10(error_I1_trapez) : -10; // -10 jako wartość zastępcza dla log(0)
    double log_error_I1_simpson = (error_I1_simpson > 0) ? log10(error_I1_simpson) : -10;
    double log_error_I2_trapez = (error_I2_trapez > 0) ? log10(error_I2_trapez) : -10;
    double log_error_I2_simpson = (error_I2_simpson > 0) ? log10(error_I2_simpson) : -10;

    file1 << n << "," << log_error_I1_trapez << "," << log_error_I1_simpson << endl;
    file2 << n << "," << log_error_I2_trapez << "," << log_error_I2_simpson << endl;
  }

  file1.close();
  file2.close();

  return 0;
}