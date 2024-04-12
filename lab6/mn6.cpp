
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// Funkcje do interpolacji
double f1(double x) { return exp(-x * x); }
double f2(double x) { return cos(2 * x); }

// Wielomian interpolacyjny Lagrange'a
double lagrangeInterpolation(const std::vector<double> &x, const std::vector<double> &y, double x2)
{
  double result = 0.0;
  for (int i = 0; i < x.size(); i++)
  {
    double term = y[i];
    for (int j = 0; j < x.size(); j++)
    {
      if (j != i)
        term = term * (x2 - x[j]) / double(x[i] - x[j]);
    }
    result += term;
  }
  return result;
}

// Generowanie węzłów interpolacji z wykorzystaniem wielomianu Czebyszewa
std::vector<double> generateChebyshevNodes(double xmin, double xmax, int n)
{
  std::vector<double> nodes(n);
  for (int i = 0; i < n; i++)
  {
    nodes[i] = 0.5 * ((xmax - xmin) * cos(M_PI * (2.0 * i + 1) / (2.0 * n + 2)) + (xmax + xmin));
  }
  return nodes;
}

// Generowanie węzłów interpolacji
std::vector<double> generateNodes(double xmin, double xmax, int n)
{
  std::vector<double> nodes(n);
  double step = (xmax - xmin) / (n - 1);
  for (int i = 0; i < n; i++)
    nodes[i] = xmin + i * step;
  return nodes;
}

// Obliczanie wartości funkcji w węzłach
std::vector<double> calculateValues(const std::vector<double> &nodes, double (*func)(double))
{
  std::vector<double> values(nodes.size());
  for (int i = 0; i < nodes.size(); i++)
    values[i] = func(nodes[i]);
  return values;
}

// Wykonanie interpolacji dla różnych wartości n
void performInterpolation(double (*func)(double), const std::string &filename)
{
  std::ofstream file(filename);
  file << "x,y_actual,y_interpolated\n";

  for (int n : {5, 10, 15, 20})
  {
    std::vector<double> nodes = generateChebyshevNodes(-5, 5, n);
    std::vector<double> values = calculateValues(nodes, func);
    for (double x = -5; x <= 5; x += 0.1)
    {
      double y_actual = func(x);
      double y_interpolated = lagrangeInterpolation(nodes, values, x);
      file << x << "," << y_actual << "," << y_interpolated << "\n";
    }
    file << "NEXT N\n";
  }

  file.close();
}

int main()
{
  performInterpolation(f1, "interpolationCzebyszew_f1.csv");
  performInterpolation(f2, "interpolationCzebyszew_f2.csv");
  return 0;
}
