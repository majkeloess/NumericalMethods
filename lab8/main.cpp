#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

double ro(double h, double T, double P0)
{
  return 0.0289644 / (8.314446 * (T - 0.0065 * h)) * P0 * pow((1 - 0.0065 * h / T), 5.255);
}

void saveDataToFile(const std::string &filename, const std::vector<double> &x, const std::vector<double> &y)
{
  std::ofstream file;
  file.open(filename);
  for (int i = 0; i < x.size(); i++)
  {
    file << x[i] << " " << y[i] << "\n";
  }
  file.close();
}

double F_drag(double v, double Cx, double Sm, double h, double T, double P0)
{
  return Cx * ro(h, T, P0) * Sm * pow(v, 2) / 2;
}

double wziuuum(double alfa)
{
  double x, y, v, Time;
  double pi, g, Cx, diam, Sm, M, T, P0;
  double dt, v_tmp;

  pi = 3.14159265358979;
  g = 9.807;

  Cx = 0.187;

  diam = 155;
  Sm = pi * pow((diam / 2 / 1000), 2);

  M = 48;
  v = 958.0;

  T = 273;
  P0 = 101300;

  alfa = alfa / 180 * pi;

  dt = 0.01;

  x = 0.0;
  y = 1.0E-5;
  Time = 0.0;

  while (y > 0.0)
  {
    x += v * cos(alfa) * dt;
    y += v * sin(alfa) * dt;

    v_tmp = v + (-F_drag(v, Cx, Sm, y, T, P0) / M - g * sin(alfa)) * dt;
    alfa -= g * cos(alfa) / v * dt;
    v = v_tmp;
    Time += dt;
  }

  return x;
}

double goldenSectionSearch(double a, double b, double epsilon)
{
  double phi = (sqrt(5.0) - 1.0) / 2.0;
  double x1, x2;

  while (fabs(b - a) > epsilon)
  {
    x1 = b - (b - a) * phi;
    x2 = a + (b - a) * phi;

    if (-wziuuum(x1) >= -wziuuum(x2))
      a = x1;
    else
      b = x2;
  }

  return (a + b) / 2.0;
}

int main()
{
  double a = 20.0;       // dolna granica
  double b = 70.0;       // górna granica
  double epsilon = 1e-6; // dokładność

  double max_angle = goldenSectionSearch(a, b, epsilon);

  std::cout << "Kąt ostrzału, który daje maksymalny zasięg wynosi: " << max_angle << " stopni\n";

  // Przygotowanie danych do zapisania do pliku
  std::vector<double> angles, ranges;
  for (double angle = a; angle <= b; angle += 0.1)
  {
    angles.push_back(angle);
    ranges.push_back(wziuuum(angle));
  }

  // Zapisanie danych do pliku
  saveDataToFile("ranges.csv", angles, ranges);

  return 0;
}
