#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>

constexpr double angle = 50.8057242086663;

double optimize(double (*f)(double), double xa, double xb, double lambda1, double lambda2, double epsilon, std::ofstream &file)
{
  if (file.is_open())
  {
    file << "iter\tlog\n";
  }

  size_t iter = 1;
  while (std::abs(xb - xa) >= epsilon)
  {
    double x1 = xa + lambda1 * (xb - xa);
    double x2 = xa + lambda2 * (xb - xa);

    if (f(x1) > f(x2))
    {
      xa = x1;
    }
    else
    {
      xb = x1;
    }

    if (file.is_open())
    {
      file << iter << "\t" << std::log(std::abs((xa + xb) / 2 - angle)) << "\n";
    }

    ++iter;
  }

  return (xa + xb) / 2;
}

double ro(double h, double T, double P0)
{
  return 0.0289644 / (8.31446 * (T - 0.0065 * h)) * P0 * std::pow((1 - 0.0065 * h / T), 5.255);
}

double F_drag(double v, double Cx, double Sm, double h, double T, double P0)
{
  return Cx * ro(h, T, P0) * Sm * std::pow(v, 2) / 2;
}

double wziuuum(double alfa, std::ofstream &file)
{
  double x, y, v, Time;
  double pi, g, Cx, diam, Sm, M, T, P0;
  double dt, v_tmp;

  pi = 3.14159265358979;
  g = 9.807;

  Cx = 0.187;
  diam = 155;
  Sm = pi * std::pow((diam / 2 / 1000), 2);

  M = 48;
  v = 958.0;

  T = 273;
  P0 = 101300;

  alfa = alfa / 180 * pi;

  x = 0.0;
  y = 1.0E-5;
  Time = 0.0;

  if (file.is_open())
  {
    file << "x\ty\n";
    file << x << "\t" << y << "\n";
  }

  while (y > 0.0)
  {
    dt = 1.0 / v;

    x += v * std::cos(alfa) * dt;
    y += v * std::sin(alfa) * dt;

    v_tmp = v + (-F_drag(v, Cx, Sm, y, T, P0) / M - g * std::sin(alfa)) * dt;
    alfa -= g * std::cos(alfa) / v * dt;
    v = v_tmp;
    Time += dt;

    if (file.is_open())
    {
      file << x << "\t" << y << "\n";
    }
  }

  return x;
}

double negative_wziuuum(double alpha)
{
  std::ofstream empty_file;
  return -wziuuum(alpha, empty_file);
}

double wziuuum_i_bum(double alfa, double distance, double elevation, std::ofstream &file)
{

  double x, y, v, Time;
  double pi, g, Cx, diam, Sm, M, T, P0;
  double dt, v_tmp;

  pi = 3.14159265358979;
  g = 9.807;

  Cx = 0.187;

  diam = 155;
  Sm = pi * std::pow((diam / 2 / 1000), 2);

  M = 48;
  v = 958.0;

  T = 273;
  P0 = 101300;

  alfa = alfa / 180 * pi;

  x = 0.0;
  y = 1.0E-5;
  Time = 0.0;

  if (file.is_open())
  {
    file << "x\ty\n";
    file << x << "\t" << y << "\n";
  }

  while (y > 0.0 && x < distance)
  {
    dt = 1.0 / v;

    x += v * std::cos(alfa) * dt;
    y += v * std::sin(alfa) * dt;

    v_tmp = v + (-F_drag(v, Cx, Sm, y, T, P0) / M - g * std::sin(alfa)) * dt;
    alfa -= g * std::cos(alfa) / v * dt;
    v = v_tmp;
    Time += dt;

    if (file.is_open())
    {
      file << x << "\t" << y << "\n";
    }
  }

  return std::sqrt(std::pow(x - distance, 2) + std::pow(y - elevation, 2));
}

auto bum0 = [](double alpha)
{
  std::ofstream empty_file;
  return wziuuum_i_bum(alpha, 30000, 0, empty_file);
};

auto bum300 = [](double alpha)
{
  std::ofstream empty_file;
  return wziuuum_i_bum(alpha, 30000, 300, empty_file);
};

int main()
{
  double epsilon = 1e-6;
  double r = (std::sqrt(5) - 1) / 2;

  std::ofstream empty_file;
  double optimal_angle = optimize(negative_wziuuum, 20, 70, std::pow(r, 2), r, epsilon, empty_file);
  double optimal_angle_range = wziuuum(optimal_angle, empty_file);
  std::cout << "Maksymalny zasieg: " << optimal_angle << ", wynosi: " << optimal_angle_range << "\n";

  std::ofstream optimal_angle_file("optimal_angle.csv");
  optimal_angle_file << "angle\trange\n";
  optimal_angle_file << optimal_angle << "\t" << optimal_angle_range << "\n";
  optimal_angle_file.close();

  std::ofstream exer3a("3a.csv");
  exer3a << "angle\trange\n";
  for (double angle = 1; angle <= 90; angle += 0.5)
  {
    exer3a << angle << "\t" << wziuuum(angle, empty_file) << "\n";
  }
  exer3a.close();

  std::ofstream exer3b("3b.csv");
  wziuuum(optimal_angle, exer3b);
  exer3b.close();

  std::ofstream exer3c_golden("3c_golden.csv");
  optimize(negative_wziuuum, 20, 70, std::pow(r, 2), r, epsilon, exer3c_golden);
  exer3c_golden.close();

  std::ofstream exer3c_b("3c_b.csv");
  optimize(negative_wziuuum, 20, 70, 1.0 / 3.0, 2.0 / 3.0, epsilon, exer3c_b);
  exer3c_b.close();

  std::ofstream exer4a("4a.csv");
  double ans0 = optimize(bum0, 20, 70, std::pow(r, 2), r, epsilon, exer4a);
  wziuuum_i_bum(ans0, 30000, 0, exer4a);
  exer4a.close();

  std::ofstream exer4b("4b.csv");
  double ans300 = optimize(bum300, 20, 70, std::pow(r, 2), r, epsilon, exer4b);
  wziuuum_i_bum(ans300, 30000, 300, exer4b);
  exer4b.close();

  return 0;
}
