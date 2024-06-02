#include <iostream>
#include <fstream>
#include <cmath>

long int U1()
{
  static long int X = 10;
  const long int a = 17;
  const long int m = pow(2, 13) - 1;
  X = (a * X) % m;
  return X;
}

long int U2()
{
  static long int X = 10;
  const long int a = 85;
  const long int m = pow(2, 13) - 1;
  X = (a * X) % m;
  return X;
}

long int U3()
{
  static long int X[3] = {10, 10, 10};
  const long int a = 1176;
  const long int b = 1476;
  const long int c = 1176;
  const long int m = pow(2, 32) - 5;
  X[2] = (a * X[2] + b * X[1] + c * X[0]) % m;
  X[0] = X[1];
  X[1] = X[2];
  return X[2];
}

double U1_normalized()
{
  static long int X = 10;
  const long int a = 17;
  const long int m = pow(2, 13) - 1;
  X = (a * X) % m;
  return static_cast<double>(X) / (m - 1);
}

double U2_normalized()
{
  static long int X = 10;
  const long int a = 85;
  const long int m = pow(2, 13) - 1;
  X = (a * X) % m;
  return static_cast<double>(X) / (m - 1);
}

double U3_normalized()
{
  static long int X[3] = {10, 10, 10};
  const long int a = 1176;
  const long int b = 1476;
  const long int c = 1176;
  const long int m = pow(2, 32) - 5;
  X[2] = (a * X[2] + b * X[1] + c * X[0]) % m;
  X[0] = X[1];
  X[1] = X[2];
  return static_cast<double>(X[2]) / (m - 1);
}

int main()
{
  std::ofstream fileU1("U1.csv");
  std::ofstream fileU2("U2.csv");
  std::ofstream fileU3("U3.csv");

  // Inicjalizacja pierwszych wartości (dla wykresów x_i = f(x_{i-2}) i x_i = f(x_{i-3}))
  U1();
  U1();
  U2();
  U2();
  U3();
  U3();

  for (int i = 1; i <= 1000; i++)
  {
    long int x1 = U1();
    long int x2 = U2();
    long int x3 = U3();

    fileU1 << i << "\t" << x1 << std::endl;
    fileU2 << i << "\t" << x2 << std::endl;
    fileU3 << i << "\t" << x3 << std::endl;

    // Dodatkowe kolumny dla wykresów x_i = f(x_{i-2}) i x_i = f(x_{i-3})
    if (i > 1)
    {
      fileU1 << i << "\t" << U1() << std::endl; // x_i = f(x_{i-2})
      fileU2 << i << "\t" << U2() << std::endl;
      fileU3 << i << "\t" << U3() << std::endl;
    }
    if (i > 2)
    {
      fileU1 << i << "\t" << U1() << std::endl; // x_i = f(x_{i-3})
      fileU2 << i << "\t" << U2() << std::endl;
      fileU3 << i << "\t" << U3() << std::endl;
    }
  }

  fileU1.close();
  fileU2.close();
  fileU3.close();

  // const int n = 20000;
  // const double l = 1.0;
  // const double r = 1.0;
  // int insideCircle = 0;
  // double piApproximation;

  // filePi << "Iteracja,Pi ,Błąd" << std::endl;

  // for (int i = 1; i <= n; i++)
  // {
  //   double x = l * U3_normalized();
  //   double y = l * U3_normalized();

  //   if (x * x + y * y <= r * r)
  //   {
  //     insideCircle++;
  //   }

  //   if (i % 100 == 0)
  //   {
  //     piApproximation = 4.0 * static_cast<double>(insideCircle) / i;
  //     double error = std::abs(M_PI - piApproximation);
  //     filePi << i << "," << piApproximation << "," << error << std::endl;
  //   }
  // }

  fileU1.close();
  fileU2.close();
  fileU3.close();
  // filePi.close();

  return 0;
}
