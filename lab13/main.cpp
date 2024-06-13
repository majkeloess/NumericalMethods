#include <iostream>
#include <fstream>
#include <vector>

void volterraLotka(double t, const std::vector<double> &y, std::vector<double> &dydt)
{
  const double a = 1.0;
  const double b = 1.0;
  const double c = 1.0;
  const double d = 1.0;

  dydt[0] = a * y[0] - b * y[0] * y[1]; // dK/dt
  dydt[1] = c * y[0] * y[1] - d * y[1]; // dZ/dt
}

// Metoda Eulera
void eulerMethod(double t0, double tK, double dt, const std::vector<double> &y0, std::ofstream &outfile)
{
  std::vector<double> y = y0;
  double t = t0;

  while (t <= tK)
  {
    outfile << t << "," << y[0] << "," << y[1] << std::endl;

    std::vector<double> dydt(2);
    volterraLotka(t, y, dydt);

    for (int i = 0; i < 2; ++i)
    {
      y[i] += dt * dydt[i];
    }
    t += dt;
  }
}

void rungeKutta2(double t0, double tK, double dt, const std::vector<double> &y0, std::ofstream &outfile)
{
  std::vector<double> y = y0;
  double t = t0;

  while (t <= tK)
  {
    outfile << t << "," << y[0] << "," << y[1] << std::endl;

    std::vector<double> k1(2), k2(2);

    // k1
    volterraLotka(t, y, k1);

    // k2
    for (int i = 0; i < 2; ++i)
    {
      k2[i] = y[i] + dt * k1[i];
    }
    volterraLotka(t + dt, k2, k2);

    // Aktualizacja y
    for (int i = 0; i < 2; ++i)
    {
      y[i] += (dt / 2.0) * (k1[i] + k2[i]);
    }

    t += dt;
  }
}

int main()
{
  double t0 = 0.0;
  double tK = 40.0;
  std::vector<double> dtValues = {0.1, 0.01, 0.001};
  std::vector<double> y0 = {1.0, 0.5}; // Warunki początkowe

  for (double dt : dtValues)
  {
    std::ofstream outfileEuler("euler_" + std::to_string(dt) + ".csv");
    outfileEuler << "t,K,Z" << std::endl;
    eulerMethod(t0, tK, dt, y0, outfileEuler);
    outfileEuler.close();

    std::ofstream outfileRK4("rk4_" + std::to_string(dt) + ".csv");
    outfileRK4 << "t,K,Z" << std::endl;
    rungeKutta2(t0, tK, dt, y0, outfileRK4);
    outfileRK4.close();
  }

  return 0;
}
