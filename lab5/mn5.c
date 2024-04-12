#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double K = 0.25;     // Stała reakcji
double initialGuess = 1.0; // Początkowa wartość D

// Funkcja f(x) i jej pochodna df(x)
double f(double x, double c0)
{
  return (x / ((c0 - x) * (c0 - x))) - K;
}

double df(double x, double c0, double h)
{
  return (f(x + h, c0) - f(x, c0)) / h;
}

// Metoda Newtona
double newtonMethod(double c0, double tolerance, double h)
{
  double x = initialGuess;
  double step;

  do
  {
    step = f(x, c0) / df(x, c0, h);
    x -= step;
  } while (fabs(step) > tolerance);

  return x;
}

int main()
{
  FILE *fp = fopen("dane.txt", "w");
  if (fp == NULL)
  {
    printf("Nie można otworzyć pliku do zapisu.\n");
    return 1;
  }

  double h = 1e-7;         // Krok do obliczania pochodnej
  double tolerance = 1e-6; // Tolerancja dla kroku metody
  double D, optimalC0 = 0;
  double minDifference = 1e6; // Duża wartość początkowa

  // Pętla dla różnych wartości c0
  for (double c0 = 2.0; c0 <= 10.0; c0 += 0.1)
  {
    D = newtonMethod(c0, tolerance, h);
    fprintf(fp, "%f %f %f\n", c0, D);

    // Zadanie 3 - poszukiwanie c0 dla którego D = B
    if (fabs(2 * D - c0) < minDifference)
    {
      minDifference = fabs(2 * D - c0);
      optimalC0 = c0;
    }
  }

  fclose(fp);
  printf("Wyniki zostały zapisane do pliku dane.txt\n");
  if (minDifference < tolerance)
  {
    printf("Szukane c0, dla którego D=B, wynosi około: %f\n", optimalC0);
  }
  else
  {
    printf("Nie znaleziono c0 spełniającego warunek D=B w danym zakresie.\n");
  }

  return 0;
