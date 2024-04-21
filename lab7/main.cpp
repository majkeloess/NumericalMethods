#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <cmath>
#include <string>

double func_a(double x) { return exp(-x * x); }
double func_b(double x) { return x < 0 ? -1.0 : 1.0; }
double func_c(double x) { return cos(2 * x); }

void interpolation(std::function<double(double)> f, int n, const std::string &func_name)
{
    std::string filename = func_name + ".csv";
    std::ofstream file(filename);

    double xi, yi, x_interp, y_interp, y_deriv2;
    double x_min = -5.0, x_max = 5.0;
    double h = (x_max - x_min) / (n - 1);

    std::vector<double> x(n), y(n), y2(n);
    for (int i = 0; i < n; i++)
    {
        xi = x_min + i * h;
        x[i] = xi;
        y[i] = f(xi);
    }

    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(spline, x.data(), y.data(), n);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    for (int i = 0; i < n; i++)
    {
        y2[i] = gsl_spline_eval_deriv2(spline, x[i], acc);
    }

    file << "x_interp\t y_inter\t y_deriv2\t\n";
    for (x_interp = x_min; x_interp <= x_max; x_interp += 0.1)
    {
        y_interp = gsl_spline_eval(spline, x_interp, acc);
        y_deriv2 = gsl_spline_eval_deriv2(spline, x_interp, acc);
        file << x_interp << "\t" << y_interp << "\t" << y_deriv2 << "\t\n";
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

int main()
{
    interpolation(func_a, 5, "exp5");
    interpolation(func_a, 10, "exp10");
    interpolation(func_a, 20, "exp20");

    interpolation(func_b, 5, "lin5");
    interpolation(func_b, 10, "lin10");
    interpolation(func_b, 20, "lin20");

    interpolation(func_c, 5, "cos5");
    interpolation(func_c, 10, "cos10");
    interpolation(func_c, 20, "cos20");

    return 0;
}
