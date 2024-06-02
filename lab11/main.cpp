
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

using namespace std;

void generate_signal(double *signal, int N, double omega)
{
  for (int i = 0; i < N; i++)
  {
    signal[2 * i] = sin(omega * i) + sin(2 * omega * i) + sin(3 * omega * i);
    signal[2 * i + 1] = 0.0;
  }
}

void add_noise(double *signal, int N)
{
  for (int i = 0; i < N; i++)
  {
    double delta = 2.0 * (rand() / (double)RAND_MAX) - 1.0;
    signal[2 * i] += delta;
  }
}

void save_to_csv(const string &filename, double *signal, int N)
{
  ofstream file(filename);
  if (!file.is_open())
  {
    cerr << "Failed to open file: " << filename << endl;
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < N; i++)
  {
    file << i << ", " << signal[2 * i] << ", " << signal[2 * i + 1] << endl;
  }
  file.close();
}

void save_magnitude_to_csv(const char *filename, double *signal, int N)
{
  FILE *file = fopen(filename, "w");
  if (!file)
  {
    perror("Failed to open file");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < N; i++)
  {
    double real = signal[2 * i];
    double imag = signal[2 * i + 1];
    double magnitude = sqrt(real * real + imag * imag);
    fprintf(file, "%d, %f\n", i, magnitude);
  }
  fclose(file);
}

void plot_with_gnuplot(const char *filename, const char *title, const char *outputfile, const char *xlabel, const char *ylabel, const char *linecolor = "blue")
{
  char command[512];
  sprintf(command, "gnuplot -e \"set terminal png; set output '%s'; set title '%s'; set xlabel '%s'; set ylabel '%s'; "
                   "set key inside right top; plot '%s' using 1:2 with lines linecolor '%s' title '%s'\"",
          outputfile, title, xlabel, ylabel, filename, linecolor, title); // Dodajemy linecolor
  system(command);
}

void plot_signals_with_gnuplot(const char *noisy_filename, const char *denoised_filename, const char *clean_filename, const char *title, const char *outputfile)
{
  char command[512];
  sprintf(command, "gnuplot -e \"set terminal png; set output '%s'; set title '%s'; set xlabel 'Index'; set ylabel 'Value'; "
                   "set key inside right top; set pointsize 0.5; "
                   "plot '%s' using 1:2 with points pointtype 7 linecolor 'red' title 'zaburzony', " // Czerwony dla zaszumionego
                   "'%s' using 1:2 with lines linecolor 'green' title 'odszumiony', "                // Zielony dla odszumionego
                   "'%s' using 1:2 with points pointtype 7 linecolor 'blue' title 'niezaburzony'\"", // Niebieski dla oryginalnego
          outputfile, title, noisy_filename, denoised_filename, clean_filename);
  system(command);
}

void plot_fft_components_with_gnuplot(const char *fft_filename, const char *magnitude_filename, const char *outputfile_real_imag, const char *outputfile_magnitude)
{
  char command[512];
  sprintf(command, "gnuplot -e \"set terminal png; set output '%s'; set title 'FFT'; set xlabel 'Index'; set ylabel 'Value'; "
                   "set key inside right top; plot '%s' using 1:2 with lines linecolor 'orange' title 'Real part', '%s' using 1:3 with lines linecolor 'purple' title 'Imaginary part'\"",
          outputfile_real_imag, fft_filename, fft_filename); // Pomarańczowy dla części rzeczywistej, fioletowy dla urojonej
  system(command);

  sprintf(command, "gnuplot -e \"set terminal png; set output '%s'; set title 'Magnitude of FFT'; set xlabel 'Index'; set ylabel 'Magnitude'; "
                   "set key inside right top; plot '%s' using 1:2 with lines linecolor 'brown' title 'Magnitude'\"",
          outputfile_magnitude, magnitude_filename); // Brązowy dla modułów
  system(command);
}

int main()
{
  int ks[] = {8, 10, 12};
  for (int idx = 0; idx < 3; idx++)
  {
    int k = ks[idx];
    int N = 1 << k;
    double omega = 2 * M_PI / N;
    double *signal = (double *)malloc(2 * N * sizeof(double));
    double *clean_signal = (double *)malloc(2 * N * sizeof(double));
    if (!signal || !clean_signal)
    {
      perror("Failed to allocate memory");
      exit(EXIT_FAILURE);
    }

    generate_signal(clean_signal, N, omega);
    char clean_filename[128];
    sprintf(clean_filename, "clean_k%d.csv", k);
    save_to_csv(clean_filename, clean_signal, N);

    generate_signal(signal, N, omega);
    add_noise(signal, N);
    char noisy_filename[128];
    sprintf(noisy_filename, "noisy_k%d.csv", k);
    save_to_csv(noisy_filename, signal, N);

    gsl_fft_complex_radix2_forward(signal, 1, N);

    char fft_filename[128];
    sprintf(fft_filename, "fft_k%d.csv", k);
    save_to_csv(fft_filename, signal, N);

    char magnitude_filename[128];
    sprintf(magnitude_filename, "magnitude_fft_k%d.csv", k);
    save_magnitude_to_csv(magnitude_filename, signal, N);

    char plot_filename[128];
    sprintf(plot_filename, "magnitude_fft_k%d.png", k);
    plot_with_gnuplot(magnitude_filename, "Magnitude FFT", plot_filename, "Index", "Magnitude");

    double max_c = 0.0;
    for (int i = 0; i < N; i++)
    {
      double real = signal[2 * i];
      double imag = signal[2 * i + 1];
      double magnitude = sqrt(real * real + imag * imag);
      if (magnitude > max_c)
      {
        max_c = magnitude;
      }
    }
    double threshold = max_c / 2.0;
    for (int i = 0; i < N; i++)
    {
      double real = signal[2 * i];
      double imag = signal[2 * i + 1];
      double magnitude = sqrt(real * real + imag * imag);
      if (magnitude < threshold)
      {
        signal[2 * i] = 0.0;
        signal[2 * i + 1] = 0.0;
      }
    }

    gsl_fft_complex_radix2_backward(signal, 1, N);

    for (int i = 0; i < N; i++)
    {
      signal[2 * i] /= N;
      signal[2 * i + 1] /= N;
    }

    char denoised_filename[128];
    sprintf(denoised_filename, "denoised_signal_k%d.csv", k);
    save_to_csv(denoised_filename, signal, N);

    char signal_plot_filename[128];
    sprintf(signal_plot_filename, "signal_plot_k%d.png", k);
    plot_signals_with_gnuplot(noisy_filename, denoised_filename, clean_filename, "Porównanie sygnałów", signal_plot_filename);

    if (k == 8)
    {
      char fft_real_imag_plot_filename[128];
      sprintf(fft_real_imag_plot_filename, "fft_real_imag_k%d.png", k);
      char magnitude_fft_plot_filename[128];
      sprintf(magnitude_fft_plot_filename, "magnitude_fft_k%d.png", k);
      plot_fft_components_with_gnuplot(fft_filename, magnitude_filename, fft_real_imag_plot_filename, magnitude_fft_plot_filename);
    }

    free(signal);
    free(clean_signal);
  }

  return 0;
}