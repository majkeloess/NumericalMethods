int main() {
    int A[10][10];
    int sum = 0;

    for (int i = 0; i < 10; i++) {
        for (int j = i+1; j < 10; j++) {
            sum += A[i][j];
        }
  }
  
  printf("%d", sum);

  return 0;
}