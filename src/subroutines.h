
void SWP( double **X, int k, int size);
void dinv(double **X, int size, double **X_inv);
void dcholdc(double **X, int size, double **L);
double ddet(double **X, int size, int give_log);
double mean(double *X, int n);
double var(double *X, int n, int unbiased);
double max(double *X, int n);
