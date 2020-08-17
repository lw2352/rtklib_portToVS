
extern "C" int test_lambda(int n, int m, double* a, double* Q, double* F, double* s);
extern "C" int test_filter(double* x, double* P, double* H, double* v, double* R, int n, int m, double* xp, double* Pp);
extern "C" int resultFilter(double* x, double* y, double* z, int n);
extern "C" void test();
