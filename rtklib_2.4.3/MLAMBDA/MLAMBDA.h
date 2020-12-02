
extern "C" int test_lambda(int n, int m, double* a, double* Q, double* F, double* s);
extern "C" int test_filter(double* x, double* P, double* H, double* v, double* R, int n, int m, double* xp, double* Pp);
extern "C" int resultFilter(double* x, double* y, double* z, int n);
extern "C" void test();
extern "C" void test_udPos(int n, double* F_in, double* x_in, double* P_in, double* xp_out, double* Pp_out);

extern "C" int resultSTD(double* rr, int n, double* out);
