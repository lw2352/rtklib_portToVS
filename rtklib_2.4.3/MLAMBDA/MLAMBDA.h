
extern "C" int test_lambda(int n, int m, double* a, double* Q, double* F, double* s);
extern "C" int test_filter(double* x, double* P, double* H, double* v, double* R, int n, int m, double* xp, double* Pp);
extern "C" int resultFilter(double* x, double* y, double* z, int n);
extern "C" void test();
extern "C" void test_udPos(int n, double* F_in, double* x_in, double* P_in, double* xp_out, double* Pp_out);

extern "C" int resultSTD(double* rr, int n, double* out);
extern "C" int test_lsq(double* A_in, double* y_in, int n, int m, double* x_in, double* Q_in);
extern "C" void testVel(double* Ir_in, double* Ib_in, double* Vs_in, double* Fr_in, double* Fb_in, double* x_in, int n, double lam);
