// MLAMBDA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "ARLambda.h"
#include "MLAMBDA.h"

using namespace std;
using namespace Eigen;

// Shows how to utilize MLAMBDA-Eigen
int main0()
{

	// Lets assume you have the float ambiguities and corresponding covariance matrix
	// I have a sample here, taken from a GPS-only Relative kinematic processing routine,
	// The ambiguites belong to double-differenced carrier phase observations
	// NOTE : The units are in cycles
	// NOTE : Covariance values are large because it is taken from initial stages of processing
	int n = 3;
	VectorXd floatAmb = VectorXd::Zero(n);
	MatrixXd floatAmbCov = MatrixXd::Zero(n, n);

	/*floatAmb << -9.75792, 22.1086, -1.98908, 3.36186, 23.2148, 7.75073;

	floatAmbCov << 0.0977961, 0.0161137, 0.0468261, 0.0320695, 0.080857, 0.0376408,
		0.0161137, 0.0208976, 0.0185378, 0.00290225, 0.0111409, 0.0247762,
		0.0468261, 0.0185378, 0.0435412, 0.0227732, 0.0383208, 0.0382978,
		0.0320695, 0.00290225, 0.0227732, 0.0161712, 0.0273471, 0.0154774,
		0.080857, 0.0111409, 0.0383208, 0.0273471, 0.0672121, 0.0294637,
		0.0376408, 0.0247762, 0.0382978, 0.0154774, 0.0294637, 0.0392536;*/

	floatAmb << 5.45, 3.1, 2.97;
	floatAmbCov <<  6.290, 5.978, 0.544,
					5.978, 6.292, 2.340,
					0.544, 2.340, 6.288;

	// Try Ambiguity Resolution...
	
	ARLambda AR;
	VectorXd intAmb(floatAmb.size()); intAmb.setZero();
	intAmb = AR.resolveIntegerAmbiguity(floatAmb, floatAmbCov);

	cout << "Float Ambiguity: \t" << floatAmb.transpose() << "\n";
	cout << "Integer Ambiguity: \t" << intAmb.transpose() << "\n";
	return 0;
}

/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
int test_lambda(int n, int m, double* a, double* Q, double* F, double* s)
{
	//m固定是2，F是n*2的矩阵，s是两个元素数组
	//VectorXd floatAmb = VectorXd::Zero(n);//a
	//MatrixXd floatAmbCov = MatrixXd::Zero(n, n);//Q
	MatrixXd F_in;
	VectorXd S_in;
	int ret;

	/*floatAmb << 5.45, 3.1, 2.97;
	floatAmbCov << 6.290, 5.978, 0.544,
		5.978, 6.292, 2.340,
		0.544, 2.340, 6.288;*/

	VectorXd floatAmb = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(a, n, 1);
	MatrixXd floatAmbCov = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(Q, n, n);

	ARLambda AR;
	VectorXd intAmb(floatAmb.size()); intAmb.setZero();
	intAmb = AR.resolveIA(floatAmb, floatAmbCov, F_in, S_in, &ret);
	//cout << "Integer Ambiguity: \t" << intAmb.transpose() << "\n";
	/*for(int i=0;i<n*2;i++)
	{
		F[i] = F_in[i];
	}*/
	Map<MatrixXd>(F, n, 2) = F_in;
	Map<MatrixXd>(s, 2, 1) = S_in;
	return ret;
}

/* kalman filter ---------------------------------------------------------------
* kalman filter state update as follows:
*
*   K=P*H*(H'*P*H+R)^-1, xp=x+K*v, Pp=(I-K*H')*P
*
* args   : double *x        I   states vector (n x 1)
*          double *P        I   covariance matrix of states (n x n)
*          double *H        I   transpose of design matrix (n x m)
*          double *v        I   innovation (measurement - model) (m x 1)
*          double *R        I   covariance matrix of measurement error (m x m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   states vector after update (n x 1)
*          double *Pp       O   covariance matrix of states after update (n x n)
* return : status (0:ok,<0:error)
* notes  : matirix stored by column-major order (fortran convention)
*          if state x[i]==0.0, not updates state x[i]/P[i+i*n]
*-----------------------------------------------------------------------------*/
int test_filter_(double* x_in, double* P_in, double* H_in,
 double* v_in, double* R_in, int n, int m,
	double* xp_out, double* Pp_out)
{
	if (x_in[0] == 0.0)
		return -1;
	VectorXd x= Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(x_in, n, 1);
	MatrixXd P = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(P_in, n, n);
	MatrixXd H = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(H_in, n, m);
	VectorXd v = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(v_in, m, 1);
	MatrixXd R = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(R_in, m, m);

	double x1[8] = {0}, p1[8 * 8] = { 0 }, h1[8 * 8] = { 0 }, v1[8] = { 0 }, r1[8 * 8] = { 0 };
	
	Map<MatrixXd>(x1, n, 1) = x;
	Map<MatrixXd>(p1, n, n) = P;
	Map<MatrixXd>(h1, n, m) = H;
	Map<MatrixXd>(v1, m, 1) = v;
	Map<MatrixXd>(r1, m, m) = R;
	
	MatrixXd K = P * H * ((H.transpose() * P * H + R).inverse());
	VectorXd xp = x + K * v;
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x.size(), x.size());
	MatrixXd Pp = (I - K * H.transpose()) * P;

	Map<MatrixXd>(xp_out, n, 1) = xp;
	Map<MatrixXd>(Pp_out, n, n) = Pp;

	return 0;
}

