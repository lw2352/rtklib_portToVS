// MLAMBDA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "ARLambda.h"
#include "MLAMBDA.h"
#include <iomanip>
#include <vector>
#include<fstream>
#include<sstream>
#include <map>
#include <omp.h>

using namespace Eigen;
using namespace std;
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

//时间更新 x=F*x, P=F*P*F+Q(Q在本函数之后的语句中加入)
void test_udPos(int n, double* F_in, double* x_in, double* P_in, double* xp_out, double* Pp_out)
{
	MatrixXd F = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(F_in, n, n);
	MatrixXd x = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(x_in, n, 1);
	MatrixXd P = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(P_in, n, n);

	double s = 0.995;
	MatrixXd xp = F*x;
	MatrixXd Pp = F*((1.0/s)*P)*F.transpose();

	Map<MatrixXd>(xp_out, n, 1) = xp;
	Map<MatrixXd>(Pp_out, n, n) = Pp;
}


/* kalman filter ---------------------------------------------------------------
* kalman filter state update as follows:
*
*   K=P*H*(H'*P*H+R)^-1, xp=x+K*v, Pp=(I-K*H')*P
*
* args   : double *x        I   states vector (n x 1)
*          double *P        I   covariance matrix of states (n x n)
*          double *H        I   transpose of design matrix (n x m)
*          double *v        I   innovation (measurement - model) (m x 1)//新息
*          double *R        I   covariance matrix of measurement error (m x m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   states vector after update (n x 1)
*          double *Pp       O   covariance matrix of states after update (n x n)
* return : status (0:ok,<0:error)
* notes  : matirix stored by column-major order (fortran convention)
*          if state x[i]==0.0, not updates state x[i]/P[i+i*n]
*-----------------------------------------------------------------------------*/
//状态更新
int test_filter(double* x_in, double* P_in, double* H_in,
 double* v_in, double* R_in, int n, int m,
	double* xp_out, double* Pp_out)
{
	double t1, t2;
	t1 = omp_get_wtime();
	/* m是参与解算的卫星颗数
	* x=(r,v,B),n是x矩阵的大小，n=3+3+8=14;
	*/
	if (x_in[0] == 0.0)
		return -1;
	VectorXd x= Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(x_in, n, 1);
	MatrixXd P = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(P_in, n, n);
	MatrixXd H = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(H_in, n, m);
	VectorXd v = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(v_in, m, 1);
	MatrixXd R = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(R_in, m, m);
#if 0
	cout << "x=\n"<<x.transpose() << "\n";
	cout << "P=\n" << P.transpose() << "\n";
	cout << "H=\n" << H.transpose() << "\n";
	cout << "v=\n" << v.transpose() << "\n";
	cout << "R=\n" << R.transpose() << "\n";
#endif
	MatrixXd K = P * H * ((H.transpose() * P * H + R).inverse());
	VectorXd xp = x + K * v;
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x.size(), x.size());
	MatrixXd Pp = (I - K * H.transpose()) * P;

	Map<MatrixXd>(xp_out, n, 1) = xp;
	Map<MatrixXd>(Pp_out, n, n) = Pp;
	t2 = omp_get_wtime();
	//printf("time use: %lf\n", (t2 - t1));
	return 0;
}
//这里的A已经是转置的了，和公式有点不一样
int test_lsq(double* A_in, double* y_in, int n, int m, double* x_in,
	double* Q_in)
{
	//cout << "start test_lsq\n";
	if (m < n) return -1;
	VectorXd x = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(x_in, n, 1);
	MatrixXd Q = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(Q_in, n, n);
	MatrixXd A = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(A_in, n, m);
 	VectorXd y = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(y_in, m, 1);

	Q = (A * A.transpose()).inverse();
	x = Q * A * y;
	/*cout << x << "\n\n";
	cout << Q << "\n\n";
	cout << A << "\n\n";
	cout << y << "\n\n";*/
	Map<MatrixXd>(x_in, n, 1) = x;	
	Map<MatrixXd>(Q_in, n, n) = Q;
	return 0;//ok
}

void testVel(double* Ir_in, double* Ib_in, double* Vs_in, double* Fr_in, double* Fb_in, double* x_in, int n,double lam)
{
	//cout << "start testVel\n";
	VectorXd x = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(x_in, 4, 1);
	MatrixXd Ir = Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(Ir_in, n, 3);
	MatrixXd Ib = Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(Ib_in, n, 3);
	MatrixXd Vs = Map<Matrix<double, Dynamic, Dynamic, RowMajor> >(Vs_in, n, 3);
	VectorXd Fr = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(Fr_in, n, 1);
	VectorXd Fb = Map<Matrix<double, Dynamic, Dynamic, ColMajor> >(Fb_in, n, 1);
#if 0	
	cout << "Ir=\n"<< Ir << "\n\n";

	cout << "Ib=\n" << Ib << "\n\n";

	cout << "Vs=\n" << Vs << "\n\n";
	cout << "Fr=\n" << Fr << "\n\n";
	cout << "Fb=\n" << Fb << "\n\n";
#endif

	VectorXd y = ((Ir - Ib) * (Vs.transpose())).diagonal() + (Fr - Fb)*lam;
	//cout << y << "\n\n";
	MatrixXd Q = MatrixXd::Zero(4, 4);
	MatrixXd A = Ir;
	A.conservativeResize(A.rows(), A.cols() + 1);
	VectorXd I= VectorXd::Ones(n, 1);
	A.col(A.cols()-1) = I;
	//cout << A << "\n\n";
#if 0
	Q = (A.transpose() * A).inverse();
	//cout << Q << "\n\n";
	x = Q * A.transpose() * y;
#else
	x = A.colPivHouseholderQr().solve(y);
#endif
	cout <<"x="<< x.transpose() << "\n";
	Map<MatrixXd>(x_in, 4, 1) = x;

}
/*int resultFilter(double* x, double* y, double* z)
{
	static vector<double> X, Y, Z;
	//X.resize(N);
	//Y.resize(N);
	//Z.resize(N);
	if (X.size() < N)
	{
		X.push_back(*x);
		Y.push_back(*y);
		Z.push_back(*z);
	}
	else
	{
		sort(X.begin(), X.end());
		sort(Y.begin(), Y.end());
		sort(Z.begin(), Z.end());
		
		X.erase(X.begin());
		X.erase(X.end() - 1);
		Y.erase(Y.begin());
		Y.erase(Y.end()-1);
		Z.erase(Z.begin());
		Z.erase(Z.end()-1);
		int n = (int)N - 2;
		*x = *y = *z = 0;
		for (int i = 0; i < n; i++)
		{
			*x += X[i]/n;
			*y += Y[i]/n;
			*z += Z[i]/n;
		}
		X.clear();
		Y.clear();
		Z.clear();
		return 1;
	}
	return 0;
}*/

int resultFilter(double* x, double* y, double* z, int n)
{
	static vector<double> X, Y, Z;
	//添加数据
	//if (X.size() < n)
	//{
		X.push_back(*x);
		Y.push_back(*y);
		Z.push_back(*z);
	//}
	//数据小于窗口大小时输出原始数据
	if (X.size() < n)
	{
		return 1;
	}
	else if(X.size() >= n)
	{
		//cout << setprecision(15)<< X[0] << "  " << setprecision(15) << X[1] << "  " << setprecision(15) << X[2] << "  " << endl;
		*x = *y = *z = 0;
		int i = X.size() - 1;
		int j = 0;
		//对最新的n个数据求平均值
		do
		{
			*x += X[i] / n;
			*y += Y[i] / n;
			*z += Z[i] / n;
			i--;
			j++;
		} while (j<n);

		X.erase(X.begin());
		Y.erase(Y.begin());
		Z.erase(Z.begin());

		return 1;
	}
}

double transformData(string data)
{
	char ret[16];
	int n = data.size();
	int i = 0,j=0;
	while (j <n)
	{
		if (data[j] >= '0' && data[j] <= '9' || data[j] == '.')
		{
			ret[i] = data[j];
			i++;
		}
		j++;
	}
	return atof(ret);
}

int resultSTD(double* rr, int n,double *out)
{
	static vector<double> X, Y, Z;
	double sumX=0, sumY=0, sumZ=0;
	//添加数据
	//if (X.size() < n)
	//{
	X.push_back(rr[0]);
	Y.push_back(rr[1]);
	Z.push_back(rr[2]);
	//}
	//数据小于窗口大小时输出原始数据
	if (X.size() < n)
	{
		return 1;
	}
	else if (X.size() >= n)
	{
		//cout << setprecision(15)<< X[0] << "  " << setprecision(15) << X[1] << "  " << setprecision(15) << X[2] << "  " << endl;
		int i = X.size() - 1;
		int j = 0;
		//对最新的n个数据求平均值
		do
		{
			//*x += X[i] / n;
			//*y += Y[i] / n;
			//*z += Z[i] / n;
			sumX += X[i];
			sumY += Y[i];
			sumZ += Z[i];
			i--;
			j++;
		} while (j < n);
		//求均值
		double meanX = sumX / n;
		double meanY = sumY / n;
		double meanZ = sumZ / n;
		double accumX = 0.0, accumY = 0.0, accumZ= 0.0;
		for(int j=0,i= X.size() - 1; j<n; j++,i--)
		{
			accumX += (X[i] - meanX) * (X[i] - meanX);
			accumY += (Y[i] - meanY) * (Y[i] - meanY);
			accumZ += (Z[i] - meanZ) * (Z[i] - meanZ);
		}
		//标准差
		out[0] = sqrt(accumX / (double)(n));
		out[1] = sqrt(accumY / (double)(n));
		out[2] = sqrt(accumZ / (double)(n));
		X.erase(X.begin());
		Y.erase(Y.begin());
		Z.erase(Z.begin());

		return 1;
	}
}

class OBS
{
public:
	string P, L, D, S;
	OBS(string p,string l,string d,string s);
};
OBS::OBS(string p, string l, string d, string s)
{
	P = p; L = l; D = d; S = s;
}
void getObs()
{
	ifstream fin("D:\\data\\3-19\\7.obs");
	string line_info;
	ofstream fout("D:\\data\\3-19\\7Result.txt");
	vector<OBS> vecObsG10;
	vector<OBS> vecObsG22;
	vector<OBS> vecObsG25;
	vector<OBS> vecObsG26;
	vector<OBS> vecObsG29;
	vector<OBS> vecObsG31;
	vector<OBS> vecObsG32;
	
	if (fin) // 有该文件
	{
		while (getline(fin, line_info)) // line中不包括每行的换行符
		{
			string input_result;
			vector<string> vecString;
			stringstream input(line_info);
			//依次输出到input_result中，并存入vectorString中
			//cout << "line_info: " << line_info << endl;
			while (input >> input_result)
				vecString.push_back(input_result);//分段保存结果
			
			if (vecString[0] == "C12")
			{
				OBS obs(vecString[1], vecString[2], vecString[3], vecString[4]);
				vecObsG10.push_back(obs);
			}
			if (vecString[0] == "C10")
			{
				OBS obs(vecString[1], vecString[2], vecString[3], vecString[4]);
				vecObsG22.push_back(obs);
			}
			if (vecString[0] == "G13")
			{
				OBS obs(vecString[1], vecString[2], vecString[3], vecString[4]);
				vecObsG25.push_back(obs);
			}
			if (vecString[0] == "G12")
			{
				OBS obs(vecString[1], vecString[2], vecString[3], vecString[4]);
				vecObsG26.push_back(obs);
			}
			if (vecString[0] == "C21")
			{
				OBS obs(vecString[1], vecString[2], vecString[3], vecString[4]);
				vecObsG29.push_back(obs);
			}
			if (vecString[0] == "G17")
			{
				OBS obs(vecString[1], vecString[2], vecString[3], vecString[4]);
				vecObsG31.push_back(obs);
			}
			if (vecString[0] == "G19")
			{
				OBS obs(vecString[1], vecString[2], vecString[3], vecString[4]);
				vecObsG32.push_back(obs);
			}
		}
		//写入文件
		for (int i = 0; i < vecObsG10.size(); i++)
		{
			string str = vecObsG10[i].D + "," + vecObsG22[i].D + ","+vecObsG25[i].D + ","+vecObsG26[i].D + "," + vecObsG29[i].D + "," + vecObsG31[i].D + "," + vecObsG32[i].D;
			fout<<str<<endl;
		}
		fout.close();
	}
	else // 没有该文件
	{
		cout << "no such file" << endl;
	}
}
VectorXd b = VectorXd::Random(100);
VectorXd x(100);
MatrixXd A = MatrixXd::Random(100, 100);
void test()
{
	//int dim = 100;
	
	//cout << "A =" << endl << A << endl;
	//MatrixXd B = MatrixXd::Random(dim, dim);
	///*
	//x = A.inverse() * b;
	x = A.colPivHouseholderQr().solve(b);
	//cout << "The solution is:\n" << x << endl;
	//*/
	//A = A * B;
}




