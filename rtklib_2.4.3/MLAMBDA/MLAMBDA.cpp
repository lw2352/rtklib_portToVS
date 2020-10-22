// MLAMBDA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "ARLambda.h"
#include "MLAMBDA.h"
#include <iomanip>
#include <vector>
#include<fstream>
#include<sstream>
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
*          double *v        I   innovation (measurement - model) (m x 1)//新息
*          double *R        I   covariance matrix of measurement error (m x m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   states vector after update (n x 1)
*          double *Pp       O   covariance matrix of states after update (n x n)
* return : status (0:ok,<0:error)
* notes  : matirix stored by column-major order (fortran convention)
*          if state x[i]==0.0, not updates state x[i]/P[i+i*n]
*-----------------------------------------------------------------------------*/
int test_filter(double* x_in, double* P_in, double* H_in,
 double* v_in, double* R_in, int n, int m,
	double* xp_out, double* Pp_out)
{
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
//test
	double x1[80] = {0}, p1[80 * 80] = { 0 }, h1[80 * 80] = { 0 }, v1[80] = { 0 }, r1[80 * 80] = { 0 };
	
	Map<MatrixXd>(x1, n, 1) = x;
	Map<MatrixXd>(p1, n, n) = P;
	Map<MatrixXd>(h1, n, m) = H;
	Map<MatrixXd>(v1, m, 1) = v;
	Map<MatrixXd>(r1, m, m) = R;
//test end	
	MatrixXd K = P * H * ((H.transpose() * P * H + R).inverse());
	VectorXd xp = x + K * v;
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x.size(), x.size());
	MatrixXd Pp = (I - K * H.transpose()) * P;
	//https://zhuanlan.zhihu.com/p/159246989 计算机处理时可能由于截断误差导致Pp矩阵不对称。使用等效的计算公式Pp=(I-KH')P(I-KH')'+KRK' 更稳健
	//MatrixXd Pp = (I - K * H.transpose()) * P*(I - K * H.transpose()).transpose()+K*R*K.transpose();

	Map<MatrixXd>(xp_out, n, 1) = xp;
	Map<MatrixXd>(Pp_out, n, n) = Pp;

	return 0;
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

/*static double a = 6378137;//6378140;  //椭球的长半轴
//static double f = 0.00335281006247;
//static double b = a * (1 - f);
static double e2 = 0.0067394967422764;//0.006694384999588;(a * a - b * b) / pow(b,2);
static double m0 = a * (1 - e2);
static double m2 = 3.0 / 2 * e2 * m0;
static double m4 = 5.0 / 4 * e2 * m2;
static double m6 = 7.0 / 6 * e2 * m4;
static double m8 = 9.0 / 8 * e2 * m6;
static double a0 = m0 + m2 / 2 + (3.0 / 8.0) * m4 + (5.0 / 16.0) * m6 + (35.0 / 128.0) * m8;
static double a2 = m2 / 2 + m4 / 2 + 15.0 / 32 * m6 + 7.0 / 16 * m8;
static double a4 = m4 / 8 + 3.0 / 16 * m6 + 7.0 / 32 * m8;
static double a6 = m6 / 32 + m8 / 16;
static double a8 = m8 / 128;
static double xx = 0;
static double yy = 0;
static double _x = 0;
static double _y = 0;
static double BB = 0;
static double LL = 0;
static double PI = 3.1415926535897932;  /* pi */
/*static double F(double Bfi0)
{
	double ret = 0.0 - a2 * sin(2 * Bfi0) / 2.0 + a4 * sin(4 * Bfi0) / 4.0 - a6 * sin(6 * Bfi0) / 6.0;
	return ret;
}
static double hcfansuan(double pX)
{
	double Bf0 = pX / a0;
	double Bf1, Bf2;
	Bf1 = Bf0;
	Bf2 = (pX - F(Bf1)) / a0;
	while ((Bf2 - Bf1) > 1.0E-11)
	{
		Bf1 = Bf2;
		Bf2 = (pX - F(Bf1)) / a0;
	}
	return Bf1;
}
static void GaussNegative(double x, double y, double L0)
{
	double Bf, Vf, l, tf, hf2, Nf, Bmiao, Lmiao;
	int Bdu, Bfen, Ldu, Lfen;
	y = y - 500000;
	Bf = hcfansuan(x);
	Vf = sqrt(1 + e2 / (1 - e2) * cos(Bf) * cos(Bf));
	tf = tan(Bf);
	hf2 = e2 / (1 - e2) * cos(Bf) * cos(Bf);
	Nf = a / sqrt(1 - e2 * sin(Bf) * sin(Bf));
	BB = (Bf - 0.5 * Vf * Vf * tf * (pow(y / Nf, 2) - 1.0 / 12 * (5 + 3 * tf * tf + hf2 - 9 * hf2 * tf * tf) * pow(y / Nf, 4) + 1.0 / 360 * (61 + 90 * tf * tf + 45 * tf * tf) * pow(y / Nf, 6))) * 180 / PI;
	Bdu = (int)BB;
	Bfen = (int)((BB - Bdu) * 60);
	Bmiao = ((BB - Bdu) * 60 - Bfen) * 60;
	BB = Bdu + 0.01 * Bfen + 0.0001 * Bmiao;
	l = 1.0 / cos(Bf) * (y / Nf - 1.0 / 6.0 * (1 + 2 * tf * tf + hf2) * pow(y / Nf, 3) + 1.0 / 120.0 * (5 + 28 * tf * tf + 24 * pow(tf, 4) + 6 * hf2 + 8 * hf2 * tf * tf) * pow(y / Nf, 5)) * 180.0 / PI;
	LL = L0 + l;
	Ldu = (int)LL;
	Lfen = (int)((LL - Ldu) * 60);
	Lmiao = ((LL - Ldu) * 60 - Lfen) * 60;
	LL = Ldu + 0.01 * Lfen + 0.0001 * Lmiao;
}*/

//另一种程序
	 double a;//'椭球体长半轴
	 double b;// '椭球体短半轴
	 double f; //'扁率
	 double e;// '第一偏心率
	 double e1; //'第二偏心率
	 double FE;//'东偏移
	 double FN;//'北偏移
	 double L0;//'中央经度
	 double W0;//'原点纬线
	 double k0;//'比例因子
	 double PI = 3.1415926535897932;
	/**
	 * 幂函数
	 * @param e
	 * @param i
	 * @return
	 */
 double MZ(double e, int i)
	{
		return pow(e, i);
	}

 void init(int TuoqiuCanshu, double CentralMeridian, double OriginLatitude, double EastOffset, double NorthOffset)
 {

			a = 6378137;
			b = 6356752.3142;
		

		f = (a - b) / a;//扁率
						//e = sqrt(1 - MZ((b / a) ,2));//'第一偏心率
		e = sqrt(2 * f - MZ(f, 2));//'第一偏心率
										//eq = sqrt(MZ((a / b) , 2) - 1);//'第二偏心率
		e1 = e / sqrt(1 - MZ(e, 2));//'第二偏心率
		L0 = CentralMeridian;//中央经
		W0 = OriginLatitude;//原点纬线
		k0 = 1;//'比例因子
		FE = EastOffset;//东偏移
		FN = NorthOffset;//北偏移
	}

	/**
	 * 高斯投影坐标 转为 经纬度坐标
	 * @param X 高斯投影坐标X
	 * @param Y 高斯投影坐标Y
	 * @return
	 */
 double resultP[2];// = new double[2];
	void GKgetJW(double X, double Y)
	{
		//'给出高克投影坐标，转换为经纬度坐标
		
		double El1 = (1 - sqrt(1 - MZ(e, 2))) / (1 + sqrt(1 - MZ(e, 2)));
		double Mf = (Y - FN) / k0;//真实坐标值
		double Q = Mf / (a * (1 - MZ(e, 2) / 4 - 3 * MZ(e, 4) / 64 - 5 * MZ(e, 6) / 256));//角度
		double Bf = Q + (3 * El1 / 2 - 27 * MZ(El1, 3) / 32) * sin(2 * Q) + (21 * MZ(El1, 2) / 16 - 55 * MZ(El1, 4) / 32) * sin(4 * Q) + (151 * MZ(El1, 3) / 96) * sin(6 * Q) + 1097 / 512 * MZ(El1, 4) * sin(8 * Q);
		double Rf = a * (1 - MZ(e, 2)) / sqrt(MZ((1 - MZ((e * sin(Bf)), 2)), 3));
		double Nf = a / sqrt(1 - MZ((e * sin(Bf)), 2));//卯酉圈曲率半径
		double Tf = MZ((tan(Bf)), 2);
		double D = (X - FE) / (k0 * Nf);
		double Cf = MZ(e1, 2) * MZ((cos(Bf)), 2);
		double B = Bf - Nf * tan(Bf) / Rf * (MZ(D, 2) / 2 - (5 + 3 * Tf + 10 * Cf - 9 * Tf * Cf - 4 * MZ(Cf, 2) - 9 * MZ(e1, 2)) * MZ(D, 4) / 24 + (61 + 90 * Tf + 45 * MZ(Tf, 2) - 256 * MZ(e1, 2) - 3 * MZ(Cf, 2)) * MZ(D, 6) / 720);
		double L = L0 * PI / 180 + 1 / cos(Bf) * (D - (1 + 2 * Tf + Cf) * MZ(D, 3) / 6 + (5 - 2 * Cf + 28 * Tf - 3 * MZ(Cf, 2) + 8 * MZ(e1, 2) + 24 * MZ(Tf, 2)) * MZ(D, 5) / 120);
		double Bangle = B * 180 / PI;
		double Langle = L * 180 / PI;
		resultP[0] = Langle;//RW * 180 / Math.PI;
		resultP[1] = Bangle + W0;//RJ * 180 / Math.PI;
		//return resultP;
	}

void test()
{
	ifstream fin("C:\\Users\\204\\Desktop\\data.txt");
	string line_info;
	ofstream fout("result.txt");
	init(2, 102.5, 0, 500000, 0);
	if (fin) // 有该文件
	{
		while (getline(fin, line_info)) // line中不包括每行的换行符
		{
			string input_result;
			vector<string> vectorString;
			stringstream input(line_info);
			//依次输出到input_result中，并存入vectorString中
			//cout << "line_info: " << line_info << endl;
			while (input >> input_result)
				vectorString.push_back(input_result);//分段保存结果
			if (vectorString.size() >= 3)
			{
				//调用转换函数
				double x = transformData(vectorString[2]), y = transformData(vectorString[3]);
				//x = 3137245.1060, y = 498899.0800;
				//GaussNegative(x, y, 102.5);
				GKgetJW(y, x);
				//打印结果
				cout << setiosflags(ios::fixed) << setprecision(8) << resultP[0]<< ", " << resultP[1] << endl;
				//输出结果保存到文件
				fout << setiosflags(ios::fixed) << setprecision(8)  << resultP[0] << ", " << resultP[1] << endl;

			}
		}
		fout.close();
	}
	else // 没有该文件
	{
		cout << "no such file" << endl;
	}
}





