/***************************************************
Code written for the optimization exercises purposes
by Lukasz Sztangret, PhD
Department of Applied Computer Science and Modelling
AGH University of Science and Technology
***************************************************/

#include <iostream>
#include <random>
#include"opt_alg.h"
#include"solution.h"

int main()
{
	try
	{
		cout << "LAB NUMBER " << LAB_NO << endl;
		cout << "LAB PART " << LAB_PART << endl << endl;
#if LAB_NO==0

#elif LAB_NO==1 && LAB_PART==1
		std::uniform_real_distribution<double> unif(-100, 100);
		std::default_random_engine re;
		double alphas[3] = { 2, 5, 10};

		for (int i = 0; i < 3; i++)
		{
			double alpha = alphas[i];
			double output[100][12];
			
			for (int j = 0; j < 100; j++)
			{
				double x0 = unif(re);

				double d = 100, epsilon = 1e-5, gamma = 1e-200; // gamma = d - d_old w lagrange
				int N_max = 1000;

				solution::clear_calls();
				double* p = expansion(x0, d, alpha, N_max);

				solution::clear_calls();
				solution opt_fib = fib(p[0], p[1], epsilon);
				int opt_fib_f_calls = solution::f_calls;

				solution::clear_calls();
				solution opt_lag = lag(p[0], p[1], epsilon, gamma, N_max);
				int opt_lag_f_calls = solution::f_calls;

				output[j][0] = x0;
				output[j][1] = p[0];
				output[j][2] = p[1];
				output[j][3] = p[2];
				output[j][4] = opt_fib.x[0]();
				output[j][5] = opt_fib.y[0]();
				output[j][6] = opt_fib_f_calls;
				output[j][7] = opt_fib.x[0]() >= 62 && opt_fib.x[0]() <= 63 ? 1.0 : 0.0;
				output[j][8] = opt_lag.x[0]();
				output[j][9] = opt_lag.y[0]();
				output[j][10] = opt_lag_f_calls;
				output[j][11] = opt_lag.x[0]() >= 62 && opt_lag.x[0]() <= 63 ? 1.0 : 0.0;;
			}

			std::ofstream out(std::to_string(alpha) + "_lab_1_part_1.csv");

			for (auto& row : output) {
				for (auto col : row)
					out << col << ',';
				out << '\n';
			}
		}
		

		/*
		double output[100][12];
		
		double x0 = 62, d = 70, alpha = 0.5, epsilon = 1e-5, gamma = 1e-200; // gamma = d - d_old w lagrange
		int N_max = 1000;
		double *p = expansion(x0, d, alpha, N_max);

		std::cout << "p[0]: " << p[0] << std::endl;
		std::cout << "p[1]: " <<  p[1] << std::endl;
		std::cout << "p[2]: " <<  p[2] << std::endl;

		solution::clear_calls();
		solution opt_fib = fib(p[0], p[1], epsilon);
		int opt_fib_f_calls = solution::f_calls;
		std::cout << opt_fib << std::endl;

		solution::clear_calls();
		solution opt_lag = lag(p[0], p[1], epsilon, gamma, N_max);
		int opt_lag_f_calls = solution::f_calls;
		std::cout << opt_lag << std::endl;


		output[0][0] = x0;
		output[0][1] = p[0];
		output[0][2] = p[1];
		output[0][3] = p[2];
		output[0][4] = opt_fib.x[0]();
		output[0][5] = opt_fib.y[0]();
		output[0][6] = opt_fib_f_calls;
		output[0][7] = 1.0;
		output[0][8] = opt_lag.x[0]();
		output[0][9] = opt_lag.y[0]();
		output[0][10] = opt_lag_f_calls;
		output[0][11] = 1.0;

		std::ofstream out(std::to_string(alpha) + "_lab_1_part_1.csv");

		for (auto& row : output) {
			for (auto col : row)
				out << col << ',';
			out << '\n';
		}
		*/
		
#elif LAB_NO==1 && LAB_PART==2
		double epsilon = 1e-5, gamma = 1e-200;
		int N_max = 1000;

		matrix ab_F(1, 1, 200);
		solution::clear_calls();
		solution opt_f = fib(-100, 100, epsilon, &ab_F);
		std::cout << opt_f << std::endl;
		std::cout << solution::f_calls << std::endl;
		std::cout << ab_F << std::endl;
		
		matrix ab_L(1, 1, 200);
		solution::clear_calls();
		solution opt_l = lag(-100, 100, epsilon, gamma, N_max, &ab_L);
		std::cout << opt_l << std::endl;
		std::cout << ab_L << std::endl;

#elif LAB_NO==1 && LAB_PART==3
		double epsilon = 1e-5, gamma = 1e-200;
		int N_max = 1000;
		
		int row = 1001;
		int col = 3;

		matrix ab_F(row, col);
		solution::clear_calls();
		solution opt_f = fib(0.0001, 0.01, epsilon, &ab_F);
		std::cout << opt_f << std::endl;
		
		std::ofstream out_f("lab_1_part_3_fib.csv");

		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
				if (j < (col - 1)) {
					out_f << ab_F(i, j) << ",";
				}
				else if (j == (col - 1)) {
					out_f << ab_F(i, j) << "\n";
				}
		}

		matrix ab_L(row, col);
		solution::clear_calls();
		solution opt_l = lag(0.0001, 0.01, epsilon, gamma, N_max, &ab_L);
		std::cout << opt_l << std::endl;

		std::ofstream out_l("lab_1_part_3_lag.csv");

		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
				if (j < (col - 1)) {
					out_l << ab_L(i, j) << ",";
				}
				else if (j == (col - 1)) {
					out_l << ab_L(i, j) << "\n";
				}
		}
		
#elif LAB_NO==2 && LAB_PART==1
		double epsilon = 1e-3, alpha_HJ = 0.5, aplha_R = 2, beta = 0.5;
		int Nmax = 5000;

		double steps[3] = { 0.1, 0.5, 1.0 };
		auto isGlobalMin = [epsilon](double x1, double x2) {
			x1 = abs(x1);
			x2 = abs(x2);
			return (x1 >= 0 && x1 < epsilon) && (x2 >= 0 && x2 < epsilon);
		};

		for (int i = 0; i < 3; i++)
		{
			double s = steps[i];
			matrix s0 = matrix(2, 1, s);
			double output[100][12];

			for (int j = 0; j < 100; j++)
			{
				matrix x0 = 2 * rand_mat(2, 1) - 1;

				solution::clear_calls();
				solution opt_HJ = HJ(x0, s, alpha_HJ, epsilon, Nmax);
				int opt_hj_f_calls = solution::f_calls;

				solution::clear_calls();
				solution opt_R = Rosen(x0, s0, aplha_R, beta, epsilon, Nmax);
				int opt_rosen_f_calls = solution::f_calls;

				output[j][0] = x0(0);
				output[j][1] = x0(1);
				output[j][2] = opt_HJ.x(0);
				output[j][3] = opt_HJ.x(1);
				output[j][4] = opt_HJ.y(0);
				output[j][5] = opt_hj_f_calls;
				output[j][6] = isGlobalMin(opt_HJ.x(0), opt_HJ.x(1)) ? 1.0 : 0;
				output[j][7] = opt_R.x(0);
				output[j][8] = opt_R.x(1);
				output[j][9] = opt_R.y(0);
				output[j][10] = opt_rosen_f_calls;
				output[j][11] = isGlobalMin(opt_R.x(0), opt_R.x(1)) ? 1.0 : 0;
			}

			std::ofstream out(std::to_string(s) + "_lab_2_part_1.csv");

			for (auto& row : output) {
				for (auto col : row)
					out << col << ',';
				out << '\n';
			}
		}
		
#elif LAB_NO==2 && LAB_PART==2
		double s = 0.1, epsilon = 1e-3, alpha_HJ = 0.5, aplha_R = 2, beta = 0.5;
		int Nmax = 5000;
		matrix s0 = matrix(2, 1, s);
		matrix x0 = matrix(2, 1);
		x0(0) = -0.0432063;
		x0(1) = -0.343887;
		cout << x0 << endl << endl;
		
		solution::clear_calls();
		matrix ud_HJ(1, 2);
		solution opt_HJ = HJ(x0, s, alpha_HJ, epsilon, Nmax, &ud_HJ);
		cout << opt_HJ << endl << endl;
		cout << ud_HJ << endl << endl;
		
		solution::clear_calls();
		matrix ud_R(1, 2);
		solution opt_R = Rosen(x0, s0, aplha_R, beta, epsilon, Nmax, &ud_R);
		cout << opt_R << endl;
		cout << ud_R << endl;

#elif LAB_NO==2 && LAB_PART==3
		double s = 0.1, epsilon = 1e-3, alpha_HJ = 0.5, aplha_R = 2, beta = 0.5;
		int Nmax = 5000;
		int row = 1001;
		int col = 2;
		matrix s0 = matrix(2, 1, s);
		matrix x0 = matrix(2, new double[2]{ 1, 1 });
		cout << x0 << endl << endl;

		solution::clear_calls();
		matrix ud_HJ(row, col);
		solution opt_HJ = HJ(x0, s, alpha_HJ, epsilon, Nmax, &ud_HJ);
		cout << opt_HJ << endl << endl;

		std::ofstream out_hj("lab_2_part_3_hj.csv");

		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
				if (j < (col - 1)) {
					out_hj << ud_HJ(i, j) << ",";
				}
				else if (j == (col - 1)) {
					out_hj << ud_HJ(i, j) << "\n";
				}
		}

		solution::clear_calls();
		matrix ud_R(row, col);
		solution opt_R = Rosen(x0, s0, aplha_R, beta, epsilon, Nmax, &ud_R);
		cout << opt_R << endl;

		std::ofstream out_r("lab_2_part_3_r.csv");

		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
				if (j < (col - 1)) {
					out_r << ud_R(i, j) << ",";
				}
				else if (j == (col - 1)) {
					out_r << ud_R(i, j) << "\n";
				}
		}
#elif LAB_NO==3 && LAB_PART==1
matrix x0, a = 4; // a znajduje siê w g3

// zrobiæ losowanie tak, ¿eby by³o w obaszarze dopuszczalnym
do {
	x0 = 4 * rand_mat(2, 1) + 1;
} while (norm(x0) > a(0));

cout << x0 << endl << endl;
//x0.add_row(10);
//cout << x0(1, 0) << endl << endl;
//cout << x0(2, 0) << endl << endl;

////double c0 = 1, dc = 2, epsilon = 1e-5; // zewnetrzna
//double c0 = 10, dc = 0.5, epsilon = 1e-5; // wewnetrzna
//int Nmax = 10000;
//
//// 4,49
//
//solution opt_zew = pen(x0, c0, dc, epsilon, Nmax, &a);
//cout << opt_zew<< endl << endl;
//cout << sqrt(pow(opt_zew.x(0), 2) + pow(opt_zew.x(1), 2))  << endl << endl; // 4.000001

#elif LAB_NO==3 && LAB_PART==2
//matrix x0(2, 1, 2), c = 1;
//solution test(x0);
//test.fit_fun(nullptr, &c);
//cout << test << endl;

double c0 = 1, dc = 2, epsilon = 1e-1; // zewnetrzna
int Nmax = 2000;

std::uniform_real_distribution<double> v_unif(-10, 10);
std::uniform_real_distribution<double> w_unif(-20, 20);
std::default_random_engine v_re;
std::default_random_engine w_re;

matrix x0(2, 1), ud(1, 3);
x0(0, 0) = v_unif(v_re);
x0(1, 0) = w_unif(w_re);

cout << x0 << endl << endl;

solution opt = pen(x0, c0, dc, epsilon, Nmax, &ud);
cout << opt << endl << endl;

#elif LAB_NO==4 && LAB_PART==1
matrix x0 = 20 * rand_mat(2, 1) - 10;
double epsilon = 1e-3, h0 = -0.05;
int Nmax = 10000;

solution opt;
opt = SD(x0, h0, epsilon, Nmax);
cout << opt << endl;
solution::clear_calls();


opt = CG(x0, h0, epsilon, Nmax);
cout << opt << endl;
solution::clear_calls();

opt = Newton(x0, h0, epsilon, Nmax);
cout << opt << endl;
solution::clear_calls();
#elif LAB_NO==4 && LAB_PART==2

#elif LAB_NO==4 && LAB_PART==3
matrix x0(3, new double[3]{ -1, 0.1, 0.1 });
solution test(x0);
test.fit_fun();
test.grad();

cout << test << endl;
cout << test.g << endl;
#elif LAB_NO==5 && LAB_PART==1

#elif LAB_NO==5 && LAB_PART==2

#endif
	}
	catch (char * EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
