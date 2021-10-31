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
		double s = 0.1, epsilon = 1e-3, alpha_HJ = 0.5, aplha_R = 2, beta = 0.5;
		int Nmax = 5000;
		matrix x0 = 2 * rand_mat(2, 1) - 1, s0 = matrix(2, 1, s);
		cout << x0 << endl << endl;
		
		solution opt_HJ = HJ(x0, s, alpha_HJ, epsilon, Nmax);
		cout << opt_HJ << endl << endl;
		
		solution::clear_calls();

		solution opt_R = Rosen(x0, s0, aplha_R, beta, epsilon, Nmax);
		cout << opt_R << endl;
#elif LAB_NO==2 && LAB_PART==2

#elif LAB_NO==2 && LAB_PART==3

#elif LAB_NO==3 && LAB_PART==1

#elif LAB_NO==3 && LAB_PART==2

#elif LAB_NO==4 && LAB_PART==1

#elif LAB_NO==4 && LAB_PART==2

#elif LAB_NO==4 && LAB_PART==3

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
