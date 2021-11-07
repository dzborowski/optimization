//Do not edit the code below (unless you know what you are doing)

#include"solution.h"

int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

solution::solution(double L)
{
	x = L;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(const matrix &A)
{
	x = A;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(int n, double *A)
{
	x = matrix(n, A);
	g = NAN;
	H = NAN;
	y = NAN;
}

int get_dim(const solution &A)
{
	return get_len(A.x);
}

void solution::clear_calls()
{
	f_calls = 0;
	g_calls = 0;
	H_calls = 0;
}

ostream &operator<<(ostream &S, const solution &A)
{
	S << "x = " << A.x << endl;
	S << "y = " << A.y << endl;
	S << "f_calls = " << solution::f_calls << endl;
	if (solution::g_calls>0)
		S << "g_calls = " << solution::g_calls << endl;
	if (solution::H_calls)
		S << "H_calls = " << solution::H_calls << endl;
	return S;
}

//You can edit the following code

void solution::fit_fun(matrix *ud, matrix *ad)
{
	++f_calls;
#if LAB_NO==1 && (LAB_PART==1 || LAB_PART==2)
	y = -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(), 2);
#elif LAB_NO==1 && LAB_PART==3
	matrix Y0 = matrix(3, new double[3]{5, 1, 10}); // Va, Vb, Tb
	matrix *Y = solve_ode(0, 1, 1000, Y0, ud, &x); // od 0, krok 1, 1000 sekund symulacja
	// Y[0] wketor czasu, Y[1], 3 kolumny Va, Vb, Tb, 1001 wpisów
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);

	for (int i = 1; i < n; i++)
	{		
		if (max < Y[1](i, 2)) {
			max = Y[1](i, 2);
		}

		if (abs(x() - 0.00241025) < 1e-5)
		{
			matrix b(1, 3);
			b(0, 0) = Y[1](i, 0);
			b(0, 1) = Y[1](i, 1);
			b(0, 2) = Y[1](i, 2);
			(*ud).set_row(b, i);
		}
	}
	y = abs(max - 50);
#elif LAB_NO==2 && (LAB_PART==1 || LAB_PART==2)
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
#elif LAB_NO==2 && LAB_PART==3
	double a_ref = 3.14, w_ref = 0;
	matrix Y0(2, 1);
	matrix* Y = solve_ode(0, 0.1, 100, Y0, ud, &x);
	int n = get_len(Y[0]);
	y = 0;

	for (int i = 0; i < n; i++)
	{
		y = y + 10 * pow(a_ref - Y[1](i, 0), 2) + pow(w_ref - Y[1](i, 1), 2) + pow(x(0) * (a_ref - Y[1](i, 0)) + x(1) * (w_ref - Y[1](i, 1)), 2);
		matrix b(1, 2);
		b(0, 0) = Y[1](i, 0);
		b(0, 1) = Y[1](i, 1);
		(*ud).set_row(b, i);
	}
	y = 0.1 * y;

#elif LAB_NO==3 && LAB_PART==1

#elif LAB_NO==3 && LAB_PART==2

#elif LAB_NO==4 && (LAB_PART==1 || LAB_PART==2)

#elif LAB_NO==4 && LAB_PART==3

#elif LAB_NO==5 && LAB_PART==1

#elif LAB_NO==5 && LAB_PART==2

#endif
}

void solution::grad(matrix *ud, matrix *ad)
{
	++g_calls;
#if LAB_NO==4 && (LAB_PART==1 || LAB_PART==2)

#elif LAB_NO==4 && LAB_PART==3

#endif
}

void solution::hess(matrix *ud, matrix *ad)
{
	++H_calls;
#if LAB_NO==4 && (LAB_PART==1 || LAB_PART==2)

#endif
}
