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
	double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
	y = sin(arg) / arg;

	//y = pow(x(0), 2) + pow(x(1), 2);

	// zewnetrzna
	if ((*ad)(1) > 1)
	{
		if (-x(0) + 1 > 0) {
			y = y + (*ad)(0) * pow(-x(0) + 1, 2);
		}

		if (-x(1) + 1 > 0) {
			y = y + (*ad)(0) * pow(-x(1) + 1, 2);
		}

		if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - (*ud)(0) > 0)
			y = y + (*ad)(0) * pow(sqrt(pow(x(0), 2) + pow(x(1), 2)) - (*ud)(0), 2);
	}
	else {
		// wewnetrzna
		if (-x(0) + 1 > 0) {
			y = 1e+10;
		}
		else {
			y = y - (*ad)(0) / (-x(0) + 1);
		}

		if (-x(1) + 1 > 0) {
			y = 1e+10;
		}
		else {
			y = y - (*ad)(0) / (-x(1) + 1);
		}

		if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - (*ud)(0) > 0) {
			y = 1e+10;
		}
		else {
			y = y - (*ad)(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - (*ud)(0));
		}
	}
#elif LAB_NO==3 && LAB_PART==2
	matrix Y0 = matrix(4, new double[4]{ 0,x(0),100,0 });
	matrix* Y = solve_ode(0, 0.01, 7, Y0, &matrix(x(1)));
	int n = get_len(Y[0]);
	int i_50 = 0, i_0 = 0;
	for (int i = 0; i < n; ++i)
	{
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i_50, 2) - 50))
			i_50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i_0, 2)))
			i_0 = i;

		matrix b(1, 3);
		b(0, 0) = Y[0](i);
		b(0, 1) = Y[1](i, 0);
		b(0, 2) = Y[1](i, 2);
		(*ud).set_row(b, i);
	}
	y = -Y[1](i_0, 0);

	if (abs(x(0)) - 10 > 0)
		y = y + (*ad)(0) * pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 20 > 0)
		y = y + (*ad)(0) * pow(abs(x(1)) - 20, 2);
	if (abs(Y[1](i_50, 0) - 5) - 1 > 0)
		y = y + (*ad)(0) * pow(abs(Y[1](i_50, 0) - 5) - 1, 2);

	cout << "i_50: " << i_50 << " val x: " << Y[1](i_50, 0) << " val y: " << Y[1](i_50, 2)  << endl;
	cout << "i_0: " << i_0 << " val x: " << Y[1](i_0, 0) << " val y: " << Y[1](i_0, 2) << endl;
#elif LAB_NO==4 && (LAB_PART==1 || LAB_PART==2)
if (ad == nullptr) {
	y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
}
else
{
	solution tmp;
	tmp.x = ad[0] + x * ad[1];
	tmp.fit_fun(ud);
	y = tmp.y;
	--f_calls;
}
#elif LAB_NO==4 && LAB_PART==3
int m = 100, n = get_len(x);
static matrix X(n, m), Y(1, m);

if (f_calls == 1)
{
	ifstream S("XData.txt");
	S >> X;
	S.close();
	S.open("YData.txt");
	S >> Y;
	S.close();
}

double h;
double P = 0;
y = 0;

for (int i = 0; i < m; i++)
{
	h = (trans(x) * X[i])();
	h = 1 / (1 + exp(-h));
	y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);
	
	int h_normalized = h >= 0.5 ? 1 : 0;

	if (h_normalized == Y(0, i)) {
		P++;
	}
}

y = y / m;
P = P / m;
cout << "P: " << P << endl;
#elif LAB_NO==5 && LAB_PART==1
y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
#elif LAB_NO==5 && LAB_PART==2

#endif
}

void solution::grad(matrix *ud, matrix *ad)
{
	++g_calls;
#if LAB_NO==4 && (LAB_PART==1 || LAB_PART==2)
	g = matrix(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
#elif LAB_NO==4 && LAB_PART==3
	int m = 100, n = get_len(x);
	static matrix X(n, m), Y(1, m);
	if (g_calls == 1)
	{
		ifstream S("XData.txt");
		S >> X;
		S.close();
		S.open("YData.txt");
		S >> Y;
		S.close();
	}
	double h;
	g = matrix(n, 1);

	for (int j = 0; j < n; ++j)
	{
		for (int i = 0; i < m; i++)
		{
			h = (trans(x) * X[i])();
			h = 1 / (1 + exp(-h));
			g(j) = g(j) + X(j, i) * (h - Y(0, i));
		}
		g(j) = g(j) / m;
	}
#endif
}

void solution::hess(matrix *ud, matrix *ad)
{
	++H_calls;
#if LAB_NO==4 && (LAB_PART==1 || LAB_PART==2)
	H = matrix(2, 2);
	H(0, 0) = H(1, 1) = 10;
	H(0, 1) = H(1, 0) = 8;
#endif
}
