#include"opt_alg.h"
#if LAB_NO>0
double *expansion(double x0, double d, double alpha, int Nmax, matrix *ud, matrix *ad)
{
	double *p = new double[3];
	solution X0(x0), X1(x0 + d);
	X0.fit_fun(ud, ad);
	X1.fit_fun(ud, ad);

	if (X0.y == X1.y)
	{
		p[0] = X0.x();
		p[1] = X1.x();
		p[2] = solution::f_calls;
		return p;
	}

	if (X0.y < X1.y)
	{
		d *= -1;
		X1.x = X0.x + d;
		X1.fit_fun(ud, ad);

		if (X0.y <= X1.y)
		{
			p[0] = X1.x();
			p[1] = X0.x() - d;
			p[2] = solution::f_calls;
			return p;
		}
	}

	solution X2;
	int i = 1;
	
	while (true)
	{
		X2.x = x0 + pow(alpha, i) * d;
		X2.fit_fun(ud, ad);

		if (X2.y >= X1.y || solution::f_calls > Nmax)
			break;

		X0 = X1;
		X1 = X2;
		++i;
	}
	d > 0 ? p[0] = X0.x(), p[1] = X2.x() : (p[0] = X2.x(), p[1] = X0.x());
	p[2] = solution::f_calls;
	return p;
}

solution fib(double a, double b, double epsilon, matrix *ud, matrix *ad)
{
	int n = static_cast<int>(ceil(log2(sqrt(5) * (b-a) / epsilon) / log2((1+ sqrt(5)) / 2)));
	int *F = new int[n] {1, 1};

	for (int i = 2; i < n; ++i)
		F[i] = F[i - 2] + F[i - 1];
	
	solution A(a), B(b), C, D;
	C.x = B.x - 1.0 * F[n-2] / F[n-1] * (B.x - A.x);
	D.x = A.x + B.x - C.x;
	C.fit_fun(ud, ad);
	D.fit_fun(ud, ad);

	for (int i = 0; i <= n - 3; ++i)
	{
		if (C.y < D.y)
			B = D;
		else
			A = C;
		C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);;
		D.x = A.x + B.x - C.x;
		C.fit_fun(ud, ad);
		D.fit_fun(ud, ad);
#if LAB_NO==1 && LAB_PART==2
		// dorzucanie d³ugoœci przedzia³u do ud
		(*ud).add_row((B.x - A.x)());
#endif
	}
	return C;
}

solution lag(double a, double b, double epsilon, double gamma, int Nmax, matrix *ud, matrix *ad)
{
	// clear f_calls
	solution A(a), B(b), C, D, D_old(a);
	C.x = (a + b) / 2;
	A.fit_fun(ud, ad);
	B.fit_fun(ud, ad);
	C.fit_fun(ud, ad);
	double l, m;
	while (true)
	{
		l = A.y(0)*(pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0)*(pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0)*(pow(A.x(0), 2) - pow(B.x(0), 2));
		m = A.y(0)*(B.x(0) - C.x(0)) + B.y(0)*(C.x(0) - A.x(0)) + C.y(0)*(A.x(0) - B.x(0));
		if (m <= 0)
		{
			C.x = NAN;
			C.y = NAN;
			return C;
		}
		D.x = 0.5 * l/m;
		D.fit_fun(ud, ad);
		if (A.x <= D.x && D.x <= C.x)
		{
			if (D.y < C.y)
			{
				B = C;
				C = D;
			}
			else
				A = D;
		}
		else if (C.x <= D.x && D.x <= B.x)
		{
			if (D.y < C.y)
			{
				A = C;
				C = D;
			}
			else
				B = D;
		}
		else
		{
			C.x = NAN;
			C.y = NAN;
			return C;
		}
#if LAB_NO==1 && LAB_PART==2
		// dorzucanie d³ugoœci przedzia³u do ud
		(*ud).add_row((B.x - A.x)());
#endif
			if (B.x - A.x < epsilon || abs(D.x() - D_old.x()) < gamma || solution::f_calls > Nmax)
				return C;
		D_old = D;
	}
}
#endif
#if LAB_NO>1
solution HJ(matrix x0, double s, double alpha, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	solution XB(x0), XB_old, X;
	XB.fit_fun(ud, ad);
	while (true)
	{
		X = HJ_trial(XB, s, ud, ad);
		if (X.y < XB.y)
		{
			while (true)
			{
				XB_old = XB;
				XB = X;
#if LAB_NO==2 && LAB_PART==2
				(*ud).add_row(trans(XB.x));
#endif
					X.x = 2 * XB.x - XB_old.x;
					X.fit_fun(ud, ad);
				X = HJ_trial(X, s, ud, ad);
				if (X.y >= XB.y)
					break;
				if (solution::f_calls > Nmax)
					return XB;
			}
		}
		else
			s *= alpha;
		if (s < epsilon || solution::f_calls > Nmax)
			return XB;
	}
}

solution HJ_trial(solution XB, double s, matrix *ud, matrix *ad)
{
	int n = get_dim(XB);
	matrix D = ident_mat(n);
	solution X;
	for (int i = 0; i < n; ++i)
	{
		X.x = XB.x + s * D[i];
		X.fit_fun(ud, ad);
		if (X.y < XB.y)
			XB = X;
		else
		{
			X.x = XB.x - s * D[i];
			X.fit_fun(ud, ad);
			if (X.y < XB.y)
				XB = X;
		}
	}
	return XB;
}

solution Rosen(matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	solution X(x0), Xt;
	int n = get_dim(X);
	matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
	X.fit_fun(ud, ad);
	while (true)
	{
		for (int i = 0; i < n; ++i)
		{
			Xt.x = X.x + s(i) * D[i];
			Xt.fit_fun(ud, ad);
			if (Xt.y < X.y)
			{
				X = Xt;
				l(i) += s(i);
				s(i) *= alpha;
			}
			else
			{
				++p(i);
				s(i) *= -beta;
			}
		}
#if LAB_NO==2 && LAB_PART==2
		(*ud).add_row(trans(X.x));
#endif
			bool change = true;
		for (int i = 0; i < n; ++i)
			if (l(i) == 0 || p(i) == 0)
			{
				change = false;
				break;
			}
		if (change)
		{
			matrix Q(n, n), v(n, 1);
			for (int i = 0; i < n; ++i)
				for (int j = 0; j <= i; ++j)
					Q(i, j) = l(i);
			Q = D * Q;
			v = Q[0] / norm(Q[0]);
			D.set_col(v, 0);
			for (int i = 1; i < n; ++i) // i = j na wyk³adzie
			{
				matrix temp(n, 1);
				for (int j = 0; j < i; ++j) { // j = k na wyk³adzie
					temp = temp + trans(Q[i]) * D[j] * D[j];
				}
				v = (Q[i] - temp) / norm(Q[i] - temp);
				D.set_col(v, i);
			}
			s = s0;
			l = matrix(n, 1);
			p = matrix(n, 1);
		}
		double max_s = abs(s(0));
		for (int i = 1; i < n; ++i)
			if (max_s < abs(s(i)))
				max_s = abs(s(i));
		if (max_s < epsilon || solution::f_calls > Nmax)
			return X;
	}
}
#endif
#if LAB_NO>2
solution pen(matrix x0, double c0, double dc, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
	solution X(x0), X1;
	matrix c(2, new double[2]{ c0,dc });
	while (true)
	{
		X1 = sym_NM(X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud, &c);
		if (solution::f_calls > Nmax || norm(X1.x - X.x) < epsilon)
			return X1;
		X = X1;
		c(0) = c(0) * dc; // dc > 1 -> kara zewnêtrzna, < wewnêtrzna
	}
}

solution sym_NM(matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	matrix D = ident_mat(n);
	int N = n + 1; // liczba wierzcho³ków, wymiar problemu + 1
	solution *S = new solution[N];
	S[0].x = x0;
	S[0].fit_fun(ud, ad);
	for (int i = 1; i < N; ++i)
	{
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun(ud, ad);
	}
	solution PR, PE, PN;
	matrix pc;
	int i_min, i_max;
	while (true)
	{
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i)
		{
			if (S[i_min].y > S[i].y)
				i_min = i;
			if (S[i_max].y < S[i].y)
				i_max = i;
		}
		pc = matrix(n, 1);
		for (int i = 0; i < N; ++i)
			if (i != i_max)
				pc = pc + S[i].x;
				pc = pc / (N - 1);
				PR.x = pc + alpha * (pc - S[i_max].x);
			PR.fit_fun(ud, ad);
			if (S[i_min].y <= PR.y && PR.y < S[i_max].y)
				S[i_max] = PR;
		else if (PR.y < S[i_min].y)
		{
			PE.x = pc + gamma * (PR.x - pc);
			PE.fit_fun(ud, ad);
			if (PE.y < PR.y)
				S[i_max] = PE;
			else
				S[i_max] = PR;
		}
		else
		{
			PN.x = pc + beta * (S[i_max].x - pc);
			PN.fit_fun(ud, ad);
			if (PN.y < S[i_max].y)
				S[i_max] = PN;
			else
			{
				for (int i = 0; i < N; ++i)
					if (i != i_min)
					{
						S[i].x = delta * S[i].x + S[i_min].x;
						S[i].fit_fun(ud, ad);
					}
			}
		}
			double max_s = norm(S[0].x - S[i_min].x);
			for (int i = 1; i < N; ++i)
				if (max_s < norm(S[i].x - S[i_min].x))
					max_s = norm(S[i].x - S[i_min].x);
					if (max_s < epsilon || solution::f_calls > Nmax)
						return S[i_min];
	}
}
#endif
#if LAB_NO>3
solution SD(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n, 1), *P = new matrix[2]; // P[0] = x, P[1] = d
	solution h; // dlugosc kroku, h0 poczatkowa
	double *ab;
	while (true)
	{
		X.grad();
		d = -X.g; // kierunek = -gradient 
		if (h0 < 0) // zmienny krok
		{
			P[0] = X.x;
			P[1] = d;
			ab = expansion(0, 1, 1.2, Nmax, ud, P);
			h = golden(ab[0], ab[1], epsilon, Nmax, ud, P);
			X1.x = X.x + h.x * d;
		}
		else
			X1.x = X.x + h0 * d;
#if LAB_NO==4 && LAB_PART==2
		(*ud).add_row(trans(X1.x));
#endif
				if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
				{
					// dla h0 tylko tu liczymy fit_fun
					X1.fit_fun(ud, ad);
					return X1;
				}
				X = X1;
	}
}

solution CG(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n, 1), *P = new matrix[2];
	solution h;
	double *ab, beta;
	X.grad();
	d = -X.g;
		while (true)
		{
			if (h0<0)
			{
				P[0] = X.x;
				P[1] = d;
				ab = expansion(0, 1, 1.2, Nmax, ud, P);
				h = golden(ab[0], ab[1], epsilon, Nmax, ud, P);
				X1.x = X.x + h.x * d;
			}
			else
				X1.x = X.x + h0 * d;
#if LAB_NO==4 && LAB_PART==2
			(*ud).add_row(trans(X1.x));
#endif
				if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
				{
					X1.fit_fun(ud);
					return X1;
				}
			X1.grad();
			beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
			d = -X1.g + beta * d;
			X = X1;
		}
}

solution Newton(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n, 1), *P = new matrix[2];
	solution h;
	double *ab;
	while (true)
	{
		X.grad();
		X.hess();
		d = -inv(X.H) * X.g;
		if (h0 < 0) // zmienny krok
		{
			P[0] = X.x;
			P[1] = d;
			ab = expansion(0, 1, 1.2, Nmax, ud, P);
			h = golden(ab[0], ab[1], epsilon, Nmax, ud, P);
			X1.x = X.x + h.x * d;
		}
		else
			X1.x = X.x + h0 * d;
#if LAB_NO==4 && LAB_PART==2
		(*ud).add_row(trans(X1.x));
#endif
			if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax)
			{
				// dla h0 tylko tu liczymy fit_fun
				X1.fit_fun(ud, ad);
				return X1;
			}
		X = X1;
	}
}

solution golden(double a, double b, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	double alfa = (sqrt(5) - 1) / 2;
	solution A, B, C, D;
	A.x = a;
	B.x = b;
	C.x = B.x - alfa * (B.x - A.x);
	C.fit_fun(ud, ad);
	D.x = A.x + alfa *(B.x - A.x);
	D.fit_fun(ud, ad);
	while (true)
	{
		if (C.y < D.y)
		{
			B = D;
			D = C;
			C.x = B.x - alfa * (B.x - A.x);
				C.fit_fun(ud, ad);
		}
		else
		{
			A = C;
			C = D;
			D.x = A.x + alfa * (B.x - A.x);
				D.fit_fun(ud, ad);
		}
		if (B.x - A.x < epsilon || solution::f_calls > Nmax)
		{
			A.x = (A.x + B.x) / 2;
				A.fit_fun(ud, ad);
			return A;
		}
	}
}

#endif
#if LAB_NO>4
solution EA(int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	//mi liczebnosc populacji
	solution *P = new solution[mi + lambda];
	solution *Pm = new solution[mi];
	//random_device rd; - zbedne
	default_random_engine gen;
	gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
	normal_distribution<double> distr(0.0, 1.0);
	matrix IFF(mi, 1), temp(N, 2); //IFF macierz z przystosowaniami //temp - kopia osobnikia
	double r, s, s_IFF;
	double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5); //tau,tau1 - mutacja
	int j_min; //indeks najlepszego rozwiazania
	for (int i = 0; i < mi; ++i)
	{
		P[i].x = matrix(N, 2);
		for (int j = 0; j < N; ++j)
		{
			P[i].x(j, 0) = (limits(j, 1) - limits(j, 0))*rand_mat(1, 1)() + limits(j, 0);
			P[i].x(j, 1) = sigma0(j);
		}
		P[i].fit_fun(ud, ad);
		if (P[i].y < epsilon)
			return P[i];
	}
	while (true)
	{
		s_IFF = 0;
		//ruletka
		for (int i = 0; i<mi; ++i)
		{
			IFF(i) = 1 / P[i].y();
			s_IFF += IFF(i);
		}
		for (int i = 0; i<lambda; ++i)
		{
			r = s_IFF * rand_mat(1, 1)();
			s = 0;
			for (int j = 0; j<mi; ++j)
			{
				s += IFF(j);
				if (r <= s)
				{
					P[mi + i] = P[j]; // j - wylosowany osobnik
					break;
				}
			}
		}
		//mutajca
		for (int i = 0; i<lambda; ++i)
		{
			r = distr(gen);
			for (int j = 0; j<N; ++j)
			{
				P[mi + i].x(j, 1) *= exp(tau1*r + tau*distr(gen));
				P[mi + i].x(j, 0) += P[mi + i].x(j, 1)*distr(gen);
			}
		}
		//krzyzowanie
		for (int i = 0; i<lambda; i += 2)
		{
			r = rand_mat(1, 1)();
			temp = P[mi + i].x;  //jeden z rodzicow
			P[mi + i].x = r*P[mi + i].x + (1 - r)*P[mi + i + 1].x;  //pierwszy potomek
			P[mi + i + 1].x = r*P[mi + i + 1].x + (1 - r)*temp;  //drugi potomek
		}
		//ocena osobnikow
		for (int i = 0; i<lambda; ++i)
		{
			P[mi + i].fit_fun(ud, ad);
			if (P[mi + i].y < epsilon)  //ocena rozwiazania
				return P[mi + i];
		}
		//wskazanie najelpszych osobnikow
		for (int i = 0; i<mi; ++i)
		{
			j_min = 0;
			for (int j = 1; j<mi + lambda; ++j)
				if (P[j_min].y>P[j].y)
					j_min = j;
			Pm[i] = P[j_min];
			P[j_min].y = 1e10;
		}
		for (int i = 0; i < mi; ++i)
			P[i] = Pm[i];  //P[i] najlepsza populacja
		if (solution::f_calls>Nmax)
			return P[0];
	}
}
#endif
