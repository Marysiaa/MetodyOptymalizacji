#include"opt_alg.h"
#include<fstream>

// lab_2

double *expansion(double x0, double d, double alpha, int Nmax) {
	double *p = new double[2];
	solution X0(x0), X1(x0 + d);
	X0.fit_fun();
	X1.fit_fun();

	if (X0.y == X1.y) {
		p[0] = X0.x(0);
		p[1] = X1.x(0);
		return p;
	}
	else if (X0.y <= X1.y) {
		d *= -1;
		X1.x = X0.x + d;
		X1.fit_fun();
		if (X0.y < X1.y) {
			p[0] = X1.x(0);
			p[1] = X0.x(0) - d;
			return p;
		}
	}
	solution X2;
	int i = 1;
	while (true) {
		X2.x = x0 + pow(alpha, i)*d;
		X2.fit_fun();
		if (X2.y >= X1.y || solution::f_calls > Nmax) {
			break;
		}
		X0 = X1;
		X1 = X2;
		++i;
	}
	d > 0 ? p[0] = X0.x(0), p[1] = X2.x(0) : (p[0] = X2.x(0), p[1] = X0.x(0));
	return p;
}

solution fib(double a, double b, double epsilon) {
	int n = static_cast<int>(ceil(log2(sqrt(5)*(b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
	int *F = new int[n] {1, 1};
	for (int i = 2; i < n; ++i) {
		F[i] = F[i - 2] + F[i - 1];
	}
	solution A(a), B(b), C, D;
	C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
	D.x = A.x + B.x - C.x;
	C.fit_fun();
	D.fit_fun();

	//ofstream file("fibonacci_iterations.csv");
	for (int i = 0; i < n - 3; ++i) {
		//file << A.x << "\t" << B.x << endl;
		if (C.y < D.y) {
			B = D;
		}
		else {
			A = C;
		}
		C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
		D.x = A.x + B.x - C.x;
		C.fit_fun();
		D.fit_fun();
	}
	//file.close();
	return C;
}

solution lag(double a, double b, double epsilon, int Nmax, double gamma) {
	solution A(a), B(b), C, D;
	C.x = (a + b) / 2;
	A.fit_fun();
	B.fit_fun();
	C.fit_fun();
	double d0 = NAN;

	//ofstream file("lagrange_iterations.csv");
	while (true) {
		//file << A.x << "\t" << B.x << endl;
		matrix As(
			new double*[3]{
				new double[3]{ A.x(0)*A.x(0), A.x(0), 1 },
				new double[3]{ B.x(0)*B.x(0), B.x(0), 1 },
				new double[3]{ C.x(0)*C.x(0), C.x(0), 1 } }, 3, 3);
		matrix bs(new double[3]{ A.y(0), B.y(0), C.y(0) }, 3);
		matrix xs(3, 1);

		try {
			xs = inv(As)*bs;
		}
		catch (char*) {
			C = NAN;
			//file.close();
			return C;
		}
		if (xs(0) <= 0) {
			C = NAN;
			//file.close();
			return 0;
		}
		D.x = -xs(1) / (2 * xs(0));
		if (D.x < A.x || D.x > B.x) {
			C = NAN;
			//file.close();
			return C;
		}
		D.fit_fun();
		if (A.x < D.x && D.x < C.x)
		{
			if (D.y < C.y)
			{
				solution temp(C);
				C = D;
				B = temp;
			}
			else
			{
				A = D;
			}
		}
		else if (C.x < D.x && D.x < B.x)
		{
			if (D.y < C.y)
			{
				A = C;
				C = D;
			}
			else
			{
				B = D;
			}
		}
		if ((B.x - A.x) < epsilon || solution::f_calls > Nmax) {
			//file.close();
			return D;
		}
		// test stacjonarnosci
		if (d0 != NAN && (abs(D.x(0) - d0) < gamma)) {
			//file.close();
			return D;
		}
		d0 = D.x(0);
	}
}

// lab_3

solution HJ(matrix x0, double s, double alpha, double epsilon, int Nmax) {
	solution XB, XB_old, X;
	XB.x = x0;
	XB.fit_fun();
	while (true) {
		X = HJ_trial(XB, s);
		// etap roboczy przynosi poprawe
		if (X.y < XB.y) {
			while (true) {
				XB_old = XB;
				XB = X;
				X.x = 2.0 * XB.x - XB_old.x;
				X.fit_fun();
				X = HJ_trial(X, s);
				if (X.y >= XB.y) {
					break;
				}
				if (solution::f_calls > Nmax) {
					return XB;
				}
			}
		}
		// etap roboczy nie przynosi poprawy
		else {
			s *= alpha;
			if (s < epsilon || solution::f_calls > Nmax) {
				return XB;
			}
		}
	}
}

solution HJ_trial(solution XB, double s) {
	int *n = get_size(XB.x);
	matrix D = unit_mat(n[0]);
	solution X;
	for (int i = 0; i < n[0]; i++) {
		X.x = XB.x + s * D[i];
		X.fit_fun();
		if (X.y < XB.y) {
			XB = X;
		}
		else {
			X.x = XB.x - s * D[i];
			X.fit_fun();
			if (X.y < XB.y) {
				XB = X;
			}
		}
	}
	delete[] n;
	return XB;
}

solution Rosenbrock(matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax) {
	int *n = get_size(x0);
	matrix l(n[0], 1), p(n[0], 1), s(s0), D = unit_mat(n[0]);
	solution X, Xt;
	X.x = x0;
	X.fit_fun();
	while (true) {
		// sprawdzamy wszystkie kierunki
		for (int i = 0; i < n[0]; i++) {
			Xt.x = X.x + s(i) * D[i];
			Xt.fit_fun();
			if (Xt.y < X.y) {
				X = Xt;
				l(i) += s(i); // przesuniecie wzgledne;
				s(i) *= alpha;
			}
			else {
				++p(i);
				s(i) *= -beta;
			}
		}
		// sprawdzenie czy mamy wykonac obrot bazy kierunkow
		bool change = true;
		for (int i = 0; i < n[0]; i++) {
			if (l(i) == 0 || p(i) == 0) {
				change = false;
				break;
			}
		}
		if (change) {
			matrix Q(n[0], n[0]), v(n[0], 1);
			for (int i = 0; i < n[0]; i++) {
				for (int j = 0; j <= i; j++) {
					Q(i, j) = l(i);
				}
			}
			Q = D * Q;
			v = Q[0] / norm(Q[0]);
			D = set_col(D, v, 0);
			for (int i = 1; i < n[0]; i++) {
				matrix temp(n[0], 1);
				for (int j = 0; j < i; j++) {
					temp = temp + trans(Q[i]) * D[j] * D[j];
				}
				v = (Q[i] - temp) / norm(Q[i] - temp);
				D = set_col(D, v, i);
			}
			s = s0;
			l = matrix(n[0], 1);
			p = matrix(n[0], 1);
		}
		// sprawdzenie warunkow stopu
		double maxS = abs(s(0));
		for (int i = 1; i < n[0]; i++) {
			if (maxS < abs(s(i))) {
				maxS = abs(s(i));
			}
		}
		if (maxS < epsilon || solution::f_calls > Nmax) {
			delete[] n;
			return X;
		}
	}
}

// lab_4

solution sym_NM(matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix A, bool wewnetrzna) {
	int *n = get_size(x0);
	matrix D = unit_mat(n[0]);
	//liczba wierzcholkow
	int N = n[0] + 1;
	solution *S = new solution[N];
	S[0].x = x0;
	S[0].fit_fun(A, wewnetrzna);

	//generowanie reszty wierzvholkow
	for (int i = 1; i < N; ++i) {
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun(A, wewnetrzna);
	}

	solution p_o, p_e, p_z;
	matrix p_sr;
	int i_min, i_max;

	while (true) {
		//wyznaczenie max min i srodka ciezkowsci
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i) {
			if (S[i_min].y > S[i].y) {
				i_min = i;
			}
			if (S[i_max].y < S[i].y) {
				i_max = i;
			}
		}
		p_sr = matrix(n[0], 1);
		for (int i = 0; i < N; ++i) {
			if (i != i_max) {
				p_sr = p_sr + S[i].x;
			}
		}
		p_sr = p_sr / (N - 1);

		//odbicie
		p_o.x = p_sr + alpha * (p_sr - S[i_max].x);
		p_o.fit_fun(A, wewnetrzna);
		if (S[i_min].y <= p_o.y  && p_o.y < S[i_max].y) {
			S[i_max] = p_o;
		}
		else if (p_o.y < S[i_min].y)
		{
			//ekspansja
			p_e.x = p_sr + gamma * (p_o.x - p_sr);
			p_e.fit_fun(A, wewnetrzna);

			if (p_e.y < p_o.y) {
				S[i_max] = p_e;
			}
			else {
				S[i_max] = p_o;
			}
		}
		else {
			//zawezenie
			p_z.x = p_sr + beta * (S[i_max].x - p_sr);
			p_z.fit_fun(A, wewnetrzna);

			if (p_z.y < S[i_max].y) {
				S[i_max] = p_z;
			}
			else
			{
				//redukcja
				for (int i = 0; i < N; ++i) {
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x + S[i_min].x);
						S[i].fit_fun(A, wewnetrzna);
					}
				}
			}
		}

		//kryteria stopu
		double max_s = norm(S[0].x - S[i_min].x);
		for (int i = 1; i < N; i++) {
			if (max_s < norm(S[i].x - S[i_min].x)) {
				max_s = norm(S[i].x - S[i_min].x);
			}
		}
		if (max_s<epsilon || solution::f_calls>Nmax) {
			return S[i_min];
		}
	}
}

solution penalty(matrix x0, double c, double a, double epsilon, int Nmax, bool wewnetrzna) {
	// wspolczynniki metody sympleks
	double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
	matrix A(new double[2]{ c, a }, 2);

	solution X, X1;
	X.x = x0;
	while (true) {
		X1 = sym_NM(X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, A, wewnetrzna);
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax) {
			return X1;
		}

		if (!wewnetrzna) {
			A(0) *= 2;		// kara zewnetrzna 
		}
		else {
			A(0) *= 0.1;	// kara wewnetrzna
		}

		X = X1;
	}
}

// lab_5

solution golden(double a, double b, double epsilon, int Nmax, matrix P) {
	double alfa = (sqrt(5.0) - 1) / 2;
	solution A, B, C, D;
	A.x = a;
	B.x = b;
	C.x = B.x - alfa * (B.x - A.x);
	C.fit_fun(P);
	D.x = A.x + alfa * (B.x - A.x);
	D.fit_fun(P);

	while (true) {
		if (C.y < D.y) {
			B = D;
			D = C;
			C.x = B.x - alfa * (B.x - A.x);
			C.fit_fun(P);
		}
		else {
			A = C;
			C = D;
			D.x = A.x + alfa * (B.x - A.x);
			D.fit_fun(P);
		}

		// kryteria stopu 
		if (B.x - A.x < epsilon || solution::f_calls>Nmax) {
			A.x = (A.x + B.x) / 2.0;
			A.fit_fun(P);
			return A;
		}
	}
}

//GRADIENTOWE
//metoda najwiekszego spadku
solution SD(matrix x0, double epsilon, int Nmax, double h0) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	double b;

	while (true) {
		//liczymy gradient z X 
		X.grad();
		d = -X.g;

		if (h0 == 0) {
			// wstawiamy do macierzy punkt stardowy i wektor kierunku
			P = set_col(P, X.x, 0);
			//w przyszlosci trzeba dodac wyliczenie dla dolnego ograniczenia "a"
			P = set_col(P, d, 1);
			// liczone na podstawie X i kierunku 
			// b gorna granica kroku
			b = compute_b(X.x, d);
			h = golden(0, b, epsilon, Nmax, P);
		}
		else {
			h = h0;
		}
		// przesuwamy sie o h.x*d
		X1.x = X.x + h.x*d;
		//kryteria stopu 
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
			X1.fit_fun();
			return X1;
		}
		// jesli warunek stopu nie dziala 
		X = X1;
	}
}

solution Newton(matrix x0, double epsilon, int Nmax, double h0) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	double b;

	while (true) {
		X.grad();
		X.hess();
		d = -inv(X.H)*X.g;

		if (h0 == 0) {
			P = set_col(P, X.x, 0);  // wstawiamy do macierzy punkt stardowy i wektor kierunku
			P = set_col(P, d, 1);
			//w przyszlosci trzeba dodac wyliczenie dla dolnego ograniczenia "a"
			b = compute_b(X.x, d);  // liczone na podstawie X i kierunku 
			h = golden(0, b, epsilon, Nmax, P);  // b gorna granica kroku
		}
		else {
			h = h0;
		}

		X1.x = X.x + h.x*d; // przesuwamy sie o h.x*d
							//kryteria stopu 
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
			X1.fit_fun();
			return X1;
		}
		X = X1; // jesli warunek stopu nie dziala 
	}
}

//greadienty sprzezone
solution CG(matrix x0, double epsilon, int Nmax, double h0) {
	int *n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2);
	solution h;
	double b, beta;

	X.grad(); //liczymy gradient z X 
	d = -X.g;
	while (true) {

		if (h0 == 0) {
			P = set_col(P, X.x, 0);  // wstawiamy do macierzy punkt stardowy i wektor kierunku
			P = set_col(P, d, 1);
			//w przyszlosci trzeba dodac wyliczenie dla dolnego ograniczenia "a"
			b = compute_b(X.x, d);  // liczone na podstawie X i kierunku 
			h = golden(0, b, epsilon, Nmax, P);  // b gorna granica kroku
		}
		else {
			h = h0;
		}

		X1.x = X.x + h.x*d; // przesuwamy sie o h.x*d
							//kryteria stopu 
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
			X1.fit_fun();
			return X1;
		}
		X1.grad();

		beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
		d = -X1.g + beta * d;
		X = X1; // jesli warunek stopu nie dziala 
	}
}

// a = 0 - nie mozemy pojsc w druga strone niz wskazuje kierunek d
// obliczenie b - drugiego punktu do zawezania przedzialu, b > 0
//wyliczenie najmniejszego mozliwego b w kazdym kierunku
double compute_b(matrix x, matrix d) {
	int *n = get_size(x);
	double b = 1e9, bi; // maksymalne ogranicznie w itym kierunku 

	for (int i = 0; i < n[0]; ++i) {
		if (d(1) == 0)
			bi = 1e9;
		else if (d(i) > 0)
			bi = (10 - x(i)) / d(i);  // skıd 10 ? gorne ograniczenie  (mozna byloby je przekazywaĉ przez argumenty) xd
		else
			bi = (-10 - x(i)) / d(i); // -10 dolne ogranicznie 
									  //bi zawsze jest dodatnie po tej operacji
		if (b > bi)
			b = bi;
	}
	return b;
}

// lab_6

solution Powell(matrix x0, double epsilon, int Nmax, matrix w) {
	//w - wspolczynniki wagowe
	int *n = get_size(x0);
	matrix D = unit_mat(n[0]);
	//A - zawira punkt startowy P.x, kierunek D[i] oraz wektor [w 0 0 ..]
	matrix A(n[0], 3);
	A(0, 2) = w(0);
	double * ab;
	solution X, P, h;
	X.x = x0;
	while (true) {
		P = X;
		for (int i = 0; i < n[0]; ++i) {
			A = set_col(A, P.x, 0);
			A = set_col(A, D[i], 1);
			ab = compute_ab(P.x, D[i]);
			h = golden(ab[0], ab[1], epsilon, Nmax, A);
			P.x = P.x + h.x*D[i];
		}
		if (norm(P.x - X.x) < epsilon || solution::f_calls > Nmax) {
			P.fit_fun();
			return P;
		}
		for (int i = 0; i < n[0] - 1; ++i) {
			D = set_col(D, D[i + 1], i);
		}
		D = set_col(D, P.x - X.x, n[0] - 1);
		A = set_col(A, P.x, 0);
		A = set_col(A, D[n[0] - 1], 1);
		ab = compute_ab(P.x, D[n[0] - 1]);
		h = golden(ab[0], ab[1], epsilon, Nmax, A);
		X.x = P.x + h.x*D[n[0] - 1];
	}
}

// funkcja liczy ograniczenia a i b 
double * compute_ab(matrix x, matrix d) {
	int * n = get_size(x);
	double * ab = new double[2]{ -1e9, 1e9 };
	double ai, bi;
	// dolne ograniczenia
	double lb[2] = { 0.2, 0.01 };
	// gorne ograniczenia
	double ub[2] = { 1, 0.05 };
	for (int i = 0; i < n[0]; i++) {
		//brak ograniczen
		if (d(i) == 0) {
			ai = -1e9;
			bi = 1e9;
		}
		//kirunek jest dodatni
		else if (d(i) > 0) {
			ai = (lb[i] - x(i)) / d(i);
			bi = (ub[i] - x(i)) / d(i);
		}
		//kierunek jest ujemny
		else {
			ai = (ub[i] - x(i)) / d(i);
			bi = (lb[i] - x(i)) / d(i);
		}
		if (ab[0] < ai) {
			ab[0] = ai;
		}
		if (ab[1] > bi) {
			ab[1] = bi;
		}
	}
	return ab;
}

// lab_7

//lb - dolne ograniczenie na wszystkie zmienne
//ub - gorne ograniczenie na wszystkie zmienne
solution EA(int N, matrix lb, matrix ub, double epsilon, int Nmax, double sigma = 1) {
	int mi = 20, lambda = 40;
	// populacja dlugosci mi + lambda
	solution *P = new solution[mi + lambda];
	// kopia populacji dlugosci mi
	solution*Pm = new solution[mi];
	random_device rd;
	default_random_engine gen;
	gen.seed(chrono::system_clock::now().time_since_epoch().count());
	// rozklad normalny
	normal_distribution<double> distr(0.0, 1.0);
	matrix IFF(mi, 1), temp(N, 2);
	double r, s, s_IFF;
	// tau = alpha, tau1 = beta
	// 2N - liczba zmiennych decyzyjnych (problem dwuwymiarowy)
	double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5);
	int j_min;

	// generowanie populacji poczatkowej
	for (int i = 0; i < mi; ++i) {
		P[i].x = matrix(N, 2);
		for (int j = 0; j < N; ++j) {
			P[i].x(j, 0) = (ub(j) - lb(j))*rd() / rd.max() + lb(j);
			P[i].x(j, 1) = sigma;
		}
		P[i].fit_fun();
		if (P[i].y < epsilon) {
			return P[i];
		}
	}
	while (true) {
		// selekcja - zastosowana regula kola ruletki
		// s_IFF - suma funkcji przystosowania calej populacji
		s_IFF = 0;
		for (int i = 0; i < mi; ++i) {
			// IFF - wartosc przystosowania osobnika
			IFF(i) = 1.0 / P[i].y(0);
			s_IFF += IFF(i);
		}
		for (int i = 0; i < lambda; ++i) {
			// losujemy liczbe z przedzialu od 0 do s_IFF 
			// prawdopodobienstwo wylosowania osobnika = IFF(i) / s_IFF
			r = s_IFF * rd() / rd.max();
			s = 0;
			for (int j = 0; j < mi; j++) {
				s += IFF(j);
				if (s >= r) {
					// kopiujemy j-ty osobnik do populacji tymczasowej (znajdujacej sie 'pod' wlasciwa populacja)
					P[i + mi] = P[j];
					break;
				}
			}
		}
		// na kazdym z lambda osobnikow wylosowanych w etapie selekcji przeprowadzamy mutacje i krzyzowanie
		// mutacja
		for (int i = 0; i < lambda; ++i) {
			r = distr(gen);
			for (int j = 0; j < N; ++j) {
				// modyfikacja zakresu mutacji (zmutowanie sigmy)
				P[mi + i].x(j, 1) *= exp(tau1*r + tau * distr(gen));
				// mutacja x-ow
				P[mi + i].x(j, 0) += P[mi + i].x(j, 1)*distr(gen);
			}
		}
		// krzyzowanie (zakladamy ze liczba osobnikow w populacji mi jest parzysta)
		for (int i = 0; i < lambda; i += 2) {
			// krzyzowanie usredniajace
			r = 1.0* rd() / rd.max();
			// kopia osobnika (ktorego zmodyfikujemy jako pierwszego)
			temp = P[mi + i].x;
			// modyfikacja osobnikow
			P[mi + i].x = r * P[mi + i].x + (1 - r)*P[mi + i + 1].x;
			P[mi + i + 1].x = r * P[mi + i + 1].x + (1 - r)*temp;
		}
		for (int i = 0; i < lambda; ++i) {
			// ocena osobnikow
			P[mi + i].fit_fun();
			if (P[mi + i].y < epsilon) {
				return P[mi + i];
			}
		}
		// wybieramy najlepszych osobnikow
		// na pewno rozwiazanie najlepsze nie bedzie utracone
		for (int i = 0; i < mi; ++i) {
			// najlepszy osobnik
			j_min = 0;
			for (int j = 1; j < mi + lambda; ++j) {
				if (P[j_min].y > P[j].y)
					j_min = j;
			}
			Pm[i] = P[j_min];
			P[j_min].y = 1e9;
		}
		// przepisanie pomulacji Pm do P
		for (int i = 0; i < mi; ++i) {
			P[i] = Pm[i];
		}
		if (solution::f_calls > Nmax) {
			return P[0];
		}
	}
}