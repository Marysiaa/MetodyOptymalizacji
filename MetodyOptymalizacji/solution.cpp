//Do not edit the code below (unless you know what you are doing)

#include"solution.h"

int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

solution::solution(double L)
{
	x = matrix(L);
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(double *A, int n)
{
	x = matrix(A, n);
	g = NAN;
	H = NAN;
	y = NAN;
}

//You can edit the following code

void solution::fit_fun(matrix A, bool wewnetrzna) // wzor funkcji celu
{
	// przyklad
	//y = x + inv(x*x);

	// zadanie
	// lab_2

	/*
	matrix Y0(new double[3]{ 5,1,10 }, 3);
	matrix *Y = solve_ode(0, 1, 1000, Y0, x);
	int *w = get_size(Y[1]);
	double max = Y[1](0, 2);
	for (int i = 1; i < w[0]; ++i) {
		if (max < Y[1](i, 2)) {
		max = Y[1](i, 2);
		}
	}
	y = abs(max - 50);
	++f_calls;
	*/

	// lab_3
	/*
	double a_ref = 3.14, o_ref = 0;
	matrix Y0(2, 1);
	matrix *Y = solve_ode(0, 0.1, 100, Y0, x);

	y(0) = 0;
	int *n = get_size(Y[1]);

	for (int i = 0; i < n[0]; i++) {
		y = y + 10 * pow(a_ref - Y[1](i, 0), 2) +
			pow(o_ref - Y[1](i, 1), 2) +
			pow(x(0)*(a_ref - Y[1](i, 0)) + x(1)*(o_ref - Y[1](i, 1)), 2);
	}
	y = y * 0.1; // dt = 0.1
	++f_calls;
	delete[] n;
	delete[] Y;
	*/

	// lab_4
	/*
	double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
	y = sin(arg) / arg;

	if (!wewnetrzna) {
		// kara zewnetrzna
		if (-x(0) + 1 > 0) {
			y = y + A(0)*pow(-x(0) + 1, 2);
		}
		if (-x(1) + 1 > 0) {
			y = y + A(0)*pow(-x(1) + 1, 2);
		}
		if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1) > 0) {
			y = y + A(0) * pow(sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1), 2);
		}
	}
	else {
		// kara wewnetrzna
		if (-x(0) + 1 > 0) {
			y = 1e10;
		}
		else {
			y = y - A(0) / (-x(0) + 1);
		}
		if (-x(1) + 1 > 0) {
			y = 1e10;
		}
		else {
			y = y - A(0) / (-x(1) + 1);
		}
		if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1) > 0) {
			y = 1e10;
		}
		else {
			y = y - A(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - A(1));
		}
	}

	++f_calls;
	*/

	// lab_5
	/*
	int *n = get_size(A);
	if (n[1] == 1) {
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
		++f_calls;
	}
	else {
		solution T;
		T.x = A[0] + x*A[1];
		T.fit_fun();
		y = T.y;
	}
	*/

	// lab_6
	/*
	int *n = get_size(A);
	double ro = 7800, P = 1e3, E = 207e9;
	if (n[1] == 1) {
		y = matrix(3, 1);
		//masa walca
		y(0) = ro*x(0)*3.14*pow(x(1), 2) / 4;
		//ugiecie
		y(1) = 64 * P * pow(x(0), 3) / (3 * E * 3.14 * (pow(x(1), 4)));
		//sigma - aktualne naprezenie
		y(2) = 32 * P*x(0) / (3.14*pow(x(1), 3));
		++f_calls;
	}
	else {
		//wywolanie z metody zlotego podzialu
		solution T;
		T.x = A[0] + x*A[1];
		T.fit_fun();

		//NORMALIZACJA MASY I UGIECIA DO PRZEDZIALU 0-1
		//T.y(0) => y0  i  T.y(1) => y1
		double y0 = (T.y(0) - 0.122) / (15.315 - 0.122);
		double y1 = T.y(1) / 0.005;

		y = A(0, 2)*y0 + (1 - A(0, 2))*y1;

		// funkcja kary przy naruszeniu ograniczenia
		if (T.y(1) > 0.005) {
			// naruszone ugiecie - ugiecie wieksze niz dopuszczalne
			y = y + 1e6*pow(T.y(1) - 0.005, 2);
		}
		if (T.y(2) > 300e6) {
			// naruszone naprezenie - naprezenie wieksze niz dopuszczalne
			y = y + 1e6*pow(T.y(2) - 300e6, 2);
		}
	}
	*/

	// lab_7

	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5*3.14*x(0)) - cos(2.5*3.14*x(1)) + 2;
	++f_calls;

}

void solution::grad()
{
	//g = NAN;

	// lab_5
	g = matrix(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
	++g_calls;
}

void solution::hess()
{
	//H = NAN;

	// lab_5
	H = matrix(2, 2);
	H(0, 0) = 10;
	H(1, 0) = H(0, 1) = 8;
	H(1, 1) = 10;
	++H_calls;
}

