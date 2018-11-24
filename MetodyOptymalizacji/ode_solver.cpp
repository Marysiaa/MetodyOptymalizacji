//Do not edit the code below (unless you know what you are doing)

#include"ode_solver.h"

matrix *solve_ode(double t0, double dt, double tend, matrix &Y0, matrix &C)
{
	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	if (N < 2)
		throw "The time interval is defined incorrectly";
	int *s = get_size(Y0);
	if (s[1] != 1)
		throw "Initial condition must be a vector";
	int n = s[0];
	delete[]s;
	matrix *S = new matrix[2]{ matrix(N), matrix(n,N) };
	S[0](0) = t0;
	for (int i = 0; i < n; ++i)
		S[1](i, 0) = Y0(i);
	matrix k1(n), k2(n), k3(n), k4(n);
	for (int i = 1; i < N; ++i)
	{
		S[0](i) = S[0](i - 1) + dt;
		k1 = dt*diff(S[0](i - 1), S[1][i - 1], C);
		k2 = dt*diff(S[0](i - 1) + 0.5*dt, S[1][i - 1] + 0.5*k1, C);
		k3 = dt*diff(S[0](i - 1) + 0.5*dt, S[1][i - 1] + 0.5*k2, C);
		k4 = dt*diff(S[0](i - 1) + dt, S[1][i - 1] + k3, C);
		for (int j = 0; j < n; ++j)
			S[1](j, i) = S[1](j, i - 1) + (k1(j) + 2 * k2(j) + 2 * k3(j) + k4(j)) / 6;
	}
	S[1] = trans(S[1]);
	return S;
	// S[O] - wektor zawierajay kolejne chwile czasowe
	// S[1] - macierz zawierajaca rozwiazania poszczegolnych rownan rozniczkowych (kolumny) w danych chwilach czasowych (wiersze)
}

//You can edit the following code

matrix diff(double t, matrix &Y, matrix &C)
{
	// rownanie rozniczowe
	// Y - aktualny stan
	// C - parametr który optymalizujemy
	matrix dY(Y);

	double PA = 1, PB = 1, a = 0.98, b = 0.63, g = 9.81, Tin = 10, Fin = 0.01, DA = 0.005, DB = 0.00365665, TA = 90;

	dY(0) = -(Y(0) > 0 ? a*b*DA*pow(2 * g*Y(0) / PA, 0.5) : 0);
	dY(1) = Fin - dY(0) - (Y(1) > 0 ? a*b*DB*pow(2 * g*Y(1) / PB, 0.5) : 0);
	dY(2) = Fin / Y(1) * (Tin - Y(2)) - dY(0) / Y(1) * (TA - Y(2));

	return dY;
}