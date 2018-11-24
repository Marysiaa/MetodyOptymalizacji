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

void solution::fit_fun() // wzor funkcji celu
{
	y = NAN; 
	++f_calls;
}

void solution::grad()
{
	g = NAN;
	++g_calls;
}

void solution::hess()
{
	H = NAN;
	++H_calls;
}

