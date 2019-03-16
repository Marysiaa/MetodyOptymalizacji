//Do not edit this file (unless you know what you are doing)

#ifndef SOLUTION_H
#define SOLUTION_H

#include"matrix.h"
#include"ode_solver.h"

class solution
{
public:
	matrix x;
	matrix g; // gradient
	matrix H; // hesjan
	matrix y;
	static int f_calls;
	static int g_calls;
	static int H_calls;
	solution(double = NAN);
	solution(double *, int);
	void fit_fun(matrix = 0.0, bool = false);
	void grad();
	void hess();
};

#endif