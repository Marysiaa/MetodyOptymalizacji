#include<random>
#include<fstream>
#include"opt_alg.h"
#include"ode_solver.h"

int main()
{
	try
	{
		cout << "Hello" << endl << endl;

		matrix Y0(new double[3]{ 5, 1, 10 }, 3); // warunki poczatkowe: 5 m^3 w zbiorniku A, 1 m^3 w zbiorniku B, 10 oC w zbiorniku B
		matrix C(1); // parametry, ktore bedziemy optymalizowac - na razie nie istotne
		
		matrix *S = solve_ode(0, 1, 1000, Y0, C);

		/*
		ofstream f1("time.csv");
		f1 << S[0] << endl;
		ofstream f2("solution.csv");
		f2 << S[1] << endl;
		f1.close();
		f2.close();
		*/
	}
	catch (char * EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
