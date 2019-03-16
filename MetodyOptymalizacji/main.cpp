#include<random>
#include<fstream>
#include"opt_alg.h"
#include"ode_solver.h"

void calculationForSigma(double sigma, ofstream & file) {
	double epsilon = 1e-3;
	int N = 2, Nmax = 50000;
	matrix lb(N, 1), ub(N, 1);
	lb(0) = lb(1) = -1;
	ub(0) = ub(1) = 1;

	for (int i = 0; i < 100; i++) {
		solution::f_calls = 0;
		solution sol = EA(N, lb, ub, epsilon, Nmax, sigma);
		file << sigma << "," << sol.x(0) << "," << sol.x(1) << "," << sol.y(0) << "," << solution::f_calls << "," << ((sol.y[0] < epsilon) ? "tak" : "nie") << endl;
	}
}

int main()
{
	//algorytmy ewolucyjne i genetyczne 

	try
	{
		//cout << "LAB_1" << endl << endl;
		/*
		//matrix Y0(new double[3]{ 5, 1, 10 }, 3); // warunki poczatkowe: 5 m^3 w zbiorniku A, 1 m^3 w zbiorniku B, 10 oC w zbiorniku B
		//matrix C(1); // parametry, ktore bedziemy optymalizowac - na razie nie istotne
		//
		//matrix *S = solve_ode(0, 1, 1000, Y0, C);

		//ofstream f1("time.csv");
		//f1 << S[0] << endl;
		//ofstream f2("solution.csv");
		//f2 << S[1] << endl;
		//f1.close();
		//f2.close();
		*/

		//cout << "LAB_2" << endl << endl;
		/*
		double x0,
			a = 1e-4,
			b = 1e-2,
			d = 1e-5,
			alpha = 1.2,
			epsilon = 1e-6,
			gamma = 1e-7;
		int Nmax = 1000;

		random_device rd;
		solution opt;
		double * p;
		*/

		/*
		x0 = (b - a)*rd() / rd.max() + a;
		cout.precision(10);
		cout << x0 << endl << endl;
		p = expansion(x0, d, alpha, Nmax);
		cout << p[0] << '\t' << p[1] << endl;
		cout << solution::f_calls << endl << endl;

		solution::f_calls = 0;
		opt = fib(p[0], p[1], epsilon);
		cout << "fib" << endl << opt.x << endl << opt.y << endl << solution::f_calls << endl << endl;

		solution::f_calls = 0;
		opt = lag(p[0], p[1], epsilon, Nmax, gamma);
		cout << "lag" << endl << opt.x << endl << opt.y << endl << solution::f_calls << endl;
		*/

		/*
		ofstream fExp("solution_Exp.csv");
		ofstream fFib("solution_Fib.csv");
		ofstream fLag("solution_Lag.csv");

		alpha = 2.2; //1.2; //1.8;
		for (int i = 0; i < 100; i++) {
			x0 = (b - a)*rd() / rd.max() + a;

			solution::f_calls = 0;
			p = expansion(x0, d, alpha, Nmax);
			fExp << x0 << "\t" << solution::f_calls << "\t" << p[0] << "\t" << p[1]<<endl;

			solution::f_calls = 0;
			opt = fib(p[0], p[1], epsilon);
			fFib << opt.x << "\t" << opt.y << "\t" << solution::f_calls << endl;

			solution::f_calls = 0;
			opt = lag(p[0], p[1], epsilon, Nmax, gamma);
			fLag << opt.x << "\t" << opt.y << "\t" << solution::f_calls << endl;

			cout << i << endl;
		}
		fExp.close();
		fFib.close();
		fLag.close();
		*/

		/*
		// bez zawezenia przedzialu
		solution::f_calls = 0;
		opt = fib(1, 100, epsilon);
		cout << opt.x << "\t" << opt.y << "\t" << solution::f_calls << endl;

		solution::f_calls = 0;
		opt = lag(1, 100, epsilon, Nmax, gamma);
		cout << opt.x << "\t" << opt.y << "\t" << solution::f_calls << endl;
		*/

		//cout << "LAB_3" << endl << endl;
		/*
		double alpha, beta, epsilon = 1e-3, s;
		int Nmax = 1000;
		matrix x0(2, 1);
		random_device rd;
		x0(0) = 10.0*rd() / rd.max();
		x0(1) = 10.0*rd() / rd.max();
		cout << "x0: " << endl << x0 << endl << endl;

		s = 0.5;
		alpha = 0.5;
		solution opt_hj = HJ(x0, s, alpha, epsilon, Nmax);
		cout << "opt_hf = " << opt_hj.x << "\t" << opt_hj.y << endl;
		cout << "f_calls = " << solution::f_calls << endl << endl;

		matrix s0(new double[2]{ s, s }, 2);
		alpha = 2;
		beta = 0.5;
		solution::f_calls = 0;
		solution opt_r = Rosenbrock(x0, s0, alpha, beta, epsilon, Nmax);
		cout << "opt_r = " << opt_r.x << "\t" << opt_r.y << endl;
		cout << "f_calls = " << solution::f_calls << endl << endl;
		*/

		/*
		double alpha_hj = 0.5, alpha_ros = 2.0, beta = 0.5, epsilon = 1e-3;
		int Nmax = 1000;
		matrix x0(2, 1);
		random_device rd;

		ofstream fStartPoints05("solution_start_points4.csv");
		ofstream fHJ05("solution_hj4.csv");
		ofstream fRos05("solution_ros4.csv");

		double s = 4;
		matrix s0(new double[2]{ s, s }, 2);

		solution opt;
		for (int i = 0; i < 100; i++) {

			x0(0) = 10.0*rd() / rd.max();
			x0(1) = 10.0*rd() / rd.max();
			fStartPoints05 << x0(0) << "," << x0(1) << endl;

			solution::f_calls = 0;
			opt = HJ(x0, s, alpha_hj, epsilon, Nmax);
			fHJ05 << opt.x << opt.y(0) << "," << solution::f_calls << endl;

			solution::f_calls = 0;
			opt = Rosenbrock(x0, s0, alpha_ros, beta, epsilon, Nmax);
			fRos05 << opt.x << opt.y(0) << "," << solution::f_calls << endl;

			cout << i << endl;
		}
		fStartPoints05.close();
		fHJ05.close();
		fRos05.close();

		ofstream fStartPoints("solution_start_points.csv");
		ofstream fHJ("solution_hj.csv");
		ofstream fRos("solution_ros.csv");
		*/

		//cout << "LAB_4" << endl << endl;
		/*
		double epsilon = 1e-4, c, a;
		int Nmax = 10000;
		matrix x0(2, 1);
		random_device rd;


		//do {
		//	x0(0) = 4.0 *rd() / rd.max() + 1;
		//	x0(1) = 4.0 *rd() / rd.max() + 1;
		//	// punkt startowy musi znajdowac sie w obszarze rozwiazan dopuszczalnych
		//} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) - a > 0);

		//solution opt = penalty(x0, c, a, epsilon, Nmax);
		//cout << "start point:\n" << x0 << endl << endl;
		//cout << "opt:\nx:\n" << opt.x << endl << "y:\n" << opt.y << endl << endl;
		//cout << "r:\t" << sqrt(pow(opt.x(0), 2) + pow(opt.x(1), 2)) << endl;
		//cout << "f_calls:\t" << solution::f_calls << endl;


		// ---------------------------------------------------------------------------------
		// a=4; - nagly wzrost wywolan funckji celu gdy jestesmy poza obszarem dopuszczalnym
		// a = 4, 4.4934, 5

		a = 5;
		// kara zewnetrzna
		double c_zew = 1;
		// kara wewnetrzna
		double c_wew = 1000;

		ofstream fStartPoints("solution_start_points5.csv");
		ofstream zewnetrzna("solution_zewnetrzna5.csv");
		ofstream wewnetrzna("solution_wewnetrzna5.csv");

		solution opt;
		for (int i = 0; i < 100; i++) {

			do {
				x0(0) = 4.0 *rd() / rd.max() + 1;
				x0(1) = 4.0 *rd() / rd.max() + 1;
			} while (sqrt(pow(x0(0), 2) + pow(x0(1), 2)) - a > 0);
			fStartPoints << x0(0) << "," << x0(1) << endl;

			solution::f_calls = 0;
			opt = penalty(x0, c_zew, a, epsilon, Nmax, false);
			zewnetrzna << opt.x << sqrt(pow(opt.x(0), 2) + pow(opt.x(1), 2)) << "," << opt.y(0) << "," << solution::f_calls << endl;

			solution::f_calls = 0;
			opt = penalty(x0, c_wew, a, epsilon, Nmax, true);
			wewnetrzna << opt.x << sqrt(pow(opt.x(0), 2) + pow(opt.x(1), 2)) << "," << opt.y(0) << "," << solution::f_calls << endl;

			cout << i << endl;
		}

		fStartPoints.close();
		zewnetrzna.close();
		wewnetrzna.close();
		*/


		//cout << "LAB_5" << endl << endl;
		/*
		double epsilon = 1e3;
		int Nmax = 1000;
		matrix x0(2, 1);
		random_device rd;
		x0(0) = 20.0 * rd() / rd.max() - 10;
		x0(1) = 20.0 * rd() / rd.max() - 10;

		cout << "start point:" << endl << x0 << endl;

		solution optSD = SD(x0, epsilon, Nmax);

		cout << "SD" << endl;
		cout << optSD.x << endl << optSD.y << endl
			<< solution::f_calls << endl << solution::g_calls << endl;

		solution::f_calls = 0;
		solution::g_calls = 0;
		solution optN = Newton(x0, epsilon, Nmax);

		cout << "Newton" << endl;
		cout << optN.x << endl << optN.y << endl
			<< solution::f_calls << endl << solution::g_calls << endl << solution::H_calls << endl;

		solution::f_calls = 0;
		solution::g_calls = 0;
		solution optCG = CG(x0, epsilon, Nmax);

		cout << "CG" << endl;
		cout << optCG.x << endl << optCG.y << endl
			<< solution::f_calls << endl << solution::g_calls << endl;
		*/

		/*
		//STALOKROKOWA - 2 WARIANTY 
		//ZMIENNOKROKOWA 

		//0.05 0.12 stalo krokowe i zmienno krokowa
		//dla stalych - wywalic p i h - h okreslic na wartosc (0.05, albo 0.12)
		double epsilon = 1e-3;
		int Nmax = 1000;
		matrix x0(2, 1);
		random_device rd;

		solution opt;
		ofstream points("start_points.csv");
		ofstream results("results.csv");

		for (int i = 0; i < 100; i++) {
			x0(0) = 20.0*rd() / rd.max() - 10;
			x0(1) = 20.0*rd() / rd.max() - 10;
			points << endl << x0(0) << "," << x0(1) << endl << endl;

			solution::f_calls = solution::g_calls = 0;
			opt = SD(x0, epsilon, Nmax, 0.05);
			results << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << "," << solution::g_calls << ",";

			solution::H_calls = solution::f_calls = solution::g_calls = 0;
			opt = CG(x0, epsilon, Nmax, 0.05);
			results << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << "," << solution::g_calls << ",";

			solution::H_calls = solution::f_calls = solution::g_calls = 0;
			opt = Newton(x0, epsilon, Nmax, 0.05);
			results << opt.x(0) << ", " << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << "," << solution::g_calls << "," << solution::H_calls << ",";

			results << endl;

			solution::f_calls = solution::g_calls = 0;
			opt = SD(x0, epsilon, Nmax, 0.12);
			results << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << "," << solution::g_calls << ",";

			solution::H_calls = solution::f_calls = solution::g_calls = 0;
			opt = CG(x0, epsilon, Nmax, 0.12);
			results << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << "," << solution::g_calls << ",";

			solution::H_calls = solution::f_calls = solution::g_calls = 0;
			opt = Newton(x0, epsilon, Nmax, 0.12);
			results << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << "," << solution::g_calls << "," << solution::H_calls << ",";

			results << endl;

			solution::f_calls = solution::g_calls = 0;
			opt = SD(x0, epsilon, Nmax);
			results << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << "," << solution::g_calls << ",";

			solution::H_calls = solution::f_calls = solution::g_calls = 0;
			opt = CG(x0, epsilon, Nmax);
			results << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << "," << solution::g_calls << ",";

			solution::H_calls = solution::f_calls = solution::g_calls = 0;
			opt = Newton(x0, epsilon, Nmax);
			results << opt.x(0) << "," << opt.x(1) << "," << opt.y(0) << "," << solution::f_calls << "," << solution::g_calls << "," << solution::H_calls << ",";

			results << endl;
		}

		points.close();
		results.close();
		*/

		//cout << "LAB_6" << endl << endl;
		/*
		//f1-masa
		//f2-ugiecie
		//f = w*f1 + (1-w)*f2
		// zmienne decyzyjne x = [ l(dlugosc) d(srednica) ] , l C [0.2 , 1], d C [0.01, 0.05], C = nalezy do

		double epsilon = 1e-6;
		int Nmax = 1000;
		matrix x0(2, 1);
		random_device rd;
		x0(0) = 0.8 * rd() / rd.max() + 0.2;
		x0(1) = 0.04 * rd() / rd.max() + 0.01;

		cout << "x0:" << endl << x0 << endl << endl;

		double w = 1;
		solution opt = Powell(x0, epsilon, Nmax, w);

		cout << "opt:" << endl
			<< "x:" << endl << opt.x << endl
			<< "y:" << endl << opt.y << endl
			<< "f_calls:" << endl << solution::f_calls << endl;

		// dla w = 1 nie bierzemy pod uwage kryterium f2 ( nie uwzgledniamy ugiecia )
		// dla w = 0 nie bierzemy pod uwage kryterium f1

		//x:
		//dlugosc
		//srednica

		//y:
		//masa
		//ugiecie
		//naprezenie2d
		*/
		/*
		ofstream file("result_lab_6.csv");

		double epsilon = 1e-6;
		int NMax = 1000;
		matrix x0(2, 1);
		random_device rd;
		int k = 0;
		for (double w = 0; w <= 1.01; w += 0.01) {
			x0(0) = 0.8*rd() / rd.max() + 0.2;
			x0(1) = 0.04*rd() / rd.max() + 0.01;
			solution::f_calls = 0;
			solution opt = Powell(x0, epsilon, NMax, w);
			file << x0(0) * 1000 << "," << x0(1) * 1000 << ","
				<< opt.x(0) * 1000 << "," << opt.x(1) * 1000 << ","
				<< opt.y(0) << "," << opt.y(1) * 1000 << ","
				<< solution::f_calls << endl;
		}
		file.close();
		*/

		//cout << "LAB_7" << endl << endl;
		/*
		double epsilon = 1e-3;
		int N = 2;
		int Nmax = 5000;
		matrix lb(N, 1), ub(N, 1);
		lb(0) = lb(1) = 0;
		ub(0) = ub(1) = 1;
		solution opt = EA(N, lb, ub, epsilon, Nmax);

		cout << "x = "<< endl << opt.x << endl
			<< "y = " << endl << opt.y << endl << endl
			<< "f_calls = " << solution::f_calls << endl;

		// znalzenione rozwiazanie optymalne jest gdy rozwiazanie y jest mniejsze niz epsilon
		// do podsumowania bierzemy tylko znalezione rozwiazania globalne
		// wpisujemy do excela rozwiazania lokalne, zaznaczamy je ale nie uwzgledniamy jej do liczenia sredniej
		// f_calls: 20 - pierwsza iteracja, 40 - wywolan w kazdej iteracji

		// pierwsza wspolrzedna; sigma
		// druga wspolrzedna; sigma
		// y
		// f_calls
*/

		ofstream output("output_file.csv");

		calculationForSigma(0.01, output);
		calculationForSigma(0.1, output);
		calculationForSigma(1, output);
		calculationForSigma(10, output);
		calculationForSigma(100, output);
		output.close();
	}
	catch (char * EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
