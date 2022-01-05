/***************************************************
Code written for the optimization exercises purposes
by Lukasz Sztangret, PhD
Department of Applied Computer Science and Modelling
AGH University of Science and Technology
***************************************************/

#include"opt_alg.h"
#include<random>

int main()
{
	try
	{
		std::cout << "LAB NUMBER " << LAB_NO << endl;
		std::cout << "LAB PART " << LAB_PART << endl << endl;
#if LAB_NO==0

#elif LAB_NO==1 && LAB_PART==1
		double t0 = 0, dt = 0.1, tend = 50;
		matrix Y0 = matrix(2, new double[2]{ 0,0 }); //pierwsza wartosc x1 od 0, druga wartosc predkosc od 0(??)
		matrix* Y = solve_ode(t0, dt, tend, Y0); //funkcja solve_ode zwraca dwie macierze. Pierwsza to czas, druga rozwi�zania w kroku czasowym
		matrix out = hcat(Y[0], Y[1]);
		ofstream sout("wyniki.csv");
		sout << out;
		sout.close();


#elif LAB_NO==1 && LAB_PART==2

#elif LAB_NO==2 && LAB_PART==1
		double x0 = -20, d = 1, alpha = 2, epsilon = 1e-5, gamma = 1e-200;
		int Nmax = 1000;
		double* p = expansion(x0, d, alpha, Nmax);
		cout << p[0] << " " << p[1] << endl;
		solution::clear_calls();
		solution opt_f = fib(p[0], p[1], epsilon);
		cout << "Fib{" << endl << opt_f << "}\n" << endl;
		solution::clear_calls();
		solution opt_lag = lag(p[0], p[1], epsilon, gamma, Nmax);
		cout << "Lag{\n" << opt_lag << "}" << endl;
		solution::clear_calls();
#elif LAB_NO==2 && LAB_PART==2
		double d = 1, alpha = 2.137, epsilon = 1e-5, gamma = 1e-200;
		int Nmax = 10000;
		ofstream sout("wyniki2_137.csv");
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> unif(-100.0, 100.0);
		double startingPoint = 0;
		for (int i = 0; i < 100; i++) {
			startingPoint = unif(gen);
			double* p = expansion(startingPoint, d, alpha, Nmax);
			int exp_fcall = solution::f_calls;
			solution::clear_calls();
			matrix ab_F(1, 1, 200);
			solution opt_f = fib(p[0], p[1], epsilon, &ab_F);
			int fib_fcall = solution::f_calls;
			solution::clear_calls();
			matrix ab_L(1, 1, 200);
			solution opt_lag = lag(p[0], p[1], epsilon, gamma, Nmax, &ab_L);
			int lag_fcall = solution::f_calls;
			solution::clear_calls();
			sout << startingPoint << ";" << p[0] << ";" << p[1] << ";" << exp_fcall << ";" << "fib" << ";" << opt_f.x << ";" << opt_f.y << ";" << fib_fcall << ";"
				<< "lag" << ";" << opt_lag.x << ";" << opt_lag.y << ";" << lag_fcall << std::endl;
		}
		sout.close();

		/*matrix ab_F(1, 1, 200);
		fib(-100, 100, 1e-5, &ab_F);
		std::cout << ab_F << endl;

		matrix ab_L(1, 1, 200);
		lag(-100, 100, 1e-5, 1e-200, 1000, &ab_L);
		std::cout << endl << ab_L << endl;*/
#elif LAB_NO==2 && LAB_PART==3
		double x0 = 0.005, d = 0.0001, alpha = 1.5, epsilon = 1e-5, gamma = 1e-200;
		int Nmax = 1000;
		ofstream souta("Lag.csv");
		ofstream soutb("Fib.csv");
		matrix* Ud = new matrix();
		//double* p = expansion(x0, d, alpha, Nmax, Ud);

		//matrix* UdF = new matrix();
		//test.fit_fun(Ud);
		solution opt_lag = lag(0.0001, 0.01, epsilon, gamma, Nmax, Ud);

		//cout << p[0] << "0" << p[1] << endl;

		//solution test(opt_lag.x);
		cout << opt_lag;
		opt_lag.fit_fun(Ud);
		souta << Ud[0];
		souta.close();
		solution::clear_calls();
		//std::cout << test << endl;
		cout << "\nTo drugie" << endl;
		solution opt_f = fib(0.0001, 0.01, epsilon, Ud);
		cout << opt_f;
		opt_f.fit_fun(Ud);
		soutb << Ud[0];
		soutb.close();
		delete Ud;
#elif LAB_NO==3 && LAB_PART==1
		//Testowa funkcja celu
		double s = 0.2137, alphaHJ = 0.5, alphaR = 2, beta = 0.5, epsilon = 1e-3;
		int Nmax = 5000;
		ofstream sout("HJ_02137.csv");
		ofstream sout2("Ros_02137.csv");

		matrix s0(2, 1, s);
		//matrix x0(2, 1, -0.1);
		for (int i = 0; i < 100; i++) {
			matrix x0 = 2 * rand_mat(2, 1) - 1;
			solution opt_HJ = HJ(x0, s, alphaHJ, epsilon, Nmax);
			sout << "x1;" << x0(0) << ";x2;" << x0(1) << ";HJ;" << opt_HJ.x(0) << ";" << opt_HJ.x(1) << ";" << opt_HJ.y
				<< ";" << solution::f_calls << std::endl;
			solution::clear_calls();
			solution opt_R = Rosen(x0, s0, alphaR, beta, epsilon, Nmax);
			sout2 << "x1;" << x0(0) << ";x2;" << x0(1) << ";Ros;" << opt_R.x(0) << ";" << opt_R.x(1) << ";" << opt_R.y
				<< ";" << solution::f_calls << std::endl;
			solution::clear_calls();
		}
		sout.close();
		sout2.close();
#elif LAB_NO==3 && LAB_PART==2
		double s = 0.1337, alphaHJ = 0.5, alphaR = 2, beta = 0.5, epsilon = 1e-3;
		int Nmax = 5000;
		ofstream sout("K2Wykres.csv");
		//matrix x0 = 2 * rand_mat(2, 1) - 1 ;
		matrix s0(2, 1, s);
		matrix x0(2, new double[2]{ -0.330306, 0.020948 });
		//matrix x0(2, 1, -0.1);
		std::cout << x0 << endl;
		matrix XS_HJ = trans(x0);
		matrix XS_R = trans(x0);

		solution opt_HJ = HJ(x0, s, alphaHJ, epsilon, Nmax, &XS_HJ);
		std::cout << opt_HJ;
		std::cout << XS_HJ;
		solution::clear_calls();
		solution opt_R = Rosen(x0, s0, alphaR, beta, epsilon, Nmax, &XS_R);
		std::cout << endl << opt_R;
		std::cout << XS_R;
		sout << "HJ\n" << XS_HJ << endl << "Ros\n" << XS_R;
		sout.close();
#elif LAB_NO==3 && LAB_PART==3
		double s = 0.5, alphaHJ = 0.5, alphaR = 2, beta = 0.5, epsilon = 1e-4;
		int Nmax = 100000;
		/*ofstream souta("K2RP_HJ.csv");
		ofstream soutb("K2RP_Ros.csv");*/

		matrix x0(2, 1, 5);
		matrix s0(2, 1, s);
		solution opt_HJ = HJ(x0, s, alphaHJ, epsilon, Nmax, nullptr);
		matrix Y0H(2, 1);
		matrix* Result_H = solve_ode(0, 0.1, 100, Y0H, nullptr, &opt_HJ.x);
		//souta << opt_HJ << endl;
		cout << opt_HJ << endl;
		//cout << Result_H[1] << endl;
		solution::clear_calls();
		solution opt_R = Rosen(x0, s0, alphaR, beta, epsilon, Nmax, nullptr);
		matrix Y0R(2, 1);
		matrix* Result_R = solve_ode(0, 0.1, 100, Y0R, nullptr, &opt_R.x);
		cout << opt_R << endl;
		//cout << Result_R[1] << endl;
		//soutb << opt_R << endl;
		solution::clear_calls();

#elif LAB_NO==4 && LAB_PART==1
		double c0 = 1, dc_out = 2, dc_in = 0.5, epsilon = 1e-3;
		int Nmax = 10000;
		matrix x0(2, 1), a = 5;
		ofstream sout1("K3_opt_zew_5_1.csv");
		ofstream sout2("K3_opt_wew_5_1.csv");
		for (int i = 0; i < 100; i++) {
			do {
				x0 = a * rand_mat(2, 1) + 1;
			} while (norm(x0) > a);
			cout << x0 << endl << endl;
			solution opt_zew = pen(x0, c0, dc_out, epsilon, Nmax, &a);
			sout1 << x0(0) << ";" << x0(1) << ";" << opt_zew.x(0) << ";" << opt_zew.x(1) << ";" << opt_zew.y << ";" << norm(opt_zew.x) << ";" << opt_zew.f_calls << endl;
			solution::clear_calls();

			solution opt_wew = pen(x0, c0, dc_in, epsilon, Nmax, &a);
			sout2 << x0(0) << ";" << x0(1) << ";" << opt_wew.x(0) << ";" << opt_wew.x(1) << ";" << opt_wew.y << ";" << norm(opt_wew.x) << ";" << opt_wew.f_calls << endl;
			solution::clear_calls();
		}

#elif LAB_NO==4 && LAB_PART==2
		double c0 = 1, dc = 4, epsilon = 1e-5;
		int Nmax = 10000;
		matrix x0(2, new double[2]{ 0, 0 });
		solution opt_zew = pen(x0, c0, dc, epsilon, Nmax);
		ofstream sout("K3_Real_Problem.csv");
		sout << opt_zew << endl;
		//std::cout << opt_zew.x(0) << std::endl;
		//std::cout << opt_zew.x(1) << std::endl;
		matrix Y0R(4, new double[4]{ 0,opt_zew.x(0), 100,0 });
		matrix ud(opt_zew.x(1));
		matrix* R = solve_ode(0, 0.01, 7, Y0R, &ud);
		sout << R[1];
		/*	matrix x0(2, 1, 2), c = 1;
		solution test(x0);
		test.fit_fun(nullptr, &c);
		cout << test << endl;*/


#elif LAB_NO==5 && LAB_PART==1
		
		// double h0 = 0.05;
		// double h0 = 0.12;
		double h0 = -1;
		double epsilon = 1e-5;
		int Nmax = 100000;

		ofstream soutSD("sd.xlsx");
		ofstream soutCG("cg.xlsx");
		ofstream soutN("newton.xlsx");

		for(int i = 0; i < 100; i++) {
			matrix x0 = 20 * rand_mat(2, 1) - 10;
			
			solution optSD, optCG, optN;
			optSD = SD(x0, h0, epsilon, Nmax);
			// cout << optSD << endl;
			// cout << "solution.x(0) " << optSD.x(0) << "\tsolution.x(1) " << optSD.x(1) << endl;
			// cout << "solution.y " << optSD.y << endl;
			// cout << "solution.f_calls " << optSD.f_calls << "\tsolution.g_calls " << optSD.g_calls << endl;
			soutSD << x0(0) << ";" << x0(1) << ";" << optSD.x(0) << ";" << optSD.x(1) << ";" << optSD.y << optSD.f_calls << ";" << optSD.g_calls << ";\n"; 
			// cout << x0(0) << ";" << x0(1) << ";" << optSD.x(0) << ";" << optSD.x(1) << ";" << optSD.y << optSD.f_calls << ";" << optSD.g_calls << ";\n"; 
			solution::clear_calls();

			optCG = CG(x0, h0, epsilon, Nmax);
			soutCG << optCG.x(0) << ";" << optCG.x(1) << ";" << optCG.y << optCG.f_calls << ";" << optCG.g_calls << ";\n";
			cout << optCG.x(0) << ";" << optCG.x(1) << ";" << optCG.y << optCG.f_calls << ";" << optCG.g_calls << ";\n";
			solution::clear_calls();

			optN = Newton(x0, h0, epsilon, Nmax);
			soutN << optN.x(0) << ";" << optN.x(1) << ";" << optN.y << optN.f_calls << ";" << optN.g_calls << ";" << optN.H_calls << ";\n";
			// cout << optN.x(0) << ";" << optN.x(1) << ";" << optN.y << optN.f_calls << ";" << optN.g_calls << ";" << optN.H_calls << ";\n";
			solution::clear_calls();

		}

#elif LAB_NO==5 && LAB_PART==2
		double h0 = -0.05; 
		// double h0 = -0.12; 
		//double h0 = -1;

		double epsilon = 1e-5;
		int Nmax = 100000;

		ofstream soutSD("sd2.xlsx");
		ofstream soutCG("cg2.xlsx");
		ofstream soutN("newton2.xlsx");

		matrix x0, *ud;

		x0 = rand_mat(2, 1) * 20 - 10;
		ud = new matrix(1, 2);
		(*ud).add_row(trans(x0));

		solution optSD, optCG, optN;

		optSD = SD(x0, h0, epsilon, Nmax, ud);
		soutSD << *ud;
		solution::clear_calls();

		ud = new matrix(1, 2);
		ud->add_row(trans(x0));
		optCG = CG(x0, h0, epsilon, Nmax, ud);
		soutCG << *ud;
		solution::clear_calls();

		ud = new matrix(1, 2);
		ud->add_row(trans(x0));
		optN = Newton(x0, h0, epsilon, Nmax, ud);
		soutN << *ud;

		
#elif LAB_NO==5 && LAB_PART==3

		// double h0 = 0.01;
		// double h0 = 0.001;
		double h0 = 0.0001;

		double epsilon = 1e-5;
		int Nmax = 100000;

		matrix x0(3, 1);

		ofstream realProblemOut("cg3.xlsx");
	
		solution optCG = CG(x0, h0, epsilon, Nmax);

		int m = 100;
		matrix X(3, m), Y(1, m);

		ifstream Xinput("XData.txt");
		Xinput >> X;
		Xinput.close();

		ifstream Yinput("YData.txt");
		Yinput >> Y;
		Yinput.close();

		double h, P = 0.;
		for (int i = 0; i < m; i++)
		{
			h = 1. / (1. + exp(-(trans(optCG.x) * X[i])()));
			if (lroundf(h) == Y(0, i))
				h = 1.;
			else
				h = 0.;
			P += h;
		}
		P /= m;

		realProblemOut << optCG.x(0) << ";" << optCG.x(1) << ";" << optCG.x(2) << ";" << optCG.y << ";" << P << ";" << solution::g_calls << endl;

	
#elif LAB_NO==6 && LAB_PART==1

#elif LAB_NO==6 && LAB_PART==2

#elif LAB_NO==7 && LAB_PART==1

#elif LAB_NO==7 && LAB_PART==2

#endif
	}
	catch (char* EX_INFO)
	{
		std::cout << EX_INFO << endl;
	}
	// system("pause");
	return 0;
}
