/***************************************************
Code written for the optimization exercises purposes
by Lukasz Sztangret, PhD
Department of Applied Computer Science and Modelling
AGH University of Science and Technology
***************************************************/

#include"opt_alg.h"

int main()
{
	try
	{
		cout << "LAB NUMBER " << LAB_NO << endl;
		cout << "LAB PART " << LAB_PART << endl << endl;
#if LAB_NO==0
		
#elif LAB_NO==1 && LAB_PART==1
		double t0 = 0, dt = 0.1, tend = 50;
		matrix Y0 = matrix(2, new double[2]{ 0,0 }); //pierwsza wartosc x1 od 0, druga wartosc predkosc od 0(??)
		matrix *Y = solve_ode(t0, dt, tend, Y0); //funkcja solve_ode zwraca dwie macierze. Pierwsza to czas, druga rozwi�zania w kroku czasowym
		matrix out = hcat(Y[0], Y[1]);
		ofstream sout("wyniki.csv");
		sout << out;
		sout.close();


#elif LAB_NO==1 && LAB_PART==2

#elif LAB_NO==2 && LAB_PART==1

        ofstream sout("wynikiPart1.csv");

        srand(time(NULL));

        double x0 = -20;

        double d = 1;

        double epsilon = 1e-5;

        double alpha[3] = {1.0, 1.5, 2};

        int Nmax = 1000;

        double gamma = 1e-200;

        for(int a = 0; a < 3; a++){

            for(int i = 0; i < 100; i++) {

                double fMin = -100;

                double fMax = 100;

                double f = (double)rand() / RAND_MAX;

                double x0 = fMin + f * (fMax - fMin);

                double *p = expansion(x0, d, alpha[a], Nmax);

                sout << x0 << ";" << p[0] << ";" << p[1] << ";" << solution::f_calls << ";";

                solution::clear_calls();

                solution fibonacci = fib(p[0], p[1], epsilon);

                bool min = false;

                if(fibonacci.y < -0.5) {

                    min = true;

                }

                sout << fibonacci.x(0) << ";" << fibonacci.y(0) << ";" << solution::f_calls << ";" << min << ";";

                solution::clear_calls();

                solution lagrange = lag(p[0], p[1], epsilon, gamma, Nmax);

                min = false;

                if(lagrange.y < -0.5) {

                    min = true;

                }

                sout << lagrange.x(0) << ";" << lagrange.y(0) << ";" << solution::f_calls << ";" << min << "\n";

                solution::clear_calls();

            }

        }

        sout.close();

#elif LAB_NO==2 && LAB_PART==2

        double epsilon = 1e-5, gamma = 1e-200;

        int Nmax = 1000;

        matrix ab_F(1, 1, 200);

        solution opt_F = fib(-100, 100, epsilon, &ab_F);

        std::cout << ab_F;

        solution::clear_calls();

        cout << endl;

        matrix ab_L(1, 1, 200);

        solution opt_L = lag(-100, 100, epsilon, gamma, Nmax, &ab_L);

        std::cout << ab_L;

        solution::clear_calls();

#elif LAB_NO==2 && LAB_PART==3

#elif LAB_NO==3 && LAB_PART==1
		double s = 0.1, alphaHJ = 0.5, alphaR = 2, beta = 0.5, epsilon = 1e-3;
		int Nmax = 5000;
		//matrix x0 = 2 * rand_mat(2, 1) - 1;
		matrix s0(2, 1, s);
		matrix x0(2, 1, -0.1);

		solution optHJ = HJ(x0, s, alphaHJ, epsilon, Nmax);
		cout << optHJ << endl << endl;
		solution::clear_calls();
		solution optR = Rosen(x0, s0, alphaR, beta, epsilon, Nmax);
		cout << optR << endl << endl;
		solution::clear_calls();
		
#elif LAB_NO==3 && LAB_PART==2
		
#elif LAB_NO==3 && LAB_PART==3
		
#elif LAB_NO==4 && LAB_PART==1
		
#elif LAB_NO==4 && LAB_PART==2
		
#elif LAB_NO==5 && LAB_PART==1
		
#elif LAB_NO==5 && LAB_PART==2
		
#elif LAB_NO==5 && LAB_PART==3
		
#elif LAB_NO==6 && LAB_PART==1
		
#elif LAB_NO==6 && LAB_PART==2
		
#elif LAB_NO==7 && LAB_PART==1
		
#elif LAB_NO==7 && LAB_PART==2
		
#endif
	}
	catch (char * EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
