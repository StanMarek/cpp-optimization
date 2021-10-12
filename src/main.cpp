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
        double m = 10, b = 2.5, k = 100, F = 10;
        double t0 = 0, tend = 50, dt = 0.1;
        matrix Y0 = matrix(2, new double[2]{ 0,0 });
        matrix P = matrix(4, new double[4]{ m,b,k,F });
        matrix *Y = solve_ode(t0, dt, tend, Y0, &P);
        matrix OUT = hcat(Y[0], Y[1]);
        ofstream S("out_1_1.csv");
        S << OUT;
        S.close();
        cout << "m = " << P(0) << endl;
        cout << "b = " << P(1) << endl;
        cout << "k = " << P(2) << endl;
        cout << "F = " << P(3) << endl;


#elif LAB_NO==1 && LAB_PART==2
        double m = 10, b = 2.5, k = 100, FA = 10, Ff = 2;
double t0 = 0, tend = 50, dt = 0.1;
matrix Y0 = matrix(2, new double[2]{ 0,0 });
matrix P = matrix(5, new double[5]{ m,b,k,FA, Ff });
matrix *Y = solve_ode(t0, dt, tend, Y0, &P);
matrix OUT = hcat(Y[0], Y[1]);
ofstream S("out_1_2.csv");
S << OUT;
S.close();
cout << "m = " << P(0) << endl;
cout << "b = " << P(1) << endl;
cout << "k = " << P(2) << endl;
cout << "FA = " << P(3) << endl;
cout << "Ff = " << P(4) << endl;
#endif
    }
    catch (char * EX_INFO)
    {
        cout << EX_INFO << endl;
    }
    //system("pause");
    int c = getchar();
    return 0;
}
