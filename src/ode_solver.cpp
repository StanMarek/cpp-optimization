//Do not edit the code below (unless you know what you are doing)

#include"ode_solver.h"

matrix* solve_ode(double t0, double dt, double tend, const matrix& Y0, matrix* ud, matrix* ad)
{
	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	if (N < 2)
		throw "The time interval is defined incorrectly";
	int* s = get_size(Y0);
	if (s[1] != 1)
		throw "Initial condition must be a vector";
	int n = s[0];
	delete[]s;
	matrix* S = new matrix[2]{ matrix(N,1), matrix(n,N) };
	S[0](0) = t0;
	for (int i = 0; i < n; ++i)
		S[1](i, 0) = Y0(i);
	matrix k1(n, 1), k2(n, 1), k3(n, 1), k4(n, 1);
	for (int i = 1; i < N; ++i)
	{
		S[0](i) = S[0](i - 1) + dt;
		k1 = dt * diff(S[0](i - 1), S[1][i - 1], ud, ad);
		k2 = dt * diff(S[0](i - 1) + 0.5 * dt, S[1][i - 1] + 0.5 * k1, ud, ad);
		k3 = dt * diff(S[0](i - 1) + 0.5 * dt, S[1][i - 1] + 0.5 * k2, ud, ad);
		k4 = dt * diff(S[0](i - 1) + dt, S[1][i - 1] + k3, ud, ad);
		for (int j = 0; j < n; ++j)
			S[1](j, i) = S[1](j, i - 1) + (k1(j) + 2 * k2(j) + 2 * k3(j) + k4(j)) / 6;
	}
	S[1] = trans(S[1]);
	return S;
}

//You can edit the following code

matrix diff(double t, const matrix& Y, matrix* ud, matrix* ad)
{
#if LAB_NO==1 && LAB_PART==1
	matrix dY(2, 1);
	double m = 1, f = 2, k = 1, b = 0.5;
	dY(0) = Y(1);
	dY(1) = (f - b * Y(1) - k * Y(0)) / m;
	return dY;

#elif LAB_NO == 1 && LAB_PART == 2

#elif LAB_NO==2 && LAB_PART==3
	double a = 0.98, b = 0.63, g = 9.81, PA = 1, TA = 90, PB = 1, DB = 0.00365665, Fin = 0.01, Tin = 10, DA = (*ad)();
	double FAout = Y(0) > 0 ? a * b * DA * sqrt(2 * g * Y(0) / PA) : 0;
	double FBout = Y(1) > 0 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;
	matrix dY(3, 1);
	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = Fin / Y(1) * (Tin - Y(2)) + FAout / Y(1) * (TA - Y(2));
	return dY;
#elif LAB_NO==3 && LAB_PART==3
	double mr = 1, mc = 10, l = 0.5, b = 0.5, a_ref = 3.14, omg_ref = 0;
	double I = mr * l * l / 3 + mc * l * l, k1 = (*ad)(0), k2 = (*ad)(1);
	double M = k1 * (a_ref - Y(0)) + k2 * (omg_ref - Y(1));
	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (M - b * Y(1)) / I;
	return dY;
#elif LAB_NO==4 && LAB_PART==2
	double c = 0.47, r = 0.12, m = 0.6, ro = 1.2, g = 9.81;
	double s = 3.14 * r * r, Dx = 0.5 * c * ro * s * abs(Y(1)) * Y(1);
	double Dy = 0.5 * c * ro * s * abs(Y(3)) * Y(3);
	double FMx = 3.14 * ro * Y(3) * (*ud)(0) * pow(r, 3);
	double FMy = 3.14 * ro * Y(1) * (*ud)(0) * pow(r, 3);
	matrix dY(4, 1);
	dY(0) = Y(1);
	dY(1) = (-Dx - FMx) / m;
	dY(2) = Y(3);
	dY(3) = (-m * g - Dy - FMy) / m;
	return dY;

#elif LAB_NO==7 && LAB_PART==2

#else
	matrix dY;
	return dY;
#endif
}