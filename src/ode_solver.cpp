//Do not edit the code below (unless you know what you are doing)

#include "ode_solver.h"

matrix *solve_ode(double t0, double dt, double tend, const matrix &Y0, matrix *ud, matrix *ad)
{
	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	if (N < 2)
		throw "The time interval is defined incorrectly";
	int *s = get_size(Y0);
	if (s[1] != 1)
		throw "Initial condition must be a vector";
	int n = s[0];
	delete[] s;
	matrix *S = new matrix[2]{matrix(N, 1), matrix(n, N)};
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

matrix diff(double t, const matrix &Y, matrix *ud, matrix *ad)
{
#if LAB_NO == 1 && LAB_PART == 1

#elif LAB_NO == 1 && LAB_PART == 2

#elif LAB_NO == 2 && LAB_PART == 3

	matrix dY(3, 1);
	double a = 0.98, b = 0.63, g = 9.81, PA = 1.0, TA = 90, PB = 1.0, DB = 0.00365665, Fin = 0.01, Tin = 10, DA = (*ad)();
	double FAout = Y(0) > 0 ? a * b * DA * sqrt(2 * g * Y(0) / PA) : 0;
	double FBout = Y(1) > 0 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;
	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout; //check: ok z fa sie tu wlewa z fb si� wylewa i f wejsciowa
	dY(2) = Fin / Y(1) * (Tin - Y(2)) + FAout / Y(1) * (TA - Y(2));
	return dY;

#elif LAB_NO == 3 && LAB_PART == 3

	//alpha1 - szukane, pozostaly alphy - pochodne
	//masa ramienia, masa ciezarka, dlugosc ramienia, wsp tarncia, alfa referencyjne, omega referencyjne
	double mr = 1, mc = 10, l = 0.5, b = 0.5, a_ref = 3.14, o_ref = 0;
	//moment bezwladnosci
	double I = mr * l * l / 3 + mc * l * l;

	// wspo�czynniki z alg data
	double k1 = (*ad)(0), k2 = (*ad)(1);
	//moment sily przekladany do ramienia

	// Y(0) aktualna wartosc alfa
	// Y(1) predkosc katowa pochodna y(0) alfa
	double m = k1 * (a_ref - Y(0)) + k2 * (o_ref - Y(1));

	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (m - b * Y(1)) / I; //ok
	return dY;
#elif LAB_NO == 4 && LAB_PART == 2

	// dla 100 zew i wew kara , potem tylko zew
	double C = 0.47, r = 0.12, m = 0.6, ro = 1.2, g = 9.81;
	double S = 3.14 * r * r, omega = (*ud)(0);
	double DX = 0.5 * C * ro * S * Y(1) * abs(Y(1)), // y(1) -> Vx
			// si�a oporu mo�e przypiesza� pik� wi�c trzeba uwzgl�dni� znak pr�dko�ci
			DY = 0.5 * C * ro * S * Y(3) * abs(Y(3)),
				 FMx = 3.14 * ro * Y(3) * omega * pow(r, 3),
				 FMy = 3.14 * ro * Y(1) * omega * pow(r, 3);

	matrix dY(4, 1);
	dY(0) = Y(1);										 // x
	dY(1) = (-DX - FMx) / m;				 // x prim
	dY(2) = Y(3);										 // y
	dY(3) = (-m * g - DY - FMy) / m; // y prim
	return dY;
#elif LAB_NO == 7 && LAB_PART == 2
	//rownaia rozniczkowe z wozorow opisujacych ruch ciezarkow'
	double m1 = 5, m2 = 5, k1 = 1, k2 = 1, F = 1;
	double b1 = (*ud)(0), b2 = (*ud)(1);
	matrix dY(4, 1);
	dY(0) = Y(1);
	dY(1) = (-b1 * Y(1) - b2 * (Y(1) - Y(3)) - k1 * Y(0) - k2 * (Y(0) - Y(2))) / m1;
	dY(2) = Y(3);
	dY(3) = (F + b2 * (Y(1) - Y(3)) + k2 * (Y(0) - Y(2))) / m2;
	return dY;
#else
	matrix dY;
	return dY;
#endif
	return -1; // added to compile properly
}