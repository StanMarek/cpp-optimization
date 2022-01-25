//Do not edit the code below (unless you know what you are doing)

#include "solution.h"

int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

solution::solution(double L)
{
	x = L;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(const matrix &A)
{
	x = A;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(int n, double *A)
{
	x = matrix(n, A);
	g = NAN;
	H = NAN;
	y = NAN;
}

int get_dim(const solution &A)
{
	return get_len(A.x);
}

void solution::clear_calls()
{
	f_calls = 0;
	g_calls = 0;
	H_calls = 0;
}

ostream &operator<<(ostream &S, const solution &A)
{
	S << "x = " << A.x << endl;
	S << "y = " << A.y << endl;
	S << "f_calls = " << solution::f_calls << endl;
	if (solution::g_calls > 0)
		S << "g_calls = " << solution::g_calls << endl;
	if (solution::H_calls)
		S << "H_calls = " << solution::H_calls << endl;
	return S;
}

//You can edit the following code

void solution::fit_fun(matrix *ud, matrix *ad)
{
	++f_calls;
#if LAB_NO == 2 && (LAB_PART == 1 || LAB_PART == 2)

	//https://www.wolframalpha.com/input/?i=-cos%280.1*x%29+*+e%5E%28%28+-%28++0.1*x+-+2*3.14+%29%5E2+%29%29%2B0.002+*+%280.1*x%29%5E2

	y = -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(), 2);
#elif LAB_NO == 2 && LAB_PART == 3
	// objetos w pierwszym, drugim zbiorniku i temperatura w drugim
	//fit fun: majac x obliczyc y
	matrix Y0 = matrix(3, new double[3]{5, 1, 10}); // warunki początkowe
	matrix *Y = solve_ode(0, 1, 1000, Y0, ud, &x);	// x = da tj wielosc otworu w zbiorniku
	//Y[0] - time
	//Y[1] - 3x3 [VA VB TB]
	// znajdz max z tb tj najwieksza temp
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 0; i < n; i++)
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	y = abs(max - 50);

#elif LAB_NO == 3 && (LAB_PART == 1 || LAB_PART == 2)

	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;

#elif LAB_NO == 3 && LAB_PART == 3

	// liczenie całki
	// y[0] czas, y[1] -> alfa ; omega
	matrix Y0(2, 1);																//warunki poczatkowe
	matrix *Y = solve_ode(0, 0.1, 100, Y0, ud, &x); //wektor czasu, ad to coś co wykorzystujemy aby wysłeć k1 i k2
	int n = get_len(Y[0]);
	y = 0;
	double a_ref = 3.14, o_ref = 0;

	for (int i = 0; i < n; ++i)
	{
		//x12 -> k12

		// 10 we wzorze: najważniejszy dla nas jest obrót, waga -> wieksza liczba szybciej się obrócimy wiecej energii; wolniejszy obrót mniej energii mniejsza liczba

		y = y + 10 * pow(a_ref - Y[1](i, 0), 2) + pow(o_ref - Y[1](i, 1), 2) + pow(x(0) * (a_ref - Y[1](i, 0)) + x(1) * (o_ref - Y[1](i, 1)), 2);
	}
	y = 0.1 * y; // zamiast monożyc kazdorazwowo w pętli raz na koncu => dt = 0.1
#elif LAB_NO == 4 && LAB_PART == 1
	//metoda optymalizacji bez ograczen
	double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
	y = sin(arg) / arg; // sinus kardynalny
	//dodawnaie funkcji kary
	// ad(0) - c kara z mocą C
	// ad(1) -  wspolczynnik kary dc
	//return;
	if ((*ad)(1) > 1)
	{																					//zew
		if (-x(0) + 1 > 0)											//g1 naruszone - dorzucamy kare
			y = y + (*ad)(0) * pow(-x(0) + 1, 2); //ad(0) to c, dorzucamy g1^2
		if (-x(1) + 1 > 0)											//g2 naruszone - dorzucamy kare
			y = y + (*ad)(0) * pow(-x(1) + 1, 2);
		if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - (*ud)(0) > 0) //g3 - a (ze  wzoru) przesylamy poprzez user data z funkcji main
			y = y + (*ad)(0) * pow(sqrt(pow(x(0), 2) + pow(x(1), 2)) - (*ud)(), 2);
	}
	else
	{
		if (-x(0) + 1 > 0)
			y = 1e10;
		else
			y = y - (*ad)(0) / (-x(0) + 1);
		if (-x(1) + 1 > 0)
			y = 1e10;
		else
			y = y - (*ad)(0) / (-x(1) + 1);
		if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - (*ud)(0) > 0)
			y = 1e10;
		else
			y = y - (*ad)(0) / sqrt(pow(x(0), 2) + pow(x(1), 2)) - (*ud)(0);
		//dla a =4 powinno byc 3.99999999999999999 , przy zewnetrznej 4.0000000001........
	}
#elif LAB_NO == 4 && LAB_PART == 2
	//predkosc poczatkowa x
	matrix Y0(4, new double[4]{0, x(0), 100, 0});
	matrix cos = x(1);
	matrix *Y = solve_ode(0, 0.01, 7, Y0, &cos); // albo new matrix

	//ofstream file("tet.csv");

	//file << Y[1];

	//file.close();

	int n = get_len(Y[0]);
	int i50 = 0, i0 = 0;				// pokazują na nr wiersza
	for (int i = 1; i < n; ++i) // n = 701 // piersza iteracja bd zera więc od drugiego
	{
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50)) /// odleglosc mniejsza od 50 y, szukamy gdy y == 50
			i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2))) // szukamy zera
			i0 = i;
	}
	y = -Y[1](i0, 0); //wartosc funkcji celu
	//vx0 [-10, 10]
	// w [-20, 20]
	// x [4, 6]

	if (abs(x(0)) - 10 > 0)									// ograniczenie y
		y = y + (*ad)(0) * pow(x(0) - 10, 2); // F = f+ c*s , ad[0] tj c, S= potęga z naruszenia
	if (abs(x(1)) - 20 > 0)
		y = y + (*ad)(0) * pow(x(1) - 20, 2); // w
	if (abs(Y[1](i50, 0) - 5) - 1 > 0)			// x
		y = y + (*ad)(0) * pow(abs(Y[1](i50, 0) - 5) - 1, 2);

		//std::cout << "X na 50m: " << Y[1](i50, 0) << "\nX na 0m: " << Y[1](i0, 0) << endl;

		//ofstream sym_file("Sym.csv");

		//sym_file << Y[1];
		//sym_file.close();
		//std::printf("Saving... ");
		//sym_file.good() ? std::printf("Ok\n") : std::printf("Fail\n");

#elif LAB_NO == 5 && (LAB_PART == 1 || LAB_PART == 2)
	if (ad == nullptr) //nieprzekazujemy zadnego ad => czyli x ani kiernuku d cztli wzor z konspektu
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	else
	{

		solution temp;
		// g(alpha) = f(x + aplha*d)
		// ad 0 - x
		// ad 1 - d
		temp.x = ad[0] + x * ad[1];
		temp.fit_fun(ud);
		y = temp.y;
		--f_calls; // nie obl f celu a jest wykonywane na poczatku
	}
#elif LAB_NO == 5 && LAB_PART == 3
	// y = J
	// x = Q

	int m = 100; // liczba punktów uczących
	int n = get_len(x);
	static matrix X(n, m), Y(1, m);
	if (solution::f_calls == 1) // odczyt z pliku tylko przy pierwszym wywolaniu
	{
		ifstream S("XData.txt");
		S >> X;
		S.close();
		S.open("YData.txt");
		S >> Y;
		S.close();
	}
	// y = {0,1}
	// h0(xi) = [0,1]
	// log we wzorze jest po to aby łatwo liczyć gradient

	double h; // wartosc hipotezy
	y = 0;		// zerujemy y z klasy solution
	for (int i = 0; i < m; ++i)
	{
		h = (trans(x) * X[i])();															 // x -> teta(teza), X[i] -> zbior uczacy(oceny) //calosc jako wykladnik exponenty
		h = 1.0 / (1.0 + exp(-h));														 // wyliczenie hipotezy
		y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h); // zgodnie z wzorem, h zmienia sie od 0 do 1
																													 // log tj ln
	}
	y = y / m;
	// nie korzystamy z zmiennokrokowych

#elif LAB_NO == 6 && LAB_PART == 1

	//y = f(x); y[y1, y2] <- dwie f celu
	// y = g(alpha), y = w*f1(x + aplha*d) + (1-w) *f2(x+aplha*d)

	// w - ud[1]
	// a = ud[0]

	if (ad == nullptr)
	{
		y = matrix(2, 1);
		int a = ud[0]();
		y(0) = a * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
		y(1) = 1.0 / a * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
	}
	else
	{
		solution temp;
		temp.x = ad[0] + x * ad[1];
		temp.fit_fun(ud, nullptr); // zosatnie wywołana część powyższa

		double w = ud[1]();

		y = w * temp.y(0) + (1 - w) * temp.y(1);
		--f_calls;
	}

#elif LAB_NO == 6 && LAB_PART == 2
	if (ad == nullptr) //y=f(x) gdzie y to [y1//masa, y2//ugiecie, y3//naprezenie] f(x) to [l, d, ...];
	{
		y = matrix(3, 1);
		double ro = 7800, P = 1e3, E = 207e9;
		y(0) = ro * x(0) * 3.14 * pow(x(1), 2) / 4;										// masa
		y(1) = 64 * P * pow(x(0), 3) / (3 * E * 3.14 * pow(x(1), 4)); // ugiecie
		y(2) = 32 * P * x(0) / (3.14 * pow(x(1), 3));									// naprezenie
	}
	else // y=g(alpha)
	{
		solution T;
		T.x = ad[0] + x * ad[1];											// x+alfa*d - krok o x w kierunku ad[1]
		T.fit_fun(ud, nullptr);												// wskakujemy do 1 czesci
		matrix yn(2, 1);															// znormalizowane
		yn(0) = (T.y(0) - 0.12) / (3.06 - 0.12);			// masa znormalizowana
		yn(1) = (T.y(1) - 4.2e-5) / (0.026 - 4.2e-5); // znormalizowane ugiecie
		y = (*ud)() * yn(0) + (1 - (*ud)()) * yn(1);	// ud() - w w ud mamy tylko w w wcześniej w oraz a
		double c = 1e10;															// liczymy b. duże bo w f kary jak wcześniej w karze zew coraz wieksza była tak jak w ostatniej iteracji f kary
		if (T.x(0) < 0.2)															// l<200mm
		{
			y = y + c * pow(0.2 - T.x(0), 2);
		}
		if (T.x(0) > 1) // l>1m
		{
			y = y + c * pow(T.x(0) - 1, 2);
		}
		if (T.x(1) < 0.01) // d<10mm
		{
			y = y + c * pow(0.01 - T.x(1), 2);
		}
		if (T.x(1) > 0.05) // d>05mm
		{
			y = y + c * pow(T.x(1) - 0.05, 2);
		}
		if (T.y(1) > 0.005) // u>5mm
		{
			y = y + c * pow(T.y(1) - 0.005, 2);
		}
		if (T.y(2) > 300e6) // l>1m
		{
			y = y + c * pow(T.y(2) - 300e6, 2);
		}
		--f_calls;
	}
#elif LAB_NO == 7 && LAB_PART == 1

	++f_calls;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;

#elif LAB_NO == 7 && LAB_PART == 2

	int N = 1001;
	static matrix X(N, 2);
	if (solution::f_calls == 1)
	{
		ifstream S("/home/stan/Desktop/university/cpp-optimization/polozenia.txt");
		S >> X;
		S.close();
	}
	matrix Y0(4, new double[4]{0, 0, 0, 0});
	matrix *Y = solve_ode(0, 0.1, 100, Y0, &x[0]);
	y = 0;
	for (int i = 0; i < N; ++i)
	{
		y = y + abs(X(i, 0) - Y[1](i, 0)) + abs(X(i, 1) - Y[1](i, 2));
	}
	y = y / (2 * N);
#endif
}

void solution::grad(matrix *ud, matrix *ad)
{
	++g_calls;
#if LAB_NO == 5 && (LAB_PART == 1 || LAB_PART == 2)
	g = matrix(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
#elif LAB_NO == 5 && LAB_PART == 3
	int m = 100;
	int n = get_len(x);
	static matrix X(n, m), Y(1, m);
	if (solution::g_calls == 1) // odczyt z pliku tylko przy pierwszym wywolaniu - gradientu tym razem
	{
		ifstream S("XData.txt");
		S >> X;
		S.close();
		S.open("YData.txt");
		S >> Y;
		S.close();
	}
	double h; // hipoteza
	g = matrix(n, 1);
	for (int j = 0; j < n; ++j) // przechodzimy przez caly zbior liczacy
	{
		for (int i = 0; i < m; ++i)
		{
			h = (trans(x) * X[i])();
			h = 1.0 / (1.0 + exp(-h));
			g(j) = g(j) + X(j, i) * (h - Y(0, i)); // gradient
		}
		g(j) = g(j) / m;
	}
#endif
}

void solution::hess(matrix *ud, matrix *ad)
{
	++H_calls;
#if LAB_NO == 5 && (LAB_PART == 1 || LAB_PART == 2)
	H = matrix(2, 2);
	H(0, 0) = H(1, 1) = 10;
	H(0, 1) = H(1, 0) = 8;
#endif
}
