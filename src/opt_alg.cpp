#include "opt_alg.h"
#if LAB_NO > 1
double *expansion(double x0, double d, double alpha, int Nmax, matrix *ud, matrix *ad)
{
	double *p = new double[2];
	solution X0(x0), X1(x0 + d);
	X0.fit_fun(ud, ad);
	X1.fit_fun(ud, ad);
	if (X0.y == X1.y)
	{
		p[0] = X0.x();
		p[1] = X1.x();

		if (ud != nullptr)
		{
			double **Array = new double *[1];
			//Array[0] = new double[] { x0, p[0], p[1], (double)solution::f_calls };
			ud->add_row(matrix(1, 4, Array));
		}

		return p;
	}
	if (X0.y < X1.y)
	{
		d *= -1;
		X1.x = X0.x + d;
		X1.fit_fun(ud, ad);
		if (X0.y <= X1.y)
		{
			p[0] = X1.x();
			p[1] = X0.x() - d;

			//if (ud != nullptr)
			//{
			//	double** Array = new double* [1];
			//	//Array[0] = new double[] { x0, p[0], p[1], (double)solution::f_calls };
			//	ud->add_row(matrix(1, 4, Array));
			//}

			return p;
		}
	}
	solution X2;
	int i = 1;
	while (true)
	{
		X2.x = x0 + pow(alpha, i) * d; //x0 - double(faktyczny punkt startowy)
		X2.fit_fun(ud, ad);

		//if (ud != nullptr)
		//{
		//	double** Array = new double* [1];
		//	//Array[0] = new double[] { x0, X0.x(), X2.x(), (double)solution::f_calls };
		//	ud->add_row(matrix(1, 4, Array));
		//}

		if (X2.y >= X1.y || solution::f_calls > Nmax)
			break;
		X0 = X1;
		X1 = X2;
		++i;
	}
	d > 0 ? p[0] = X0.x(), p[1] = X2.x() : (p[0] = X2.x(), p[1] = X0.x());
	return p;
}

solution fib(double a, double b, double epsilon, matrix *ud, matrix *ad)
{
	int n = static_cast<int>(ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
	int *F = new int[n]{1, 1};
	for (int i = 2; i < n; ++i)
		F[i] = F[i - 2] + F[i - 1];
	solution A(a), B(b), C, D;
	C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
	D.x = A.x + B.x - C.x;
	C.fit_fun(ud, ad);
	D.fit_fun(ud, ad);
	for (int i = 0; i <= n - 3; ++i)
	{
		if (C.y < D.y)
			B = D;
		else
			A = C;
		C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
		D.x = A.x + B.x - C.x;
		C.fit_fun(ud, ad);
		D.fit_fun(ud, ad);
#if LAB_NO == 2 && LAB_PART == 2
		if (ud != nullptr)
		{
			(*ud).add_row((B.x - A.x)());
		}
#endif
	}
	return C;
}

solution lag(double a, double b, double epsilon, double gamma, int Nmax, matrix *ud, matrix *ad)
{
	solution A(a), B(b), C, D, D_old(a);
	C.x = (a + b) / 2;
	A.fit_fun(ud, ad);
	B.fit_fun(ud, ad);
	C.fit_fun(ud, ad);
	double l, m;
	while (true)
	{
		l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0) * (pow(A.x(0), 2) - pow(B.x(0), 2));
		m = A.y(0) * (B.x(0) - C.x(0)) + B.y(0) * (C.x(0) - A.x(0)) + C.y(0) * (A.x(0) - B.x(0));
		if (m <= 0) //check this: Ok, not positive
		{
			C.x = NAN;
			C.y = NAN;
			return C;
		}
		D.x = 0.5 * l / m;
		D.fit_fun(ud, ad);
		if (A.x <= D.x && D.x <= C.x) // a --- d ---- c ---- b
		{
			if (D.y < C.y)
			{
				B = C; // a --- c --- b
				C = D;
			}
			else
			{
				A = D; //a -- c --- b
			}
		}
		else if (C.x <= D.x && D.x <= B.x)
		{
			if (D.y < C.y)
			{
				A = C;
				C = D;
			}
			else
			{
				B = D;
			}
		}
		else
		{
			// punkt poza przedzialem ab
			C.x = NAN;
			C.y = NAN;
			return C;
		}
#if LAB_NO == 2 && LAB_PART == 2
		if (ud != nullptr)
		{
			(*ud).add_row((B.x - A.x)()); // jakie� dane usera
		}
#endif

		if (B.x - A.x < epsilon || abs(D.x() - D_old.x()) < gamma || solution::f_calls > Nmax)
		{
			return C;
		}
		D_old = D;
	}
}

#endif
#if LAB_NO > 2
solution HJ(matrix x0, double s, double alpha, double epsilon, int Nmax, matrix *ud, matrix *ad) //jesli etap probny konczy sie porazka to zmniejszamy s i zaczynamy szukac blizej srodka - zmiejszyc krok
{																																																 // kryterium stopu - osiagniecie dlugosci kroku mniejszej niz episilon - jesli jest to na odlegosc mniejszej niz epsilon nie ma mniejszej wartosci  | alfa < episolon ===> koniec
	solution XB(x0), XB_old, X;																																		 //XB - baza, XB_old - stara baza - odbicie symetryczne starej bazy wzgledem nowej bazy
	XB.fit_fun(ud, ad);

	if (alpha >= 1)
	{
		// alpha ma zmniejszać krok
		throw "Alpha shouldn't be more than 1";
	}

	while (true)
	{
		X = HJ_trial(XB, s, ud, ad); //etap próbny
		if (X.y < XB.y)							 //jesli konzcy sie sukcesem to kolejna petla, gdy etap probny zakonczyl sie sukcesem
		{
			while (true) // petla sukcesywnie wykonuje etapy robocze
			{
				XB_old = XB; //stara baza to XB, X_old na wykłądzie jest z podkreśleniem
				XB = X;			 //XB to jest X
#if LAB_NO == 3 && LAB_PART == 2
				if (ud != nullptr)
				{
					(*ud).add_row(trans(XB.x));
				}
#endif
				X.x = 2 * XB.x - XB_old.x; //Odbicie symetryczne xold przez xb
				X.fit_fun(ud, ad);
				X = HJ_trial(X, s, ud, ad); //etap probny w etpapie roboczym, z punktu x uruchamiam etap probny
				if (X.y >= XB.y)
				{
					break; //przerwe wykonywanie etapow roboczych gdy X bedzie gorszy niz Xb
				}
				if (solution::f_calls > Nmax)
				{
					return XB; //konczy caly algorytm, gdy liczba wywowal funkcji celu za duza
				}
			}
		}
		else
		{
			s *= alpha; //redukujemy dlugosc kroku, alfa zawsze <1
		}
		if (s < epsilon || solution::f_calls > Nmax) //sprawdzamy kryteria stopu
		{
			return XB;
		}
	}
}
// Póbkowanie - etap próbnya
solution HJ_trial(solution XB, double s, matrix *ud, matrix *ad) //moze zwrocic ten sam punkt albo lepszy
{
	int n = get_dim(XB);				//sprawdzamy rozmiar problemu
	matrix D = ident_mat(n);		//macierz zawierająca kierunki (macierz jednostkowa)
	solution X;									//rozwiązanie X
	for (int i = 0; i < n; ++i) //petla po kazdym kierunku
	{
		X.x = XB.x + s * D[i]; //XB + s w kierunku D
		X.fit_fun(ud, ad);		 //wartość funkcji celu w tym punkcie
		if (X.y < XB.y)				 //jezeli ten punkt jest lepszy
		{
			XB = X; //to w tym momencie się przesuwamy
		}
		else
		{
			X.x = XB.x - s * D[i]; //jesli nie jest lepszy - jest gorszy to stawiamy punkt trzeci z innej strony
			X.fit_fun(ud, ad);		 //krok w tym samym kierunku a przeciwnym zwrocie ^ ^ ^
			if (X.y < XB.y)				 //jesli pomaga to sie przesuwamy
			{
				// Ważne silne nierówności
				XB = X;
			}
		} // jesli nie to zostajemy w tym samym punkcie
	}
	return XB;
}
//s0 matrix - dla kadzego kroku oddzielna dlugosc kroku, //alhpa, beta - zwiekszanie/zmniejszanie dlugosci kroku
solution Rosen(matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix *ud, matrix *ad) //w tej metodzie tylko wylacznie etap probny, wystepuje sytuacje w ktorej zmienamy kierunki, moze byc sytuacja ze wykonujemy obwod
{																																																									// zwiekszamy i zmniejszamy dlugosc kroku, w HJ tylko zmniejszamy, | RS dlugosc korku zawsze inna
	solution X(x0), Xt;																																															//xt - punkt tymczasowy
	int n = get_dim(X);																																															//ilosc kierunkow - wymiar probleu
	matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);																																//macierz kierunkow, lambda - sumaryczne przesuniecie(l) - dla kazdego kroku oddzielnie, p - wektor wartosci 0,1,2,3... - liczy ile bylo porazek w danym kierunku, s - macierz dlugosci krokow kierunkow, D - kierunki
	X.fit_fun(ud, ad);
	while (true)
	{
		for (int i = 0; i < n; ++i) //proba w kazdym kierunku
		{
			Xt.x = X.x + s(i) * D[i];
			Xt.fit_fun(ud, ad);
			if (Xt.y < X.y) //jesli punkt lepszy
			{
				X = Xt;				 //dodaje do lambdy dlugosc kroku
				l(i) += s(i);	 //zwiekszam lambde o s
				s(i) *= alpha; //zwiekszam s, wazna kolejnosc, aplha powinna byc wieksza od 1
			}
			else
			{
				++p(i);
				s(i) *= -beta; //zmniejszam dlugosc kroku
			}
		}
#if LAB_NO == 3 && LAB_PART == 2
		if (ud != nullptr)
		{
			(*ud).add_row(trans(X.x));
		}
#endif
		bool change = true;
		for (int i = 0; i < n; ++i)		//zmiana bazy kierunkow, obracamy wektory
			if (l(i) == 0 || p(i) == 0) //obrot wykonujemy jezeli w kazdym kierunku jest przesuniecie i w kazdym kierunku jest porazka, wszystkie p musza byc mniejsze niz 0 i wwszystkie presuniecia(l) mniejsze od 0
			{
				change = false;
				break;
			}
		if (change) // nowe kierunki zaleza od starych kierunkow i l, musza byc prostopadle, wektory tych kierunkow musza byc znormlanizowane(dlugosc 1), uzywamy ortonoralizacji G-S,
		{
			matrix Q(n, n), v(n, 1); // wyzerowane, v wektor pionowy, macierz trójkątna dolna zawierająca wartości lambd
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j <= i; ++j) // wypelnianie macierzy Q lambdami
				{
					Q(i, j) = l(i);
				}
			}
			Q = D * Q;
			v = Q[0] / norm(Q[0]);			//pierwsza kolumna podzielona przez norme
			D.set_col(v, 0);						//wsatawiamy do macierzy D jako kierunek piewrszy - kolumne zerowa
			for (int i = 1; i < n; ++i) // nastepne kierunki liczmy w petli
			{
				matrix temp(n, 1);
				for (int j = 0; j < i; ++j)
				{
					temp = temp + trans(Q[i]) * D[j] * D[j]; //transponowana kolumna j razy D(kierunki) ktore sa juz policzone
				}
				v = (Q[i] - temp) / norm(Q[i] - temp);
				D.set_col(v, i);
			}									//koniec obrotu, wszystkie nowe kierunki
			s = s0;						//zerujemy s
			l = matrix(n, 1); //zerujemy lambdy
			p = matrix(n, 1); //zerujemy porazki
		}
		double max_s = abs(s(0));
		for (int i = 1; i < n; ++i) //kryterkium kroku gdy jedno z s za male/duze abs
		{
			if (max_s < abs(s(i)))
			{
				max_s = abs(s(i));
			}
		}
		if (max_s < epsilon || solution::f_calls > Nmax)
		{
			return X;
		}
	}
}
#endif
#if LAB_NO > 3
solution pen(matrix x0, double c0, double dc, double epsilon, int Nmax, matrix *ud, matrix *ad) //funkcja kary
{
	double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
	solution X(x0), X1;
	matrix c(2, new double[2]{c0, dc}); // c0 - waga przy karze dc - odpowiada za wybor funkcji kary; zew, wew ; wysyłana jak ad
	while (true)
	{
		X1 = sym_NM(X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud, &c); //metoda symplex, odpal metode optymalizacji bez ograniczeń
		if ((norm(X.x - X1.x) < epsilon) || solution::f_calls > Nmax)					 //gdy przestanie sie zmieniac lub max liczba iteracji
			return X1;
		c(0) *= dc;
		X = X1; //rozpoczynamy symulacje dla zmienionej wartosci- punktu ktory znalezlismy
	}
}
/// <summary>
/// sympleks
/// </summary>
/// <param name="x0">pkt poczatkowe</param>
/// <param name="s"></param>
/// <param name="alpha">Przy odbiciu</param>
/// <param name="beta">zawezenie</param>
/// <param name="gamma">ekspansja</param>
/// <param name="delta">dedukcja</param>
/// <param name="epsilon"></param>
/// <param name="Nmax"></param>
/// <param name="ud"></param>
/// <param name="ad"></param>
/// <returns></returns>
solution sym_NM(matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	matrix D = ident_mat(n);			 //wektory w macierzy D
	int N = n + 1;								 // liczba wierzcholkow symleksu; sympleks ilośc zmiennych + 1
	solution *S = new solution[N]; //symplex
	S[0].x = x0;
	S[0].fit_fun(ud, ad);
	for (int i = 1; i < N; ++i)
	{
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun(ud, ad);
	}
	solution PR, PE, PN; //reflected, ecspansed, narrowed
	matrix pc;					 //srodek ciezkosci sympleksu
	int i_min, i_max;
	while (true)
	{
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i)
		{
			if (S[i_min].y > S[i].y)
				i_min = i;
			if (S[i_max].y < S[i].y)
				i_max = i;
		}
		pc = matrix(n, 1);
		for (int i = 0; i < N; ++i)
			if (i != i_max)
				pc = pc + S[i].x;
		pc = pc / (N - 1);										 //liczenie srodka c; albo n
		PR.x = pc + alpha * (pc - S[i_max].x); //odbicie
		PR.fit_fun(ud, ad);
		if (S[i_min].y < PR.y && PR.y < S[i_max].y)
			S[i_max] = PR;
		else if (PR.y < S[i_min].y)
		{
			PE.x = pc + gamma * (PR.x - pc);
			PE.fit_fun(ud, ad);
			if (PE.y < PR.y) //aceptujemy odbicie, ackeptujemy ten ktory sprawia ze symplex bd mniejsz
				S[i_max] = PE;
			else
				S[i_max] = PR;
		}
		else
		{
			PN.x = pc + beta * (S[i_max].x - pc);
			PN.fit_fun(ud, ad);
			if (PN.y < S[i_max].y)
				S[i_max] = PN;
			else
			{
				for (int i = 0; i < N; ++i)
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x + S[i_min].x);
						S[i].fit_fun(ud, ad);
					}
			}
		}
		double max_s = norm(S[i_min].x - S[0].x); //odleglosc pomiedzy wierzcholkiem 0 i minimalnym
		for (int i = 1; i < N; ++i)								//liczymy najdluzsza krawedz symplexu
			if (max_s < norm(S[i_min].x - S[i].x))
				max_s = norm(S[i_min].x - S[i].x);
		if (max_s < epsilon || solution::f_calls > Nmax) //kryteria stopu
			return S[i_min];
	}
}
#endif
#if LAB_NO > 4
solution SD(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad) // punkt startowy, dlugosc kroku(jesli zmiennokrokowa <0)
{
	int n = get_len(x0);								// rozmiar problemu
	solution X, X1;											// rozwiazania
	X.x = x0;														// punkt startowy
	matrix d(n, 1), *P = new matrix[2]; // okreslenie kierunku
	solution h;													// dlugosc kroku ktora wyznaczamy
	double *ab;													// przedzial ab - aby podzielic ekspansje
	while (true)
	{
		X.grad();		// liczymy gradient
		d = -X.g;		// kierunek najszybszego spadku
		if (h0 < 0) //jestli to metoda zmiennokrokowa - znajdujemy dlugosc kroku
		{
			P[0] = X.x; // P zawiera X ora d - kierunek
			P[1] = d;
			ab = expansion(0, 1, 1.2, Nmax, ud, P); // wywolujemy metode ekspansji, P jako ad(wsk)
			h = golden(ab[0], ab[1], epsilon, Nmax, ud, P);
			X1.x = X.x + h.x * d; // przesuwamy sie
		}
		else // stalokrokowa
			X1.x = X.x + h0 * d;
#if LAB_NO == 5 && LAB_PART == 2
		(*ud).add_row(trans(X1.x));
#endif
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) // kryterium stopu
		{
			X1.fit_fun(ud, ad);
			return X1;
		}
		X = X1; // przesuwamy sie jesli nie ma konca
	}
}

solution CG(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n, 1), *P = new matrix[2];
	solution h;
	double *ab, beta;
	X.grad(); // liczymy gradient
	d = -X.g; // pierwszy gradient liczymy inaczej niż resztę, stąd przed pętlą
	while (true)
	{
		if (h0 < 0)
		{
			P[0] = X.x; //wysylamy P do fit_fun jako ad, potrzebne do g(alfa); wysylamy P do fit_fun jako ad, potrzebne do g(alfa)
			P[1] = d;
			ab = expansion(0, 1, 1.2, Nmax, ud, P);					//przedział początkowy - metoda ekspansji
			h = golden(ab[0], ab[1], epsilon, Nmax, ud, P); //liczymy za pomocą złotego podziału
			X1.x = X.x + h.x * d;
		}
		else // stalokrokowa
		{
			X1.x = X.x + h0 * d;
		}
#if LAB_NO == 5 && LAB_PART == 2
		(*ud).add_row(trans(X1.x));
#endif
		if (solution::f_calls > Nmax || solution::g_calls > Nmax || norm(X1.x - X.x) < epsilon) //kiedy się zatrzymać
		{
			X1.fit_fun(ud, ad);
			return X1;
		}
		X1.grad();
		beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
		d = -X1.g + beta * d; // nowy kierunek
		X = X1;								// zaczynamy z nowego punktu
	}
}

solution Newton(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n, 1), *P = new matrix[2];
	solution h;
	double *ab;
	while (true)
	{
		X.grad();						 // liczymy gradient
		X.hess();						 // liczymy hesian
		d = -inv(X.H) * X.g; // kierunek d
		if (h0 < 0)					 //jestli to metoda zmiennokrokowa - znajdujemy dlugosc kroku
		{
			P[0] = X.x; // P zawiera X ora d - kierunek
			P[1] = d;
			ab = expansion(0, 1, 1.2, Nmax, ud, P); // wywolujemy metode ekspansji, P jako ad(wsk)
			h = golden(ab[0], ab[1], epsilon, Nmax, ud, P);
			X1.x = X.x + h.x * d; // przesuwamy sie
		}
		else // stalokrokowa
			X1.x = X.x + h0 * d;
#if LAB_NO == 5 && LAB_PART == 2
		(*ud).add_row(trans(X1.x));
#endif
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) // kryterium stopu
		{
			X1.fit_fun(ud, ad);
			return X1;
		}
		X = X1; // przesuwamy sie jesli nie ma konca
	}
}

solution golden(double a, double b, double epsilon, int Nmax, matrix *ud, matrix *ad) // metoda zlotego podzialu
{
	double alfa = (sqrt(5) - 1) / 2;
	solution A, B, C, D;
	A.x = a; // poczatek
	B.x = b; // koniec przedzialu
	C.x = B.x - alfa * (B.x - A.x);
	C.fit_fun(ud, ad);
	D.x = A.x + alfa * (B.x - A.x);
	D.fit_fun(ud, ad);
	while (true)
	{
		if (C.y < D.y)
		{
			B = D;
			D = C;
			C.x = B.x - alfa * (B.x - A.x); // PW WHY ? Checkout: d6d1e5c74ab642797df877050f7274ca9d8fdfe6
			C.fit_fun(ud, ad);
		}
		else
		{
			A = C;
			C = D;
			D.x = A.x + alfa * (B.x - A.x);
			D.fit_fun(ud, ad);
		}
		if (B.x - A.x < epsilon || solution::f_calls > Nmax) //kryterium stopu
		{
			A.x = (A.x + B.x) / 2; //jeśli stop to zwracamy środek pomiędzy a i b
			A.fit_fun(ud, ad);
			return A;
		}
	}
}

#endif
#if LAB_NO > 5
solution Powell(matrix x0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	// a[0] = x
	// a[1] = d
	int n = get_len(x0);
	matrix D = ident_mat(n), *A = new matrix[2];
	solution X, P, h;
	X.x = x0;
	double *ab; // z ekspansji
	while (true)
	{
		P = X;
		for (int i = 0; i < n; ++i)
		{
			A[0] = P.x;
			A[1] = D[i];
			ab = expansion(0, 1, 1.2, Nmax, ud, A);
			h = golden(ab[0], ab[1], epsilon, Nmax, ud, A);
			P.x = P.x + h.x * D[i];
		}
		if (norm(X.x - P.x) < epsilon || solution::f_calls > Nmax)
		{
			P.fit_fun(ud);
			return P;
		}
		for (int i = 0; i < n - 1; ++i)
		{
			D.set_col(D[i + 1], i);
		}
		//ostatni kierunek daje sie inaczej
		D.set_col(P.x - X.x, n - 1);
		A[0] = P.x;
		A[1] = D[n - 1];
		ab = expansion(0, 1, 1.2, Nmax, ud, A);
		h = golden(ab[0], ab[1], epsilon, Nmax, ud, A);
		X.x = P.x + h.x * D[n - 1];
	}
}

#endif
#if LAB_NO > 6
solution EA(int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	solution *P = new solution[mi + lambda];
	solution *Pm = new solution[mi];
	random_device rd;
	default_random_engine gen;
	gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
	normal_distribution<double> distr(0.0, 1.0);
	matrix IFF(mi, 1), temp(N, 2);
	double r, s, s_IFF;
	double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5);
	int j_min;
	for (int i = 0; i < mi; ++i)
	{
		P[i].x = matrix(N, 2);
		for (int j = 0; j < N; ++j)
		{
			P[i].x(j, 0) = (limits(j, 1) - limits(j, 0)) * rand_mat(1, 1)() + limits(j, 0);
			P[i].x(j, 1) = sigma0(j);
		}
		P[i].fit_fun(ud, ad);
		cout << P[i].y << endl;
		if (P[i].y < epsilon)
		{
			return P[i];
		}
	}
	int WhileLoopCount = 0;
	while (true)
	{
		if (ud != nullptr)
		{
			for (size_t i = 0; i < mi; i++)
			{
				matrix data(1, 3, 10);
				data(0, 0) = WhileLoopCount;
				data(0, 1) = P[i].x(0, 0);
				data(0, 2) = P[i].x(1, 0);
				ud->add_row(data);
			}
		}
		WhileLoopCount++;
		cout << WhileLoopCount << endl;

		s_IFF = 0;
		for (int i = 0; i < mi; ++i)
		{
			IFF(i) = 1 / P[i].y();
			s_IFF += IFF(i);
		}
		for (int i = 0; i < lambda; ++i)
		{
			r = s_IFF * rand_mat(1, 1)();
			s = 0;
			for (int j = 0; j < mi; ++j)
			{
				s += IFF(j);
				if (r <= s)
				{
					P[mi + i] = P[j];
					break;
				}
			}
		}

		for (int i = 0; i < lambda; ++i)
		{
			r = distr(gen);
			for (int j = 0; j < N; ++j)
			{
				P[mi + i].x(j, 1) *= exp(tau1 * r + tau * distr(gen));
				P[mi + i].x(j, 0) += P[mi + i].x(j, 1) * distr(gen);
			}
		}

		for (int i = 0; i < lambda; i += 2)
		{
			r = rand_mat(1, 1)();
			temp = P[mi + i].x;
			P[mi + i].x = r * P[mi + i].x + (1 - r) * P[mi + i + 1].x;
			P[mi + i + 1].x = r * P[mi + i + 1].x + (1 - r) * temp;
		}
		for (int i = 0; i < lambda; ++i)
		{
			P[mi + i].fit_fun(ud, ad);
			if (P[mi + i].y < epsilon)
				return P[mi + i];
		}
		for (int i = 0; i < mi; ++i)
		{
			j_min = 0;
			for (int j = 1; j < mi + lambda; ++j)
				if (P[j_min].y > P[j].y)
					j_min = j;
			Pm[i] = P[j_min];
			P[j_min].y = 1e10;
		}
		for (int i = 0; i < mi; ++i)
			P[i] = Pm[i];
		if (solution::f_calls > Nmax)
			return P[0];
	}
}
#endif
