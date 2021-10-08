//Do not edit the code below (unless you know what you are doing)

#include"solution.h"

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
	if (solution::g_calls>0)
		S << "g_calls = " << solution::g_calls << endl;
	if (solution::H_calls)
		S << "H_calls = " << solution::H_calls << endl;
	return S;
}

//You can edit the following code

void solution::fit_fun(matrix *ud, matrix *ad)
{
	++f_calls;
#if LAB_NO==2 && (LAB_PART==1 || LAB_PART==2)
	
#elif LAB_NO==2 && LAB_PART==3
	
#elif LAB_NO==3 && (LAB_PART==1 || LAB_PART==2)
	
#elif LAB_NO==3 && LAB_PART==3
	
#elif LAB_NO==4 && LAB_PART==1
	
#elif LAB_NO==4 && LAB_PART==2
	
#elif LAB_NO==5 && (LAB_PART==1 || LAB_PART==2)
	
#elif LAB_NO==5 && LAB_PART==3
	
#elif LAB_NO==6 && LAB_PART==1
	
#elif LAB_NO==6 && LAB_PART==2
	
#elif LAB_NO==7 && LAB_PART==1
	
#elif LAB_NO==7 && LAB_PART==2
	
#endif
}

void solution::grad(matrix *ud, matrix *ad)
{
	++g_calls;
#if LAB_NO==5 && (LAB_PART==1 || LAB_PART==2)
	
#elif LAB_NO==5 && LAB_PART==3
	
#endif
}

void solution::hess(matrix *ud, matrix *ad)
{
	++H_calls;
#if LAB_NO==5 && (LAB_PART==1 || LAB_PART==2)
	
#endif
}
