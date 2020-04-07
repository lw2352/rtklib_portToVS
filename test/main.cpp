#include "main.h"
float ls_sqrt(int a)
{
	float X, x;
	X = 2 * a;
	x = a;
	/*y=x^2-a;
	dy/dx=2x;

	*/
	for (int i = 0; i < 10; i++)
	{
		X = (2 * x * x - (x * x - 2)) / (2 * x);
		x = X;
	}
	return X;
}
