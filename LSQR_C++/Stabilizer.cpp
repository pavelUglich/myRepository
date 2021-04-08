#include "Stabilizer.h"
#include<cmath>



void Stabilizer::SquareRoot()
{
	Diagonal[0] = sqrt(Diagonal[0]);
	UpDiagonal[0] /= Diagonal[0];
	for (size_t i = 1; i < size - 1; i++) {
		Diagonal[i] = sqrt(Diagonal[i] - UpDiagonal[i - 1] * UpDiagonal[i - 1]);
		UpDiagonal[i] /= Diagonal[i];
	}
	Diagonal[size - 1] = sqrt(Diagonal[size - 1] - UpDiagonal[size - 2] * UpDiagonal[size - 2]);
}


Stabilizer::Stabilizer()
= default;


Stabilizer::~Stabilizer()
= default;
