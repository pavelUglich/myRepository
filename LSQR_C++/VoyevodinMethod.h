#pragma once
#include<vector>
#include "MatrixSystem.h"


class VoyevodinMethod
{
	double AlphaInitialValue;//Ќачальное значение параметра регул€ризации
	double step; //длина отрезка разбиени€
	double h; //ѕогрешность оператора
	double delta;//ѕогрешность правой части
	vector<double> RightPart, p1, p2, Qtu;
	size_t size;
	vector<double> Solution;

public:


	VoyevodinMethod(const vector<vector<double>>& matrix,
		vector<double> rightpart, double Step, BoundaryCondition Left = Neumann, 
		BoundaryCondition Right = Neumann, double p = 1.0, 
		double alphaInitialValue = 0.1e-1, double H = 1.0e-4, double Delta = 1.0e-4,
		double eps = 0.1e-3);


	vector<double> solution() const;
};

