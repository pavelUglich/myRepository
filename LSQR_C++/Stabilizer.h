#pragma once
#include<vector>
using namespace std;

enum BoundaryCondition { Dirichle, Neumann };


class stabilizer
{
	///������ ������� �������������
	size_t size;
	///���������
	vector<double> diagonal_;
	///������������
	vector<double> up_diagonal_;	
	void square_root();

public:
	///������������
	stabilizer();

	stabilizer(size_t n, double step, double p, BoundaryCondition Left, BoundaryCondition Right);
	vector<double> diagonal() const { return diagonal_; }
	vector<double> Updiagonal() const { return up_diagonal_; }
};

