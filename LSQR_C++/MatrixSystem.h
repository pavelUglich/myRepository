#pragma once
#include<vector>
#include "Stabilizer.h"
//using namespace std;

double norm(const vector<double> & v);
class matrix_system
{
	//size_t size;
	size_t rows;
	size_t columns;
	vector<vector<double>> Matrix;
	vector<double> right_part;
	vector<double> p1, p2;
	stabilizer stabilizer;

	void multiply_ASinv();
	void del_col(size_t k);
	void del_row(size_t k);
	void multiply_transpose_au();
	void qpr();
	void multiply_rx();



public:
	matrix_system(
		vector<vector<double>> A,
		const vector<double>& b,
		double step,
		double p = 1,
		BoundaryCondition left = Neumann,
		BoundaryCondition right = Neumann);

	vector<double> diagonal() const;
	vector<double> up_diagonal() const;
	vector<double> rightPart() const;
	vector<double> multiply_qtu(const vector<double> & v);
	void multiply_rtx(vector<double> &u);
	void multiply_Sinv(vector<double> &u) const;
};

