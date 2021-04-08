#pragma once
#include<vector>
#include "Stabilizer.h"
//#include<boost/numeric/ublas/matrix.hpp>
using namespace std;

double norm(const vector<double> & v);
class MatrixSystem
{
	size_t size{};
	vector<vector<double>> Matrix;
	vector<double> right_part;
	vector<double> p1, p2;
	Stabilizer stabilizer;

	void multiply_ASinv();
	double DelCol(size_t k);
	double DelRow(size_t k);
	void MultiplyTransposeAu();
	void QPR();
	void multiply_Rx();



public:
	MatrixSystem();
	MatrixSystem(
		vector<vector<double>> A,
		const vector<double>& b,
		double step,
		double p = 1,
		BoundaryCondition left = Neumann,
		BoundaryCondition right = Neumann);
	~MatrixSystem();

	///<summary>
	///Диагональ двухдиагональной матрицы;
	///</summary>


	vector<double> Diagonal() const;

	///<summary>
	///Наддиагональ двухдиагональной матрицы;
	///</summary>
	vector<double> UpDiagonal() const;;

	///<summary>
	///Правая часть СЛАУ;
	///</summary>
	vector<double> rightPart() const;;

	vector<double> MultiplyQtu(const vector<double> & v);
	void multiply_Rtx(vector<double> &u);
	void multiply_Sinv(vector<double> &u);
};

