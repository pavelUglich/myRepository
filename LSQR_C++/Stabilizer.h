#pragma once
#include<vector>
using namespace std;

enum BoundaryCondition { Dirichle, Neumann };


class Stabilizer
{
	///������ ������� �������������
	size_t size{};
	///���������
	vector<double> Diagonal;
	///������������
	vector<double> UpDiagonal;
	/// ��������� ������������ ��������������� ������� �� �������� ������ ����������� �����
	void SquareRoot();

public:
	///������������
	Stabilizer();

	Stabilizer(size_t n, double step, double p, BoundaryCondition Left, BoundaryCondition Right) :size(n) {
		Diagonal = vector<double>(size);
		UpDiagonal = vector<double>(size - 1);
		const double hStab = p / step / step;
		for (size_t i = 0; i < size; i++)  Diagonal[i] = 1 + 2 * hStab;
		for (size_t i = 0; i < size - 1; i++)  UpDiagonal[i] = -hStab;
		if (Left == Dirichle) Diagonal[0] = 1 + 3 * hStab;
		else Diagonal[0] = 1 + hStab;
		if (Right == Dirichle) Diagonal[size - 1] = 1 + 3 * hStab;
		else Diagonal[size - 1] = 1 + hStab;
		SquareRoot();
	};
	~Stabilizer();
	vector<double> diagonal() const { return Diagonal; }
	vector<double> Updiagonal() const { return UpDiagonal; }
};

