#include "MatrixSystem.h"
#include<vector>
#include<cassert>
#include<cmath>
#include "vector_ops.h"

double innerprod(const vector<double>& a, const vector<double>& b)
{
	assert(a.size() == b.size());
	double sum = 0;
	for (size_t i = 0; i < a.size(); i++) sum += a[i] * b[i];
	return sum;
}


/**
 * \brief 
 */
MatrixSystem::MatrixSystem()
= default;

MatrixSystem::MatrixSystem(vector<vector<double>> A, const vector<double>& b, 
	double step, double p, BoundaryCondition left, BoundaryCondition right):
	Matrix(std::move(A)), right_part(b)
{
	size = Matrix.front().size();
	stabilizer = Stabilizer(size, step, p, left, right);
	multiply_ASinv();
	MultiplyTransposeAu();
	QPR(); 
	multiply_Rx();
}

MatrixSystem::~MatrixSystem() = default;

vector<double> MatrixSystem::Diagonal() const
{
	return p1;
}

vector<double> MatrixSystem::UpDiagonal() const
{
	return p2;
}

vector<double> MatrixSystem::rightPart() const
{
	return right_part;
}

vector<double> MatrixSystem::MultiplyQtu(const vector<double>& v)
{
	assert(size == v.size());
	auto Qtu = v;
	for (size_t i = 0; i < size; i++)
	{
		vector<double> a(size - i);
		for (size_t j = 0; j < size - i; j++) a[j] = Matrix[j + i][i];
		double sc = 0;
		for (size_t k = 0; k < size - i; k++) sc += a[k] * Qtu[k + i];
		for (size_t j = i; j < size; j++) Qtu[j] -= 2 * a[j - i] * sc;
	}
	return Qtu;
}

///<summary>
///Умножение справа на матрицу, обратную к матрице стабилизатора 
///</summary>
void MatrixSystem::multiply_ASinv()
{
	auto Diagonal = stabilizer.diagonal();
	auto UpDiagonal = stabilizer.Updiagonal();
	for (size_t i = 0; i < Matrix.size(); i++) Matrix[i][0] /= Diagonal[0];
	for (size_t i = 1; i < size; i++)
		for (size_t j = 0; j < Matrix.size(); j++)
		{
			Matrix[j][i] -= UpDiagonal[i - 1] * Matrix[j][i - 1];
			Matrix[j][i] /= Diagonal[i];
		}
}

///<summary>
///Умножение слева на матрицу отражения
///</summary>
///<param name="k">k - номер столбца</param>
double MatrixSystem::DelCol(size_t k)
{
	return delete_the_column(this->Matrix, k);
}


///<summary>
///Умножение справа на матрицу отражения
///</summary>
///<param name="k">k - номер строки</param>
double MatrixSystem::DelRow(size_t k)
{
	const auto l = size - k - 1;
	vector<double> av(l);
	for (size_t i = 0; i < l; i++) av[i] = Matrix[k][i + k + 1];
	av[0] -= norm(av);
	normalize(av);
	vector<double> vv(l);
	for (size_t i = 0; i < l; i++) { vv[i] = Matrix[k][i + k + 1]; }
	double sc = innerprod(vv, av);
	double pp = Matrix[k][k + 1] - 2 * av[0] * sc;
	for (size_t i = k + 1; i < size; i++)
	{
		for (size_t j = 0; j < l; j++) vv[j] = Matrix[i][j + k + 1];
		sc = innerprod(vv, av);
		for (size_t j = k + 1; j < size; j++)
			Matrix[i][j] -= 2 * av[j - k - 1] * sc;
	}
	for (size_t i = 0; i < l; i++) Matrix[k][i + k + 1] = av[i];
	return pp;
}

///<summary>
///Проеобразование правой части, умножение на транспонированную матрицу СЛАУ;
///</summary>
void MatrixSystem::MultiplyTransposeAu()
{
	vector<double> v(size);
	for (size_t i = 0; i < size; i++)
	{
		v[i] = 0;
		for (size_t j = 0; j < size; j++)
			v[i] += Matrix[j][i] * right_part[j];
	}
	right_part = v;
}

///<summary>
///Сведение матрицы СЛАУ к двухдиагональному виду
///</summary>
void MatrixSystem::QPR()
{
	p1.resize(size);
	p2.resize(size);
	for (size_t i = 0; i < size - 2; i++)
	{
		p1[i] = DelCol(i);
		p2[i] = DelRow(i);
	}
	p1[size - 2] = DelCol(size - 2);
	p2[size - 2] = Matrix[size - 2][size - 1];
	p1[size - 1] = DelCol(size - 1);
	p2[size - 1] = 0;
}


///<summary>
///Умножение вектора правой части на матрицу, обратную к ортогональной матрице R;
///</summary>
void MatrixSystem::multiply_Rx()
{
	for (size_t i = 0; i < size - 1; i++)
	{
		vector<double> av(size);
		for (size_t j = i + 1; j < size; j++)
			av[j] = Matrix[i][j];
		double sc = 0;
		for (size_t j = i + 1; j < size; j++)
			sc += av[j] * right_part[j];
		for (size_t j = i + 1; j < size; j++) right_part[j] -= 2 * av[j] * sc;
	}
}

///<summary>
///Умножение вектора u на матрицу, обратную к ортогональной матрице R;
///</summary>
///<param name="u">вектор правой части</param>
void MatrixSystem::multiply_Rtx(vector<double>& u)
{
	auto v = u;
	for (size_t i = 0; i < size; i++)
	{
		const size_t l = size - i;
		vector<double> a(i);
		for (size_t j = 0; j < i; j++)
			a[j] = Matrix[l - 1][j + l];
		double sc = 0;
		for (size_t j = 0; j < i; j++) sc += a[j] * v[j + l];
		for (auto j = l; j < size; j++) v[j] -= 2 * a[j - l] * sc;
	}
	u = v;
}

///<summary>
///Умножение вектора u на матрицу, обратную к матрице стабилизатора;
///</summary>
///<param name="u">вектор правой части</param>
void MatrixSystem::multiply_Sinv(vector<double>& u)
{
	auto Diagonal = stabilizer.diagonal();
	auto UpDiagonal = stabilizer.Updiagonal();
	auto x = u;
	x[size - 1] = u[size - 1] / Diagonal[size - 1];
	for (size_t i = 1; i < size; i++)
	{
		auto j = size - i - 1;
		x[j] = (u[j] - UpDiagonal[j] * x[j + 1]) / Diagonal[j];
	}
	u = x;
}
