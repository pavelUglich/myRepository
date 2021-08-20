#include "MatrixSystem.h"
#include<vector>
#include<cassert>
#include "vector_ops.h"

double innerprod(const vector<double>& a, const vector<double>& b)
{
	assert(a.size() == b.size());
	double sum = 0;
	for (size_t i = 0; i < a.size(); i++) sum += a[i] * b[i];
	return sum;
}

/**
 * \brief конструктор
 * \param A матрица
 * \param b правая часть
 * \param step шаг
 * \param p параметр стабилизатора
 * \param left граничное условие слева
 * \param right граничное условие справа
 */
matrix_system::matrix_system(vector<vector<double>> A, const vector<double>& b,
	double step, double p, BoundaryCondition left, BoundaryCondition right) :
	Matrix(std::move(A)), right_part(b)
{
	// size = Matrix.front().size();
	rows = Matrix.size();
	columns = Matrix.front().size();
	stabilizer = { columns, step, p, left, right };
	multiply_ASinv();
	multiply_transpose_au();
	qpr();
	multiply_rx();
}

/**
 * \brief диагональ (геттер)
 * \return диагональ
 */
vector<double> matrix_system::diagonal() const
{
	return p1;
}

/**
 * \brief наддиагональ
 * \return наддиагональ
 */
vector<double> matrix_system::up_diagonal() const
{
	return p2;
}

/**
 * \brief правая часть
 * \return правая часть
 */
vector<double> matrix_system::rightPart() const
{
	return right_part;
}

/**
 * \brief умножение вектора на матрицу, обратную к Q
 * \param v вектор
 * \return новое значение вектора
 */
vector<double> matrix_system::multiply_qtu(const vector<double>& v)
{
	auto qtu = v;
	for (size_t i = 0; i < columns; i++)
	{
		if (i > rows)
		{
			break;
		}
		vector<double> a(rows - i);
		for (size_t j = 0; j < rows - i; j++) a[j] = Matrix[j + i][i];
		double sc = 0;
		for (size_t k = 0; k < rows - i; k++)
		{
			sc += a[k] * qtu[k + i];
		}
		for (size_t j = i; j < qtu.size(); j++) qtu[j] -= 2 * a[j - i] * sc;
	}
	return qtu;
}

/**
 * \brief Умножение справа на матрицу, обратную к матрице стабилизатора
 */
void matrix_system::multiply_ASinv()
{
	auto diagonal = stabilizer.diagonal();
	auto up_diagonal = stabilizer.Updiagonal();
	for (auto& i : Matrix) i[0] /= diagonal[0];
	for (size_t i = 1; i < Matrix[0].size(); i++)
		for (auto& j : Matrix)
		{
			j[i] -= up_diagonal[i - 1] * j[i - 1];
			j[i] /= diagonal[i];
		}
}

/**
 * \brief обнуление поддиагональных элементов столбца путем умножения слева на
 * матрицу отражения
 * \param k номер столбца
 * \return значение диагонального элемента матрицы
 */
void matrix_system::del_col(size_t k) // 
{
	if (k >= columns || k >= rows)
	{
		return;
	}
	p1.push_back(delete_the_column(this->Matrix, k));
}

/**
 * \brief Умножение справа на матрицу отражения
 * \param k номер строки
 * \return новое значение элемента наддиагонали
 */
void matrix_system::del_row(size_t k)
{
	if (k >= columns - 1 || k >= rows)
	{
		return;
	}
	const auto l = columns - k - 1;
	vector<double> av(l);
	for (size_t i = 0; i < l; i++) av[i] = Matrix[k][i + k + 1];
	av[0] -= norm(av);
	normalize(av);
	vector<double> vv(l);
	for (size_t i = 0; i < l; i++) { vv[i] = Matrix[k][i + k + 1]; }
	double sc = innerprod(vv, av);
	const double pp = Matrix[k][k + 1] - 2 * av[0] * sc;
	for (size_t i = k + 1; i < rows; i++)
	{
		for (size_t j = 0; j < l; j++) vv[j] = Matrix[i][j + k + 1];
		sc = vv * av;
		for (size_t j = k + 1; j < columns; j++)
			Matrix[i][j] -= 2 * av[j - k - 1] * sc;
	}
	for (size_t i = 0; i < l; i++) Matrix[k][i + k + 1] = av[i];
	p2.push_back(pp);
}

/**
 * \brief Преобразование правой части, умножение на транспонированную матрицу
 * СЛАУ;
 */
void matrix_system::multiply_transpose_au()
{
	vector<double> v(columns);
	for (size_t i = 0; i < columns; i++)
	{
		v[i] = 0;
		for (size_t j = 0; j < rows; j++)
			v[i] += Matrix[j][i] * right_part[j];
	}
	right_part = v;
}

/**
 * \brief Сведение матрицы СЛАУ к двухдиагональному виду
 */
void matrix_system::qpr()
{
	const auto size = std::min(rows, columns);
	p1.clear();
	p2.clear();
	for (size_t i = 0; i < size; i++)
	{
		del_col(i);
		del_row(i);
	}
}

/**
 * \brief Умножение вектора правой части на матрицу, обратную к ортогональной
 * матрице R;
 */
void matrix_system::multiply_rx()
{
	for (size_t i = 0; i < rows - 1; i++)
	{
		vector<double> av(columns);
		for (size_t j = i + 1; j < columns; j++)
			av[j] = Matrix[i][j];
		double sc = 0;
		for (size_t j = i + 1; j < columns; j++)
			sc += av[j] * right_part[j];
		for (size_t j = i + 1; j < columns; j++) right_part[j] -= 2 * av[j] * sc;
	}
}

/**
 * \brief Умножение вектора u на матрицу, обратную к ортогональной матрице R;
 * \param u вектор правой части
 */
void matrix_system::multiply_rtx(vector<double>& u)
{
	auto v = u;
	for (size_t i = 0; i < rows; i++)
	{
		const size_t l = rows - i;
		if (columns < l)
		{
			continue;
		}
		vector<double> a(columns - l);
		for (size_t j = 0; j < a.size(); j++)
			a[j] = Matrix[l - 1][j + l];
		double sc = 0;
		for (size_t j = 0; j < a.size(); j++) sc += a[j] * v[j + l];
		for (size_t j = l; j < columns; j++) v[j] -= 2 * a[j - l] * sc;
	}
	u = v;
}

/**
 * \brief Умножение вектора u на матрицу, обратную к матрице стабилизатора;
 * \param u вектор правой части
 */
void matrix_system::multiply_Sinv(vector<double>& u) const
{
	auto diagonal = stabilizer.diagonal();
	auto up_diagonal = stabilizer.Updiagonal();
	auto x = u;
	x[columns - 1] = u[columns - 1] / diagonal[columns - 1];
	for (size_t i = 1; i < columns; i++)
	{
		const auto j = columns - i - 1;
		x[j] = (u[j] - up_diagonal[j] * x[j + 1]) / diagonal[j];
	}
	u = x;
}
