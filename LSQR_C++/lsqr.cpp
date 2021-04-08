#include "lsqr.h"
#include <numeric>
#include <stdexcept>
#include "vector_ops.h"
#include <iostream>

/**
 * \brief конструктор
 * \param vectors вектор элементов
 * \param eps погрешность
 */
lsqr::lsqr(const std::vector<std::vector<double>>& vectors, double eps) :
	rows(vectors.size()), eps(eps)
{
	if (!rows)
	{
		columns = 0;
		return;
	}
	columns = vectors[0].size();
	for (const auto& x : vectors)
	{
		if (x.size() != columns)
		{
			throw std::invalid_argument("");
		}
	}
	matrix = vectors;
}

/**
 * \brief решение СЛАУ
 * \param right_part правая часть СЛАУ
 * \return решение СЛАУ
 */
std::vector<double> lsqr::lin_solve(const std::vector<double>& right_part) const
{
	std::vector<double> u(right_part);
	double beta = normalize(u);

	auto v = multiply_matrix_T_and_vector(u);
	double alpha = normalize(v);

	std::vector<double> result(columns, 0);

	auto phi_bar = beta;
	auto rho_bar = alpha;

	auto w(v);
	auto f_norm = frobenius_norm(matrix);
	auto n_right = norm(right_part);
	for (size_t i = 0; i < ITER_MAX; i++)
	{
		u = matrix * v - alpha * u;
		beta = normalize(u);

		v = -beta * v;
		v = at_multiply_b_plus_u(u, v);
		alpha = normalize(v);

		const auto rho = sqrt(rho_bar * rho_bar + beta * beta);
		const auto c = rho_bar / rho;
		const auto s = beta / rho;
		const auto theta = s * alpha;
		rho_bar = -c * alpha;
		const auto phi = c * phi_bar;
		phi_bar *= s;

		const auto t1 = phi / rho;
		const auto t2 = -theta / rho;
		const auto x_old = result;
		for (size_t ii = 0; ii < columns; ii++)
		{
			result[ii] += t1 * w[ii];
			w[ii] = t2 * w[ii] + v[ii];
		}
		/*const auto x_norm = norm(result) * f_norm / n_right;
		std::cout << i << " ";
		if (x_norm > 1.0)
		{
			result = x_old;
			break;
		}*/
	}
	return result;
}

/**
 * \brief перемножить матрицу и вектор
 * \param vector вектор
 * \return произведение
 */
std::vector<double> lsqr::multiply_matrix_and_vector(
	const std::vector<double>& vector) const
{
	if (columns != vector.size())
	{
		throw std::invalid_argument("");
	}
	return matrix * vector;
}

/**
 * \brief перемножение транспонированной матрицы и вектора
 * \param vector вектор
 * \return произведение
 */
std::vector<double> lsqr::multiply_matrix_T_and_vector(
	const std::vector<double>& vector) const
{
	if (rows != vector.size())
	{
		throw std::invalid_argument("");
	}
	std::vector<double> result(columns, 0);
	for (size_t i = 0; i < columns; i++)
	{
		for (size_t ii = 0; ii < rows; ii++)
		{
			result[i] += matrix[ii][i] * vector[ii];
		}
	}
	return result;
}

/**
 * \brief Расчёт произведения Av+u
 * \param v вектор
 * \param u вектор
 * \return Av+u результат
 */
std::vector<double> lsqr::a_multiply_b_plus_u(const std::vector<double>& v,
	const std::vector<double>& u) const
{
	return multiply_matrix_and_vector(v) + u;
}

/**
 * \brief Расчёт произведения ATv+u
 * \param v вектор
 * \param u вектор
 * \return Av+u результат
 */
std::vector<double> lsqr::at_multiply_b_plus_u(const std::vector<double>& u,
	const std::vector<double>& v) const
{
	return multiply_matrix_T_and_vector(u) + v;
}