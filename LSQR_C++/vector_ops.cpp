#include "vector_ops.h"
#include <stdexcept>
#include <algorithm>
#include <numeric>

std::vector<double> operator+(const std::vector<double>& left,
	const std::vector<double>& right)
{
	if (left.size() != right.size())
	{
		throw std::invalid_argument("");
	}
	std::vector<double> result(left.size(), 0);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = left[i] + right[i];
	}
	return result;
}

std::vector<double> operator-(const std::vector<double>& left,
	const std::vector<double>& right)
{
	if (left.size() != right.size())
	{
		throw std::invalid_argument("");
	}
	std::vector<double> result(left.size(), 0);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = left[i] - right[i];
	}
	return result;
}

std::vector<double> operator*(double number,
	const std::vector<double>& vector)
{
	std::vector<double> result(vector);
	std::transform(result.begin(), result.end(),
		result.begin(), [number](auto x) {return x * number; });
	return result;
}

double operator*(const std::vector<double>& left, 
	const std::vector<double>& right)
{
	if (left.size() != right.size())
	{
		throw std::invalid_argument("");
	}
	double result = 0;
	for (size_t i = 0; i < left.size(); i++)
	{
		result += left[i] * right[i];
	}
	return result;
}

double norm(const std::vector<double>& vector)
{
	return sqrt(std::accumulate(vector.begin(), vector.end(),
		0.0, [](double x, double y) {return x + y * y; }));
}

double normalize(std::vector<double>& vector)
{
	auto _norm = norm(vector);
	if (abs(_norm) < DBL_EPSILON)
	{
		vector = std::vector<double>(vector.size(), 0);
		return 0.0;
	}
	std::transform(vector.begin(), vector.end(),
		vector.begin(), [_norm](auto x) {return x / _norm; });
	return _norm;
}

std::vector<double> operator*(const std::vector<std::vector<double>>& matrix,
	const std::vector<double>& vector)
{
	std::vector<double> result(matrix.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t ii = 0; ii < vector.size(); ii++)
		{
			result[i] += matrix[i][ii] * vector[ii];
		}
	}
	return result;
}

std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& matrix, double number)
{
	std::vector<std::vector<double>> result(matrix.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		std::vector<double> row(matrix[i].size());
		for (size_t ii = 0; ii < row.size(); ii++)
		{
			row[i] = matrix[i][ii] * number;
		}
		result[i] = row;
	}
	return result;
}

std::vector<double> operator*(const SquareMatrix<double>& matrix, const std::vector<double>& vector)
{
	return matrix.getElements() * vector;
}

double frobenius_norm(const std::vector<std::vector<double>>& matrix)
{
	double result = 0;
	for (const auto& x : matrix)
	{
		result += std::accumulate(std::begin(x), std::end(x), 0.0, 
			[](auto x, auto y) {return x + y * y; });
	}
	return sqrt(result);
}

double delete_the_column(std::vector<std::vector<double>>& matrix, size_t k)
{
	const auto l = matrix.size() - k;
	//av Образующий единичный вектор матрицы отражения
	std::vector<double> av(l);
	for (size_t i = 0; i < l; i++)
		av[i] = matrix[i + k][k]; //!!!
	av[0] -= norm(av);
	normalize(av);
	//vv Поддиагональная часть столбца матрицы
	std::vector<double> vv(l);
	for (size_t i = 0; i < l; i++) vv[i] = matrix[i + k][k];
	double sc = av * vv;//*/innerprod(av, vv);
	const double pp = matrix[k][k] - 2 * av[0] * sc;
	for (size_t i = k + 1; i < matrix.size(); i++)
	{
		for (size_t j = 0; j < l; j++) vv[j] = matrix[j + k][i];
		sc = av * vv;//*/ innerprod(av, vv);
		for (size_t j = k; j < matrix.size(); j++)
			matrix[j][i] -= 2 * av[j - k] * sc;
	}
	for (size_t i = 0; i < l; i++) matrix[i + k][k] = av[i];
	return pp;
}

std::vector<double> qr_algorithm(std::vector<std::vector<double>>& matrix)
{
	std::vector<double> result;
	for (size_t i = 0; i < matrix.size(); i++)
	{
		result.push_back(abs(delete_the_column(matrix, i)));
	}
	std::sort(result.begin(), result.end());
	return result;
}

