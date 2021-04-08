#pragma once
#include <vector>
const size_t ITER_MAX = 2;

class lsqr
{
	size_t rows, columns;                    // число строк и столбцов
	double eps;                              // погрешность решения
	std::vector<std::vector<double>> matrix; // матрица системы

	std::vector<double> multiply_matrix_and_vector(const std::vector<double>&) const;
	std::vector<double> multiply_matrix_T_and_vector(
		const std::vector<double>&) const;
	std::vector<double> a_multiply_b_plus_u(const std::vector<double>&v,
	                                        const std::vector<double>&u) const;
	std::vector<double> at_multiply_b_plus_u(const std::vector<double>&,
		const std::vector<double>&) const;
public:
	lsqr(const std::vector<std::vector<double>>& vectors, double eps = 0.1e-6);
	std::vector<double> lin_solve(const std::vector<double>&) const;

};

