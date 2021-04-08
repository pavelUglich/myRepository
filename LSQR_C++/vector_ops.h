#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H
#include <vector>

#include "SquareMatrix.h"

std::vector<double> operator+(const std::vector<double>& left,
                              const std::vector<double>& right);

std::vector<double> operator-(const std::vector<double>& left,
	const std::vector<double>& right);

std::vector<double> operator*(double number, 
	const std::vector<double>& vector);

double operator*(const std::vector<double>& left, 
	const std::vector<double>& right);

double norm(const std::vector<double>& vector);

double normalize(std::vector<double>& vector);

std::vector<double> operator*(
	const std::vector<std::vector<double>>& matrix,
	const std::vector<double>& vector);

std::vector<double> operator*(
	const SquareMatrix<double>& matrix,
	const std::vector<double>& vector);


double frobenius_norm(const std::vector<std::vector<double>>& matrix);


double delete_the_column(std::vector<std::vector<double>>& matrix, size_t k);

std::vector<double> qr_algorithm(std::vector<std::vector<double>>& matrix);

#endif