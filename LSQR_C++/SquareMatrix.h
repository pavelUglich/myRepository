#pragma once
#include "Matrix.h"

template<class T>
class SquareMatrix : public Matrix<T>
{
	std::vector<T> transform_the_right_part(
		const std::vector<T>& rightPart, const std::vector<size_t>& p) const;
	std::vector<T> there_and_back_again(const std::vector<std::vector<T>>& lu,
		const std::vector<T>& rp) const;


public:
	std::vector<std::vector<T>> LUDecompose(std::vector<size_t>& p) const;

	/**
	 * \brief конструктор задаёт размеры квадратной матрицы (все конструкторы
	 * должны обращаться к конструкторам базового класса)
	 * \param rows количество строк
	 */
	SquareMatrix(size_t rows);

	/**
	 * \brief конструктор задаёт размеры квадратной матрицы и заполняет её
	 * значеним value
	 * \param rows количество строк
	 * \param value значение
	 */
	SquareMatrix(size_t rows, const T& value);

	/**
	 * \brief конструктор задаёт размеры квадратной матрицы и заполняет её
	 * значениями из вектора
	 * \param rows размер матрицы
	 * \param elements вектор значений
	 */
	SquareMatrix(size_t rows, const std::vector<std::vector<T>>& elements);

	SquareMatrix(const Matrix<T>& matrix);


	/**
	 * \brief определитель матрицы
	 * \return определитель матрицы
	 */
	T det() const;

	/**
	 * \brief решение системы линейных алгебраических уравнений
	 * \param rightPart вектор правой части
	 * \return решение СЛАУ
	 */
	std::vector<T> linSolve(const std::vector<T>& rightPart) const;

	/**
	 * \brief обращение матрицы
	 * \return обратная матрица
	 */
	SquareMatrix operator~() const;

	/**
	 * \brief деление на другую матрицу (умножение справа на матрицу, обратную
	 * к другой)
	 * \param squareMatrix
	 * \return
	 */
	SquareMatrix& operator/=(const SquareMatrix<T>& squareMatrix);
	SquareMatrix operator/(const SquareMatrix<T>& squareMatrix) const;
};

template <class T>
SquareMatrix<T>::SquareMatrix(size_t rows) : Matrix<T>(rows, rows)
{
}

template <class T>
SquareMatrix<T>::SquareMatrix(size_t rows, const T& value) :
	Matrix<T>(rows, rows, value)
{
}

template <class T>
SquareMatrix<T>::SquareMatrix(size_t rows,
	const std::vector<std::vector<T>>& elements) :
	Matrix<T>(rows, rows, elements)
{
}

template<class T>
SquareMatrix<T>::SquareMatrix(const Matrix<T>& matrix) :Matrix<T>(
	matrix.getRows(), matrix.getColumns(), matrix.getElements())
{
	if (matrix.getRows() != matrix.getColumns())
	{
		throw std::invalid_argument("");
	}
}

template<class T>
T SquareMatrix<T>::det() const
{
	std::vector<size_t> p(this->getRows());
	auto lu = LUDecompose(p);
	T result = 1;
	for (size_t i = 0; i < p.size(); i++)
	{
		result *= lu[i][i];
	}
	int parity = 0;
	for (size_t i = 0; i < p.size(); i++)
	{
		for (size_t ii = i + 1; ii < p.size(); ii++)
		{
			if (p[i] > p[ii])
			{
				++parity;
			}
		}
	}
	return result * (parity % 2 == 0 ? 1 : -1);
}

template <class T>
std::vector<std::vector<T>> SquareMatrix<T>::LUDecompose(
	std::vector<size_t>& p) const
{
	auto elems = this->getElements();
	p.resize(this->getRows());
	for (size_t i = 0; i < this->getRows(); i++)
		p[i] = i;
	for (size_t i = 0; i < this->getRows(); i++) {
		T max = abs(elems[i][i]);
		size_t imax = i;
		for (size_t k = i; k < this->getRows(); k++) {
			if (abs(elems[k][i]) > max) {
				max = abs(elems[k][i]);
				imax = k;
			}
		}
		if (max < DBL_EPSILON) throw std::invalid_argument("");
		if (imax != i) {
			std::swap(p[i], p[imax]);
			std::swap(elems[i], elems[imax]);
		}
		for (size_t j = i + 1; j < this->getRows(); j++) {
			elems[j][i] /= elems[i][i];
			for (size_t k = i + 1; k < this->getRows(); k++)
				elems[j][k] -= elems[j][i] * elems[i][k];
		}
	}
	return elems;
}

template <class T>
std::vector<T> SquareMatrix<T>::transform_the_right_part(
	const std::vector<T>& rightPart, const std::vector<size_t>& p) const
{
	std::vector<T> rp(this->getRows());
	for (size_t i = 0; i < this->getRows(); i++)
	{
		rp[i] = rightPart[p[i]];
	}
	return rp;
}

template <class T>
std::vector<T> SquareMatrix<T>::there_and_back_again(
	const std::vector<std::vector<T>>& lu,
	const std::vector<T>& rp) const
{
	std::vector<T> y(this->getRows());
	for (size_t i = 0; i < this->getRows(); i++)
	{
		T sum = 0;
		for (size_t j = 0; j < i; j++)
		{
			sum += y[j] * lu[i][j];
		}
		y[i] = rp[i] - sum;
	}
	for (int i = static_cast<int>(this->getRows() - 1); i >= 0; i--)
	{
		T sum = 0;
		for (size_t j = i + 1; j < this->getRows(); j++)
		{
			sum += y[j] * lu[i][j];
		}
		y[i] = (y[i] - sum) / lu[i][i];
	}
	return y;
}

template<class T>
std::vector<T> SquareMatrix<T>::linSolve(const std::vector<T>& rightPart) const
{
	std::vector<size_t> p(this->getRows());
	auto lu = LUDecompose(p);
	auto rp = transform_the_right_part(rightPart, p);
	return there_and_back_again(lu, rp);
}

template<class T>
SquareMatrix<T> SquareMatrix<T>::operator~() const
{
	SquareMatrix<T> result(this->getRows());
	std::vector<size_t> p(this->getRows());
	auto lu = LUDecompose(p);
	for (size_t i = 0; i < this->getRows(); i++)
	{
		std::vector<T> right_part(this->getRows(), 0);
		right_part[i] = 1;
		auto rp = transform_the_right_part(right_part, p);
		auto column = there_and_back_again(lu, rp);
		for (size_t ii = 0; ii < this->getRows(); ii++)
		{
			result[ii][i] = column[ii];
		}
	}
	return result;
}

template<class T>
SquareMatrix<T>& SquareMatrix<T>::operator/=(
	const SquareMatrix<T>& squareMatrix)
{
	*this *= ~squareMatrix;
	return *this;
}

template<class T>
SquareMatrix<T> SquareMatrix<T>::operator/(
	const SquareMatrix<T>& squareMatrix) const
{
	auto temp(*this);
	return temp /= squareMatrix;
}