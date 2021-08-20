#include "Stabilizer.h"
#include<cmath>



/**
 * \brief разложение трёхдиагональной матрицы методом квадратного корня
 */
void stabilizer::square_root()
{
	diagonal_[0] = sqrt(diagonal_[0]);
	up_diagonal_[0] /= diagonal_[0];
	for (size_t i = 1; i < size - 1; i++) {
		diagonal_[i] = sqrt(diagonal_[i] - up_diagonal_[i - 1] * up_diagonal_[i - 1]);
		up_diagonal_[i] /= diagonal_[i];
	}
	diagonal_[size - 1]
		= sqrt(diagonal_[size - 1] - up_diagonal_[size - 2] * up_diagonal_[size - 2]);
}


/**
 * \brief конструктор
 */
stabilizer::stabilizer()
= default;

/**
 * \brief конструктор
 * \param n размер
 * \param step шаг
 * \param p параметр стабилизатора
 * \param Left граничное условие слева 
 * \param Right граничное условие справа
 */
stabilizer::stabilizer(size_t n, double step, double p, BoundaryCondition Left, BoundaryCondition Right): size(n)
{
	diagonal_ = vector<double>(size);
	up_diagonal_ = vector<double>(size - 1);
	const double hStab = p / step / step;
	for (size_t i = 0; i < size; i++) diagonal_[i] = 1 + 2 * hStab;
	for (size_t i = 0; i < size - 1; i++) up_diagonal_[i] = -hStab;
	if (Left == Dirichle) diagonal_[0] = 1 + 3 * hStab;
	else diagonal_[0] = 1 + hStab;
	if (Right == Dirichle) diagonal_[size - 1] = 1 + 3 * hStab;
	else diagonal_[size - 1] = 1 + hStab;
	square_root();
}
