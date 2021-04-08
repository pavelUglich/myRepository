#include "system_evaluation.h"

#include <utility>
#include "Parameters.h"

/**
 * \brief конструктор
 * \param eqs
 * \param kappa частота колебаний //!!!!
 * \param left_conditions краевые условия слева
 * \param right_conditions краевые условия справа
 * \param epsilon погрешность вычислений
 */
system_evaluation::system_evaluation(
	std::vector<std::function<double(double, double, const std::vector<double>&)>> eqs,
	double kappa,
	std::map<size_t, double> left_conditions,
	std::map<size_t, double> right_conditions, double epsilon) :
	_equations(std::move(eqs)),
	left_conditions(std::move(left_conditions)),
	right_conditions(std::move(right_conditions)),
	_epsilon(epsilon)
{
	//_equations// = create_the_system(kappa, eqs);
}

/**
 * \brief расчёт правой части СЛАУ
 * \param points вектор точек
 * \param component компонета решения
 * \return вектор правой части СЛАУ
 */
std::vector<double> system_evaluation::evaluate_the_right_part(
	const std::vector<double>& points, size_t component) const
{
	std::vector<double> result;
	for (auto x : points)
	{
		const auto kappa = x;
		const auto eq = create_the_system(kappa, this->_equations);
		boundary_value_problem bvp = { eq, left_conditions, right_conditions,
			_epsilon };
		result.push_back(bvp.solve()[component]);
	}
	return result;
}

/**
 * \brief Расчёт правой части СЛАУ
 * \param min_kappa минимальная частота
 * \param max_kappa максимальная частота
 * \param size размер вектора
 * \param component компонента решения
 * \return вектор правой части СЛАУ
 */
std::vector<double> system_evaluation::evaluate_the_right_part(
	double min_kappa, double max_kappa, size_t size, size_t component) const
{
	std::vector<double> points;
	const auto h = (max_kappa - min_kappa) / size;
	for (size_t i = 0; i <= size; i++) // !!!!
	{
		points.push_back(min_kappa + i * h);
	}
	return evaluate_the_right_part(points, component);
}


/**
 * \brief найти матрицу СЛАУ
 * \param points_kappa вектор точек по kappa
 * \param points_y вектор точек по y
 * \param expression
 * \return матрица СЛАУ
 */
std::vector<std::vector<double>> system_evaluation::evaluate_the_matrix(
	const std::vector<double>& points_kappa, const std::vector<double>& points_y,
	const std::function<double(double, const std::vector<double>&, size_t)>&
	expression) const
{
	std::vector<std::vector<double>> result(points_kappa.size(),
		std::vector<double>(points_y.size()));
	size_t i = 0;
	for (auto x : points_kappa)
	{
		const auto kappa = x;
		const auto eq = create_the_system(kappa, this->_equations);
		boundary_value_problem bvp = { eq, left_conditions, right_conditions,
			_epsilon };
		auto values = bvp.solve(points_y);
		for (size_t ii = 0; ii < points_y.size(); ii++)
		{
			result[i][ii] = expression(kappa, values[ii], ii);
		}
		++i;
	}
	return result;
}

/**
 * \brief строит вектор из значений, равноотстоящих на отрезке
 * \param min минимальное значение
 * \param max максимальное значение
 * \param size размер вектора
 * \return вектор
 */
std::vector<double> create_vector(double min, double max, size_t size)
{
	std::vector<double> result;
	const auto h = (max - min) / size;
	for (size_t i = 0; i < size; i++)
	{
		result.push_back(min + (i + 1) * h);
	}
	return result;
}

/**
 * \brief построить матрицу системы
 * \param min_kappa наименьшая частота
 * \param max_kappa наибольшая частота
 * \param size размер
 * \param expression
 * \return матрица
 */
std::vector<std::vector<double>> system_evaluation::evaluate_the_matrix(
	double min_kappa, double max_kappa, size_t size,
	const std::function<double(double, const std::vector<double>&, size_t)>&
	expression) const
{
	return evaluate_the_matrix(
		create_vector(min_kappa, max_kappa, size),
		create_vector(0, 1, size), expression);
}

double system_evaluation::equation(double x) const
{
	const auto eq = create_the_system(x, this->_equations);
	const boundary_value_problem bvp = { eq, left_conditions, right_conditions,
		_epsilon };
	return bvp.determinant();
}


double secant_method(double a, double b, const std::function<double(double)>& equation, double eps = 0.1e-6)
{
	double va = equation(a);
	double vb = equation(b);
	double cs = a;
	double cn = b;
	while (abs(cn - cs) > eps)
		if (va * vb < 0) {
			cs = cn;
			cn = b - (b - a) * vb / (vb - va);
			try
			{
				const double vc = equation(cn);
				if (va * vc < 0) {
					b = cn;
					vb = vc;
				}
				else {
					a = cn;
					va = vc;
				}
			}
			catch (...)
			{
				return cn;
			}			
		}
	return cn;
}

std::vector<double> system_evaluation::eigen_frequencies(double min,
	double max, double step) const
{
	std::vector<double> result;
	double a = min;
	double b = min + step;
	while (a < max)
	{
		if (equation(a) * equation(b) < 0)
		{
			result.push_back(secant_method(a, b, [=](double x) {return equation(x); }));
		}
		a = b;
		b += step;
	}
	return result;
}

/*
std::vector<std::vector<double>> system_evaluation::evaluate_the_matrix(double min_kappa, double max_kappa, size_t size,
																		const std::vector<const std::function<double(double, const std::
																			vector<double>&, size_t)>&>& expressions) const
{
	std::vector<std::vector<double>> result(size);
	for (const auto& expression : expressions)
	{
		auto matrix = evaluate_the_matrix(
			create_vector(min_kappa, max_kappa, size),
			create_vector(0, 1, size), expression);
		for (size_t i = 0; i < size; i++)
		{
			result[i].insert(result[i].end(), matrix[i].begin(), matrix[i].end());
		}
	}
	return result;
}*/

/**
 * \brief построение системы ОДУ
 * \param kappa частота
 * \param eqs
 * \return вектор правых частей ОДУ
 */
std::vector<std::function<double(double, const std::vector<double>&)>>
create_the_system(double kappa, const std::vector<std::function<double(double, double, const std::vector<double>&)>>& eqs)
{
	std::vector<std::function<double(double, const std::vector<double>&)>> result;
	result.reserve(eqs.size());
	for (const auto& eq : eqs)
	{
		result.emplace_back([=](double x, const std::vector<double>& vec) {return eq(x, kappa, vec); });
	}
	return result;
}
