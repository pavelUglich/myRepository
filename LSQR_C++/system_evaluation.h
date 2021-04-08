#pragma once
//#include <utility>
#include "boundary_value_problem.h"


// сделать универсальным, переделать конструктор

std::vector<std::function<double(double, const std::vector<double>&)>>
create_the_system(double kappa, const std::vector<std::function<double(double, double, const std::vector<double>&)>>& eqs);

/**
 * \brief класс, строящий СЛАУ
 */
class system_evaluation
{
	std::vector<std::function<double(double, double, const std::vector<double>&)>>
		_equations; // вектор правых частей уравнений системы
	// краевые условия на левом конце отрезка
	std::map<size_t, double> left_conditions;
	// краевые условия на правом конце отрезка
	std::map<size_t, double> right_conditions;
	// погрешность вычислений
	double _epsilon;

public:
	system_evaluation(
		std::vector<std::function<double(double, double, const std::vector<double>&)>> eqs, 
		double kappa,
		std::map<size_t, double> left_conditions,
		std::map<size_t, double> right_conditions,
		double epsilon = 0.1e-6);

	std::vector<double> evaluate_the_right_part(
		const std::vector<double>& points, size_t component) const;
	std::vector<double> evaluate_the_right_part(double min_kappa, 
		double max_kappa, size_t size, size_t component) const;
	std::vector<std::vector<double>> evaluate_the_matrix(
		const std::vector<double>& points_kappa, 
		const std::vector<double>& points_y, 
		const std::function<double(double, const std::vector<double>&, size_t)>&
		expression) const;
	std::vector<std::vector<double>> evaluate_the_matrix(
		double min_kappa, double max_kappa, size_t size,
		const std::function<double(double, const std::vector<double>&, size_t)>&
		expression) const;
	double equation(double x) const;
	std::vector<double> eigen_frequencies(double min, double max, double step) const;
	/*std::vector<std::vector<double>> evaluate_the_matrix(
		double min_kappa, double max_kappa, size_t size,
		const std::vector<const std::function<double(double, const std::
			vector<double>&, size_t)>&>& expressions) const;*/
};

