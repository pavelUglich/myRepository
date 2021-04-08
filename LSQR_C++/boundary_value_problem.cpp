#include "boundary_value_problem.h"
#include "OdeSolver.h"
#include "SquareMatrix.h"
#include "vector_ops.h"
#include "Parameters.h"

/**
 * \brief построение начальных условий вспомогательных задач Коши в виде
 * набора векторов
 * \return набор векторов
 */
std::vector<std::vector<double>>
boundary_value_problem::initial_conditions() const
{
	const auto size = _equations.size();
	const auto cols = size - left_conditions.size();
	std::vector<std::vector<double>> result(size);
	size_t counter = 0;
	for (size_t i = 0; i < result.size(); i++)
	{
		if (left_conditions.find(i) == left_conditions.end())
		{
			result[i] = std::vector<double>(cols, 0);
			result[i][counter++] = 1;
		}
		else
		{
			result[i] = std::vector<double>(cols);
			result[i][0] = left_conditions.at(i);
		}
	}
	return result;
}

/**
 * \brief построение решения вспомогательных задач Коши в виде набора векторов
 * \return набор векторов
 */
std::vector<std::vector<double>>
boundary_value_problem::cauchy_problem_solutions() const
{
	auto conditions = this->initial_conditions();
	const auto size = conditions[0].size();
	std::vector<std::vector<double>> solutions;
	OdeSolver ode_solver(_equations, this->_epsilon, RKF_78);
	for (size_t i = 0; i < size; i++)
	{
		std::vector<double> initials(_equations.size());
		for (size_t ii = 0; ii < initials.size(); ii++)
		{
			initials[ii] = conditions[ii][i];
		}
		auto solution = Parameters::kind == THIRD ?
			ode_solver.solve(Parameters::points, initials).back():
			ode_solver.solve(0, 1, initials);
		solutions.push_back(solution);
	}
	return solutions;
}

/**
 * \brief построение начальных условия для отыскания решения краевой задачи
 * \param solutions вектор решений вспомогательных задач Коши
 * \return вектор начальных условий
 */
std::vector<double> boundary_value_problem::get_the_initial_conditions(
	const std::vector<std::vector<double>>& solutions) const
{
	std::vector<std::vector<double>> matrix;
	std::vector<double> right_part;
	for (auto x : right_conditions)
	{
		right_part.push_back(x.second);
		std::vector<double> row(right_conditions.size());
		for (size_t i = 0; i < row.size(); i++)
		{
			row[i] = solutions[i][x.first];
		}
		matrix.push_back(row);
	}
	const SquareMatrix<double> slau_matrix(right_conditions.size(), matrix);
	return slau_matrix.linSolve(right_part);
}

/**
 * \brief конструктор
 * \param functions вектор правых частей системы ОДУ
 * \param left_conditions краевые условия слева
 * \param right_conditions краевые условия справа
 * \param epsilon погрешность
 */
boundary_value_problem::boundary_value_problem(
	std::vector<std::function<double(double, const std::vector<double>&)>> functions,
	std::map<size_t, double> left_conditions, 
	std::map<size_t, double> right_conditions, double epsilon):
	_equations(std::move(functions)),
	left_conditions(std::move(left_conditions)),
	right_conditions(std::move(right_conditions)),
	_epsilon(epsilon)
{
}

/**
 * \brief решение краевой задачи
 * \return  решение краевой задачи в виде вектора решения при y = 1
 */
std::vector<double> boundary_value_problem::solve() const
{
	auto solutions = cauchy_problem_solutions();
	auto initial_conditions = get_the_initial_conditions(solutions);
	std::vector<double> result(_equations.size(), 0);
	for (size_t i = 0; i < solutions.size(); i++)
	{
		result = result + initial_conditions[i] * solutions[i];
	}
	return result;
}

/**
 * \brief решение краевой задачи
 * \return решение краевой задачи в наборе точек points
 */
std::vector<std::vector<double>> boundary_value_problem::solve(
	const std::vector<double>& points) const
{
	const auto solutions = cauchy_problem_solutions();
	const auto initial_conditions = get_the_initial_conditions(solutions);
	std::vector<double> initials;
	size_t counter = 0;
	for (size_t i = 0; i < _equations.size(); i++)
	{
		if (left_conditions.find(i) == left_conditions.end())
		{
			initials.push_back(initial_conditions[counter++]);
		}
		else
		{
			initials.push_back(left_conditions.at(i));
		}
	}
	OdeSolver ode_solver(_equations, this->_epsilon, RKF_78);
	return ode_solver.solve(points, initials);
}

double boundary_value_problem::determinant() const
{
	const auto solutions = cauchy_problem_solutions();
	std::vector<std::vector<double>> matrix;
	for (auto x : right_conditions)
	{
		std::vector<double> row(right_conditions.size());
		for (size_t i = 0; i < row.size(); i++)
		{
			row[i] = solutions[i][x.first];
		}
		matrix.push_back(row);
	}
	const SquareMatrix<double> slau_matrix(right_conditions.size(), matrix);
	return slau_matrix.det();
}
