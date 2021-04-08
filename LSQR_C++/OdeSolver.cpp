#include "OdeSolver.h"
#include "vector_ops.h"

std::vector<std::vector<double>> OdeSolver::evaluate_k(
	const std::vector<double>& initial_conditions, double a, double h)
{
	std::vector<std::vector<double>> result;
	std::vector<double> vector(_equations.size());
	for (size_t i = 0; i < _equations.size(); i++)
	{
		vector[i] = h * _equations[i](a, initial_conditions);
	}
	result.push_back(vector);
	for (size_t i = 0; i < _butcher_tableau.size() - 2; i++)
	{
		const double xx = a + _butcher_tableau[i][0] * h;
		std::vector<double> uh(initial_conditions);
		for (size_t ii = 0; ii < uh.size(); ii++)
		{
			for (size_t iii = 0; iii < result.size(); iii++)
			{
				uh[ii] += _butcher_tableau[i][iii + 1] * result[iii][ii];
			}
		}
		for (size_t ii = 0; ii < _equations.size(); ii++)
		{
			vector[ii] = h * _equations[ii](xx, uh);
		}
		result.push_back(vector);
	}
	return result;
}

std::pair<std::vector<double>, std::vector<double>> OdeSolver::evaluate_uh(
	const std::vector<double>& initial_conditions, double a, double h)
{
	std::vector<double> u(initial_conditions);
	std::vector<double> uh(initial_conditions);
	std::vector<std::vector<double>> k = evaluate_k(initial_conditions, a, h);
	const size_t size = _butcher_tableau.size() - 1;
	for (size_t i = 0; i < initial_conditions.size(); i++)
	{
		for (size_t ii = 1; ii < _butcher_tableau[size].size(); ii++)
		{
			u[i] += _butcher_tableau[size - 1][ii] * k[ii - 1][i];
			uh[i] += _butcher_tableau[size][ii] * k[ii - 1][i];
		}
	}
	return { u, uh };
}

double OdeSolver::calculate_the_residual(double a, 
	const std::vector<double>& initial_conditions, double h, 
	std::vector<double>& u)
{
	const auto uh = evaluate_uh(initial_conditions, a, h);
	u = uh.first;
	return norm(uh.first - uh.second);
}

OdeSolver::OdeSolver(
	std::vector<std::function<double(double, const std::vector<double>&)>>
	functions, double epsilon, tableau t) : _equations(std::move(functions))
{
	if (epsilon < 0)
	{
		throw std::invalid_argument("");
	}
	this->_epsilon = epsilon;
	this->_butcher_tableau = get_butcher_tableau(t);
}

std::vector<double> OdeSolver::solve(double a, double b,
	const std::vector<double>& initial_conditions)
{
	if (a > b)
	{
		throw std::invalid_argument("");
	}
	double h = b - a;
	std::vector<double> initials(initial_conditions);
	while (abs(b - a) > DBL_EPSILON)
	{
		std::vector<double> u_;
		double r = calculate_the_residual(a, initials, h, u_);
		while (abs(r) > _epsilon)
		{
			h /= 2;
			r = calculate_the_residual(a, initials, h, u_);
		}
		a += h;
		initials = u_;
	}
	return initials;
}

std::vector<std::vector<double>> OdeSolver::solve(
	const std::vector<double>& points,
	const std::vector<double>& initial_conditions)
{
	std::vector<std::vector<double>> result;
	result.push_back(solve(0, points[0], initial_conditions));
	for (size_t i = 1; i < points.size(); i++)
	{
		result.push_back(solve(points[i - 1], points[i], result[i - 1]));
	}
	return result;
}
