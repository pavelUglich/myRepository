#pragma once
#include <vector>
#include <functional>
#include <map>

class boundary_value_problem
{
	std::vector<std::function<double(double, const std::vector<double>&)>>
		_equations;                            // правые части ОДУ
	std::map<size_t, double> left_conditions;  // граничные условия слева
	std::map<size_t, double> right_conditions; // граничные условия справа
	double _epsilon; // погрешность
	std::vector<std::vector<double>> initial_conditions() const;
	std::vector<std::vector<double>> cauchy_problem_solutions() const;
	std::vector<double> get_the_initial_conditions(
		const std::vector<std::vector<double>>& solutions) const;

public:
	boundary_value_problem(
		std::vector<std::function<double(double,
			const std::vector<double>&)>> functions,
		std::map<size_t, double> left_conditions,
		std::map<size_t, double> right_conditions, double epsilon = 0.1e-6);
	std::vector<double> solve() const;
	std::vector<std::vector<double>> solve(
		const std::vector<double>& points) const;
	double determinant() const;
};

