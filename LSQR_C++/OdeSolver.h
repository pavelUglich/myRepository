#pragma once
#include <utility>
#include <vector>
#include <functional>
#include "ButcherTableau.h"
#include <stdexcept>

enum tableau {
	HEUN, BOGACKI_SHAMPINE, CASH_CARP, RUNGE_KUTTA_FELDBERG,
	DORMAND_PRINCE, RKF_78
};

inline std::vector<std::vector<double>> get_butcher_tableau(tableau t)
{
	switch (t) {
	case HEUN: return heun;
	case BOGACKI_SHAMPINE: return bogacki_shampine;
	case CASH_CARP: return cash_carp;
	case RUNGE_KUTTA_FELDBERG: return runge_kutta_feldberg;
	case DORMAND_PRINCE: return dormand_prince;
	case RKF_78: return Rkf78;
	default: throw std::invalid_argument("");
	}
}

class OdeSolver
{
	std::vector<std::function<double(double, const std::vector<double>&)>>
		_equations;
	double _epsilon;
	std::vector<std::vector<double>> _butcher_tableau;

	
	std::vector<std::vector<double>> evaluate_k(
		const std::vector<double>& initial_conditions, double a, double h);
	std::pair<std::vector<double>, std::vector<double>> evaluate_uh(
		const std::vector<double>& initial_conditions, double a, double h);
	double calculate_the_residual(double a,
		const std::vector<double>& initial_conditions,
		double h, std::vector<double>& u);


public:
	OdeSolver(
		std::vector<std::function<double(double, const std::vector<double>&)>>
		functions, double epsilon, tableau t);
	std::vector<double> solve(double a, double b,
		const std::vector<double>& initial_conditions);
	std::vector<std::vector<double>> solve(const std::vector<double>& points,
		const std::vector<double>& initial_conditions);

};

