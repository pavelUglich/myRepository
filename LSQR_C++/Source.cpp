#include <vector>
#include <cstdlib>
#include <iostream>
#include "OdeSolver.h"
#include "boundary_value_problem.h"
#include "system_evaluation.h"
#include "Parameters.h"
#include "vector_ops.h"
#include "plots.h"
#include <chrono>
#include <iomanip>

#include "lsqr.h"
#include "Stabilizer.h"
#include "VoyevodinMethod.h"

//class VoyevodinMethod;
std::vector<std::function<double(double)>> Parameters::smooth_params;
std::vector<double> Parameters::const_params;
std::vector<double> Parameters::points;
std::vector<std::vector<double>> Parameters::piecewise_linear_params;
kind_of_solution Parameters::kind;

const double pi = 3.1415926538;

std::vector<std::vector<double>> rand_matrix(size_t n)
{
	std::vector<std::vector<double>> result(n, std::vector<double>(n, 0));
	for (size_t i = 0; i < n; i++)
	{
		for (size_t ii = 0; ii < n; ii++)
		{
			result[i][ii] = static_cast<double>(rand()) / RAND_MAX;
		}
	}
	return result;
}


std::vector<double> rand_vector(size_t n)
{
	std::vector<double> result(n, 0);
	for (size_t i = 0; i < n; i++)
	{
		result[i] = static_cast<double>(rand()) / RAND_MAX;
	}
	return result;
}


std::vector<double> landweber(const std::vector<std::vector<double>>& matrix,
	const std::vector<double>& vector, double omega)
{
	if (matrix.size() != vector.size())
	{
		throw std::invalid_argument("");
	}
	std::vector<double> x(vector.size(), 0);
	const auto rows = matrix[0].size();
	while (true)
	{
		auto y = matrix * x - vector;
		std::vector<double> result(y.size(), 0);
		for (size_t i = 0; i < result.size(); i++)
		{
			for (size_t ii = 0; ii < rows; ii++)
			{
				result[i] += matrix[ii][i] * y[ii];
			}
		}
		x = x - omega * result;
		if (norm(result) / norm(x) < 0.01)
		{
			break;
		}
	}
	return x;
}

std::function<double(double)> lagrange(const std::vector<double>& parameters)
{
	if (parameters.empty())
	{
		throw std::invalid_argument("");
	}
	const auto step = 1.0 / (parameters.size() - 1);
	std::vector<double> points(parameters.size());
	for (size_t i = 0; i < parameters.size(); i++)
	{
		points[i] = i * step;
	}
	return [=](double x)
	{
		double sum = 0;
		for (size_t i = 0; i < parameters.size(); i++)
		{
			double basis = 1.0;
			for (size_t ii = 0; ii < i; ii++)
			{
				basis *= (x - points[ii]) / (points[i] - points[ii]);
			}
			for (size_t ii = i + 1; ii < parameters.size(); ii++)
			{
				basis *= (x - points[ii]) / (points[i] - points[ii]);
			}
			sum += parameters[i] * basis;
		}
		return 1.0 + sum;
	};
}

std::vector<std::vector<double>> fill_up_the_matrix(size_t n, size_t m,
	const std::function<double(double, double)>& lambda,
	double min_kappa, double max_kappa, double a = 0, double b = 1)
{
	std::vector<std::vector<double>> result(n);
	const auto h_kappa = (max_kappa - min_kappa) / n;
	const auto h_x = (b - a) / m;
	for (size_t i = 0; i < n; i++)
	{
		std::vector<double> row(m);
		for (size_t ii = 0; ii < m; ii++)
		{
			row[ii] = lambda((i + 0.5) * h_kappa, (ii + 0.5) * h_x);
		}
		result[i] = row;
	}
	return result;
}

std::vector<double> fill_up_the_vector(size_t m,
	const std::function<double(double)>& lambda,
	double a = 0, double b = 1)
{
	std::vector<double> result(m);
	const auto h_x = (b - a) / m;
	for (size_t ii = 0; ii < m; ii++)
	{
		result[ii] = lambda((ii + 0.5) * h_x);
	}
	return result;
}

int main()
{
	setlocale(0, "");
	/*
		OdeSolver cauchy_problem =
		{
			{
				[](double x, const std::vector<double> & u) {return u[0]; }
			}
			,0.1e-8, RUNGE_KUTTA_FELDBERG
		};
		auto sol = cauchy_problem.solve(0, 1, { 1 });*/

	Parameters::kind = FIRST;

	// ����������� �������
	const double min_kappa = 2.5;
	const double max_kappa = 4.0;
	const double min_gamma = 2.9;
	const double max_gamma = 4.4;
	// ���������� �����
	const size_t points_x = 30;
	const size_t points_k = 10;
	const auto h_kappa = (max_kappa - min_kappa) / points_k;
	const auto h_gamma = (max_gamma - min_gamma) / points_k;
	// �������� �������
	std::vector<double> points_kappa;
	for (size_t i = 0; i <= points_k; i++)
	{
		points_kappa.push_back(min_kappa + i * h_kappa);
	}
	std::vector<double> points_gamma;
	for (size_t i = 0; i <= points_k; i++)
	{
		points_gamma.push_back(min_gamma + i * h_gamma);
	}
	std::vector<double> vv;
	const double h_x = 1.0 / points_x;
	for (size_t i = 0; i <= points_x; i++)
	{
		vv.push_back(i * h_x);
	}

	// ����� ���������
	Parameters::points = vv;
	Parameters::piecewise_linear_params = { points_x + 1, {1.0, 1.0} };
	Parameters::const_params = { 1.0, 1.0 };
	Parameters::smooth_params = {
		// 1
		//[](auto x) {return 2.0 - (x - 1) * (x - 1); },
		//[](auto x) {return 1.0 + exp(5 * (x - 1)); }
		// 2
		[](auto x) {return 2.0 - x * x; },
		[](auto x) {return 1.0 + exp(-4 * x + 0.02); }
		// 3
		//[](auto x) {return 1.0 + sin(pi * x); },
		//[](auto x) {return 2 * x * x - 2 * x + 1; }
	};

	std::vector<double> exact_solution_mu(points_x);
	std::vector<double> exact_solution_rho(points_x);
	for (size_t i = 0; i < points_x; i++)
	{
		exact_solution_mu[i] = Parameters::smooth_params[0](vv[i]);
		exact_solution_rho[i] = Parameters::smooth_params[1](vv[i]);
	}

	// ������� ������
	// ������� ���������
	std::function<double(double)> e = [](double x) { return Parameters::evaluate(x, 0); };
	std::function<double(double)> rho = [](double x) { return Parameters::evaluate(x, 1); };
	const std::vector<std::function<double(double, double, const std::vector<double>&)>> longitudinal = {
		[=](double x, double kappa, const std::vector<double>& vec) {return vec[1] / e(x); },
		[=](double x, double kappa, const std::vector<double>& vec) {return -kappa * kappa * rho(x) * vec[0]; }
	};
	const system_evaluation longitidinal_evaluation = { longitudinal, 1.0, { {0,0.0} }, { { 1,1.0 }} };

	const std::vector<std::function<double(double, double, const std::vector<double>&)>> flexural = {
		[=](double x, double kappa, const std::vector<double>& vec) {return vec[1]; },
		[=](double x, double kappa, const std::vector<double>& vec) {return vec[2] / e(x); },
		[=](double x, double kappa, const std::vector<double>& vec) {return vec[3]; },
		[=](double x, double kappa, const std::vector<double>& vec) {return kappa * kappa * kappa * kappa * rho(x) * vec[0]; }
	};
	const system_evaluation flexural_evaluation = { flexural, 1.0, { {0,0.0}, {1,0.0} }, { {2,0.0}, { 3,1.0 }} };

	// 1. ������ ������ �����
	// 1.1 ����������� ���� �����������
	const auto longitudinal_observed = longitidinal_evaluation.evaluate_the_right_part(min_kappa, max_kappa, points_k, 0);
	const auto flexural_observed = flexural_evaluation.evaluate_the_right_part(min_gamma, max_gamma, points_k, 0);
	std::vector<double> observed(longitudinal_observed);
	observed.insert(observed.end(), flexural_observed.begin(), flexural_observed.end());

	// 1.2 ��������� ���� �����������
	Parameters::kind = FIRST;
	e = [](double x) {return -x + 2; };//lagrange({ 0.26287, 0.88830 });
	rho = [](double x) {return -x + 2; };// lagrange({ 0.19171, 0.41566 });

	for (size_t i = 0; i < points_x; i++)
	{
		Parameters::piecewise_linear_params[i][0] = e(vv[i]);
		Parameters::piecewise_linear_params[i][1] = rho(vv[i]);
	}

	const std::vector<std::function<double(double, double, const std::vector<double>&)>> longitudinal_et = {
		[=](double x, double kappa, const std::vector<double>& vec) {return vec[1] / e(x); },
		[=](double x, double kappa, const std::vector<double>& vec) {return -kappa * kappa * rho(x) * vec[0]; }
	};
	const system_evaluation longitidinal_system_etalon = { longitudinal_et, 1.0, { {0,0.0} }, { { 1,1.0 }} };
	auto longitudinal_etalon = longitidinal_system_etalon.evaluate_the_right_part(min_kappa, max_kappa, points_k, 0);
	const std::vector<std::function<double(double, double, const std::vector<double>&)>> flexural_et = {
		[=](double x, double kappa, const std::vector<double>& vec) {return vec[1]; },
		[=](double x, double kappa, const std::vector<double>& vec) {return vec[2] / e(x); },
		[=](double x, double kappa, const std::vector<double>& vec) {return vec[3]; },
		[=](double x, double kappa, const std::vector<double>& vec) {return kappa * kappa * kappa * kappa * rho(x) * vec[0]; }
	};
	const system_evaluation flexural_evaluation_etalon = { flexural_et, 1.0, { {0,0.0}, {1,0.0} }, { {2,0.0}, { 3,1.0 }} };
	auto flexural_etalon = flexural_evaluation_etalon.evaluate_the_right_part(min_gamma, max_gamma, points_k, 0);
	std::vector<double> etalon(longitudinal_etalon);
	etalon.insert(etalon.end(), flexural_etalon.begin(), flexural_etalon.end());
	auto right_part = observed - etalon;



	// 2. ������ �������
	// 2.1 ������-��������� ��� �������� ����
	// ���������� ���������
	auto mu_reconstruct_longitudinal = [=](double kappa, const std::vector<double>& v, size_t idx)
	{return  -/*h_x */ v[1] * v[1] / Parameters::piecewise_linear_params[idx][0]
		/ Parameters::piecewise_linear_params[idx][0]; };

	auto rho_reconstruct_longitudinal = [=](double kappa, const std::vector<double>& v, size_t idx)
	{return  /*h_x  */ kappa * kappa * v[0] * v[0]; };

	// �������� ���������
	auto mu_reconstruct_flexural = [=](double kappa, const std::vector<double>& v, size_t idx)
	{return  /*h_x */ v[2] * v[2] / Parameters::piecewise_linear_params[idx][0]
		/ Parameters::piecewise_linear_params[idx][0]; };
	auto rho_reconstruct_flexural = [=](double kappa, const std::vector<double>& v, size_t idx)
	{return  -/*h_x */ kappa * kappa * kappa * kappa * v[0] * v[0]; };

	auto matrix = longitidinal_system_etalon.evaluate_the_matrix(points_kappa, vv, mu_reconstruct_longitudinal);
	auto flex = longitidinal_system_etalon.evaluate_the_matrix(points_kappa, vv, rho_reconstruct_longitudinal);
	for (size_t i = 0; i <= points_k; i++)
	{
		matrix[i].insert(matrix[i].end(), flex[i].begin(), flex[i].end());
	}
	auto matrix_flexural = flexural_evaluation_etalon.evaluate_the_matrix(points_gamma, vv, mu_reconstruct_flexural);
	flex = flexural_evaluation_etalon.evaluate_the_matrix(points_gamma, vv, rho_reconstruct_flexural);
	for (size_t i = 0; i <= points_k; i++)
	{
		matrix_flexural[i].insert(matrix_flexural[i].end(), flex[i].begin(), flex[i].end());
	}

	matrix.insert(matrix.end(), matrix_flexural.begin(), matrix_flexural.end());

	//lsqr _lsqr(matrix);
	//auto sol = _lsqr.lin_solve(right_part);

	VoyevodinMethod voyevodin_method(matrix, right_part, 1.0 / points_x, Dirichle, Dirichle, 1.0, 0.1, 1.0e-4, 0, 1.0e-8);
	auto sol = voyevodin_method.solution();

	Parameters::kind = THIRD;
	size_t iter = 0;
	std::vector<double> residials;
	auto norm_right_part = norm(right_part);
	residials.push_back(norm_right_part);
	while (norm_right_part > 0.1e-4)
	{
		for (size_t i = 0; i <= points_x; i++)
		{
			Parameters::piecewise_linear_params[i][0] += sol[i];
			Parameters::piecewise_linear_params[i][1] += sol[i + points_x];
		}
		if (iter > 100)
		{
			break;
		}
		longitudinal_etalon = longitidinal_evaluation.evaluate_the_right_part(min_kappa, max_kappa, points_k, 0);
		flexural_etalon = flexural_evaluation.evaluate_the_right_part(min_gamma, max_gamma, points_k, 0);
		etalon = longitudinal_etalon;
		etalon.insert(etalon.end(), flexural_etalon.begin(), flexural_etalon.end());
		right_part = observed - etalon;
		norm_right_part = norm(right_part);
		residials.push_back(norm_right_part);
		matrix = longitidinal_evaluation.evaluate_the_matrix(points_kappa, vv, mu_reconstruct_longitudinal);
		flex = longitidinal_evaluation.evaluate_the_matrix(points_kappa, vv, rho_reconstruct_longitudinal);
		for (size_t i = 0; i <= points_k; i++)
		{
			matrix[i].insert(matrix[i].end(), flex[i].begin(), flex[i].end());
		}
		matrix_flexural = flexural_evaluation.evaluate_the_matrix(points_gamma, vv, mu_reconstruct_flexural);
		flex = flexural_evaluation.evaluate_the_matrix(points_gamma, vv, rho_reconstruct_flexural);
		for (size_t i = 0; i <= points_k; i++)
		{
			matrix_flexural[i].insert(matrix_flexural[i].end(), flex[i].begin(), flex[i].end());
		}
		matrix.insert(matrix.end(), matrix_flexural.begin(), matrix_flexural.end());
		//_lsqr = { matrix };
		//sol = _lsqr.lin_solve(right_part);

		//VoyevodinMethod
		voyevodin_method = { matrix, right_part, 1.0 / points_x, Dirichle, Dirichle };
		sol = voyevodin_method.solution();
		++iter;
	}
	std::vector<double> mu_reconstrucred(points_x);
	std::vector<double> rho_reconstrucred(points_x);
	for (size_t i = 0; i < points_x; i++)
	{
		mu_reconstrucred[i] = Parameters::piecewise_linear_params[i][0];
		rho_reconstrucred[i] = Parameters::piecewise_linear_params[i][1];
	}
	plotTheWaveField({ {"black", exact_solution_mu}, {"red", mu_reconstrucred} }, "mu1.txt", h_x);
	plotTheWaveField({ {"black", exact_solution_rho}, {"red", rho_reconstrucred} }, "rho1.txt", h_x);
	plotTheWaveField({ {"black", residials} }, "res.txt", 1.0);
	//*/
	system("pause");
}
