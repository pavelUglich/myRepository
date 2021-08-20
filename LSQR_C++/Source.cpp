#include <vector>
#include <cstdlib>
#include "OdeSolver.h"
#include "boundary_value_problem.h"
#include "system_evaluation.h"
#include "Parameters.h"
#include "vector_ops.h"
#include "plots.h"
#include <chrono>
#include <iomanip>

#include "lsqr.h"
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

std::vector<std::vector<double>> fill_up_the_matrix(
	const std::function<double(double, double)>& lambda,
	const std::vector<double>& points_x, const std::vector<double>& points_kappa)
{
	std::vector<std::vector<double>> result(points_kappa.size());
	for (size_t i = 0; i < points_kappa.size(); i++)
	{
		std::vector<double> row(points_x.size());
		for (size_t ii = 0; ii < points_x.size(); ii++)
		{
			row[ii] = lambda(points_kappa[i], points_x[ii]);
		}
		result[i] = row;
	}
	return result;
}

std::vector<double> fill_up_the_vector(
	const std::function<double(double)>& lambda, 
	const std::vector<double> & points)
{
	std::vector<double> result(points.size());
	for (size_t ii = 0; ii < points.size(); ii++)
	{
		result[ii] = lambda(points[ii]);
	}
	return result;
}


std::vector<double> create_points_vector(size_t points_k, double min_kappa, 
	double max_kappa)
{
	std::vector<double> points_kappa(points_k + 1);
	const auto h_kappa = (max_kappa - min_kappa) / points_k;
	for (size_t i = 0; i <= points_k; i++)
	{
		points_kappa[i] = min_kappa + i * h_kappa;
	}
	return points_kappa;
}

system_evaluation longitudinal(const std::function<double(double)>& e, const std::function<double(double)> & rho)
{
	const std::vector<std::function<double(double, double, const std::vector<double>&)>> longitudinal = {
		[=](double x, double kappa, const std::vector<double>& vec) {return vec[1] / e(x); },
		[=](double x, double kappa, const std::vector<double>& vec) {return -kappa * kappa * rho(x) * vec[0]; }
	};
	return  { longitudinal, { {0,0.0} }, { { 1,1.0 }} };
}

system_evaluation flexural(const std::function<double(double)>& e, const std::function<double(double)>& rho)
{
	const std::vector<std::function<double(double, double, const std::vector<double>&)>> flexural = {
			[=](double x, double kappa, const std::vector<double>& vec) {return vec[1]; },
			[=](double x, double kappa, const std::vector<double>& vec) {return vec[2] / e(x); },
			[=](double x, double kappa, const std::vector<double>& vec) {return vec[3]; },
			[=](double x, double kappa, const std::vector<double>& vec) {return kappa * kappa * kappa * kappa * rho(x) * vec[0]; }
	};
	return  { flexural, { {0,0.0}, {1,0.0} }, { {2,0.0}, { 3,1.0 }} };
}




int main()
{
	setlocale(0, "");
	Parameters::kind = FIRST;

	const size_t points_k = 20;
	const double min_kappa = 2.4;
	const double max_kappa = 4.2;


	// минимальная частота
	const double min_gamma = 2.3;
	const double max_gamma = 4.1;
	// количество точек
	const size_t points_x = 20;
	// значения частоты
	const double h_x = 1.0 / points_x;
	const double h_kappa = (max_kappa - min_kappa) / points_k;
	const std::vector<double> points_kappa = create_points_vector(points_k, min_kappa, max_kappa);
	const std::vector<double> points_gamma = create_points_vector(points_k, min_gamma, max_gamma);
	const std::vector<double> vv = create_points_vector(points_x, 0, 1);
	
	// задаём параметры
	Parameters::points = vv;
	Parameters::piecewise_linear_params = { points_x + 1, {1.0, 1.0} };
	Parameters::const_params = { 1.0, 1.0 };
	Parameters::smooth_params = {
		// 1
		[](auto x) {return 1.0 + 0.1 * x;/* 2.0 - (x - 1) * (x - 1);*/ },
		[](auto x) {return 1.0 + 0.1 * x; /*1.0 + exp(5 * (x - 1));*/ }
		// 2
		//[](auto x) {return 2.0 - x * x; },
		//[](auto x) {return 1.0 + exp(-4 * x + 0.02); }
		// 3
		//[](auto x) {return 1.0 + sin(pi * x); },
		//[](auto x) {return 2 * x * x - 2 * x + 1; }

	
	};
	
	/* 20.08
	std::vector<double> exact_solution_mu = fill_up_the_vector(Parameters::smooth_params[0], vv);//(points_x);
	std::vector<double> exact_solution_rho = fill_up_the_vector(Parameters::smooth_params[1], vv);


	std::vector<double> exact_solution_mu1 = fill_up_the_vector([](auto x) {return 0.1 * x;/* 2.0 - (x - 1) * (x - 1);*/ //}, vv);//(points_x);
	//std::vector<double> exact_solution_rho1 = fill_up_the_vector([](auto x) {return 0.1 * x;/* 2.0 - (x - 1) * (x - 1);*/ }, vv);

	/*
	std::vector<double> exact_solution_mu = fill_up_the_vector([](double x) {return 1.0 + x / 10; }, vv);//(points_x);
	std::vector<double> exact_solution_rho = fill_up_the_vector([](double x) {return 1.0 + x * x / 50; }, vv);


	std::vector<double> exact_solution_mu1 = fill_up_the_vector([](auto x) {return 1.0 / 50 / x / x / x * (51 * x * x + 5 * x - 2) * sin(x) + (-55 * x * x + 2 * x) * cos(x) + 50 * x * x); , vv);
	std::vector<double> exact_solution_rho1 = fill_up_the_vector([](auto x) {return 1.0 / 600 / x / x / x * (24 * x * x + 3) * cos(x) * sin(x) + 9 * x * cos(x) * cos(x) + 617 * x * x * x - 12*x); }, vv);
	*/

	/*
	// краевые задачи
	const system_evaluation longitidinal_evaluation = longitudinal(
		[](double x) { return Parameters::evaluate(x, 0); },
		[](double x) { return Parameters::evaluate(x, 1); });
	const system_evaluation flexural_evaluation = flexural(
		[](double x) { return Parameters::evaluate(x, 0); },
		[](double x) { return Parameters::evaluate(x, 1); });
	const auto longitudinal_observed = longitidinal_evaluation.evaluate_the_right_part(min_kappa, max_kappa, points_k, 0);
	const auto flexural_observed = flexural_evaluation.evaluate_the_right_part(min_gamma, max_gamma, points_k, 0);
	std::vector<double> observed(longitudinal_observed);
	observed.insert(observed.end(), flexural_observed.begin(), flexural_observed.end());
	*/

	// 1.2 Эталонное поле перемещений
	/* 20.08
	Parameters::kind = FIRST;
	auto e = [](double x) {return 1;/* 2 - x*/ //};
	//auto rho = [](double x) {return 1; /* 2 - x */ //};
	/*for (size_t i = 0; i < points_x; i++)
	{
		Parameters::piecewise_linear_params[i][0] = e(vv[i]);
		Parameters::piecewise_linear_params[i][1] = rho(vv[i]);
	}
	const system_evaluation longitidinal_system_etalon = longitudinal(e, rho);
	auto longitudinal_etalon = longitidinal_system_etalon.evaluate_the_right_part(min_kappa, max_kappa, points_k, 0);
	const system_evaluation flexural_evaluation_etalon = flexural(e, rho);
	auto flexural_etalon = flexural_evaluation_etalon.evaluate_the_right_part(min_gamma, max_gamma, points_k, 0);
	std::vector<double> etalon(longitudinal_etalon);
	etalon.insert(etalon.end(), flexural_etalon.begin(), flexural_etalon.end());
	*/
	auto longitudinal = fill_up_the_vector([](auto x) {return 1.0 / 50 / x / x / x * ((51 * x * x + 5 * x - 2) * sin(x) + (-55 * x * x + 2 * x) * cos(x) + 50 * x * x); }, points_kappa);
	auto flexural_etalon = fill_up_the_vector([](auto x) {return 1.0 / 600 / x / x / x * ((24 * x * x + 3) * cos(x) * sin(x) + 9 * x * cos(x) * cos(x) + 617 * x * x * x - 12 * x); }, points_gamma);
	std::vector<double> right_part(longitudinal);
	right_part.insert(right_part.end(), flexural_etalon.begin(), flexural_etalon.end());


	//auto right_part = observed - etalon;


	// 2. Создаём матрицу
	// 2.1 Лямбда-выражения для создания ядер
	// продольные колебания
/*
	auto mu_reconstruct_longitudinal = [=](double kappa, const std::vector<double>& v, size_t idx)
	{return  - v[1] * v[1] / Parameters::piecewise_linear_params[idx][0]
		/ Parameters::piecewise_linear_params[idx][0]; };

	auto rho_reconstruct_longitudinal = [=](double kappa, const std::vector<double>& v, size_t idx)
	{return   kappa * kappa * v[0] * v[0]; };

	auto mu_reconstruct_flexural = [=](double kappa, const std::vector<double>& v, size_t idx)
	{return  v[2] * v[2] / Parameters::piecewise_linear_params[idx][0]
		/ Parameters::piecewise_linear_params[idx][0]; };
	auto rho_reconstruct_flexural = [=](double kappa, const std::vector<double>& v, size_t idx)
	{return  - kappa * kappa * kappa * kappa * v[0] * v[0]; };

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
	*/

	auto matrix = fill_up_the_matrix([](double x, double xi) {return sin(x * xi); }, vv, points_kappa);
	auto flex = fill_up_the_matrix([](double x, double xi) {return cos(x * xi); }, vv, points_kappa);
	for (size_t i = 0; i <= points_k; i++)
	{
		matrix[i].insert(matrix[i].end(), flex[i].begin(), flex[i].end());
	}
	auto matrix_flexural = fill_up_the_matrix([](double x, double xi) {return cos(x * xi)* cos(x * xi); }, vv, points_gamma);
	flex = fill_up_the_matrix([](double x, double xi) {return sin(x * xi) * sin(x * xi); }, vv, points_gamma);
	for (size_t i = 0; i <= points_k; i++)
	{
		matrix_flexural[i].insert(matrix_flexural[i].end(), flex[i].begin(), flex[i].end());
	}


	matrix.insert(matrix.end(), matrix_flexural.begin(), matrix_flexural.end());


	//matrix = matrix * h_x;
	//lsqr _lsqr(matrix);
	//auto sol = _lsqr.lin_solve(right_part);

/*	std::vector<double> exact(exact_solution_mu1);
	exact.insert(exact.end(), exact_solution_rho1.begin(), exact_solution_rho1.end());
	const auto right_part1 = matrix * exact;
	*/
	VoyevodinMethod voyevodin_method(matrix, right_part, h_x, Dirichle, Neumann, 1.0, 0.1, 1.0e-4, 0, 1.0e-3);
	auto sol = voyevodin_method.solution();

/*	Parameters::kind = THIRD;
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
		voyevodin_method = { matrix, right_part, h_x, Dirichle, Neumann, 1.0, 0.1, 1.0e-4, 0, 1.0e-3 };
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
	plotTheWaveField({ {"black", exact_solution_mu}, {"red", mu_reconstrucred} }, "mu513_.txt", h_x);
	plotTheWaveField({ {"black", exact_solution_rho}, {"red", rho_reconstrucred} }, "rho513_.txt", h_x);
	plotTheWaveField({ {"black", residials} }, "res.txt", 1.0);*/
	
	system("pause");
}
