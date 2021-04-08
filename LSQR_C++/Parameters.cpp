#include "Parameters.h"
#include <stdexcept>

double Parameters::get_value(double t, size_t param_num)
{
	if (param_num > const_params.size())
	{
		throw std::invalid_argument("");
	}
	size_t i = 0;
	for (; i < points.size(); i++)
	{
		if (t < points[i] || abs(t - points[points.size() - 1]) < DBL_EPSILON)
		{
			break;
		}
	}
	if (i >= points.size())
	{
		throw std::invalid_argument("");
	}
	if (!i)
	{
		return piecewise_linear_params[i][param_num];
	}
	const auto a = points[i - 1];
	const auto b = points[i];
	const auto fa = piecewise_linear_params[i - 1][param_num];
	const auto fb = piecewise_linear_params[i][param_num];
	return fa + (t - a) * (fb - fa) / (b - a);
}

double Parameters::evaluate(double t, size_t param_num)
{
	switch (kind) {
	case FIRST:
		return smooth_params[param_num](t);
	case SECOND:
		return 1;
	case THIRD:
		return get_value(t, param_num);
	default:
		throw std::invalid_argument("");
	}
}
