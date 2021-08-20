#pragma once
#include<vector>
#include "IterativeProcess.h"
using namespace std;
class iterative_process
{
	vector<double> a;
	vector<double> b;
	vector<double> c;
	vector<double> p1;
	vector<double> p2;
	vector<double> Qtu;
	vector<double> RightPart;
	vector<double> Solution;
	double alpha{};
	int Iterations{};
	size_t size{};
	double step{}, h{}, delta{}, eps{};

	void tridiag();
	vector<double> marching(const vector<double> & aa);
	double residual(double alpha_);
	void iterations_run();

public:
	iterative_process(
		const vector<double>& p1,
		vector<double> p2,
		vector<double> right_part,
		vector<double> qtu,
		double alpha,
		double step,
		double h,
		double delta,
		double eps,
		int iterations = 50
	);

	vector<double> solution() const;
};

