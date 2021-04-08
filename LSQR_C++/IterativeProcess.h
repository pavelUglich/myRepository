#pragma once
#include <utility>
#include<vector>
#include<iterator>
#include "IterativeProcess.h"
using namespace std;
class IterativeProcess
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
	vector<double> Marching(const vector<double> & aa);
	double Residual(double alpha_);
	void IterationsRun();

public:
	IterativeProcess();
	IterativeProcess(
		const vector<double>& p1,
		vector<double> p2,
		vector<double> RightPart,
		vector<double> Qtu,
		double Alpha,
		double step,
		double h,
		double delta,
		double eps,
		int iterations = 50
	);


	~IterativeProcess();
	vector<double> solution() const {
		return Solution;
	}
};

