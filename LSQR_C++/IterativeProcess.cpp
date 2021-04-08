#include "IterativeProcess.h"
#include<cmath>
#include "MatrixSystem.h"


/**
 * \brief Метод прогоник
 * \param aa главная диагональ трёхдиагональной матрицв
 * \return решение СЛАУ
 */
vector<double> IterativeProcess::Marching(const vector<double> &aa) {
	vector<double> x(size);
	vector<double> xi(size + 1);
	vector<double> eta(size + 1);
	xi[0] = 0;
	eta[0] = 0;
	for (size_t i = 0; i < size; i++) {
		xi[i + 1] = b[i] / (aa[i] - c[i] * xi[i]);
		eta[i + 1] = (c[i] * eta[i] - RightPart[i]) / (aa[i] - c[i] * xi[i]);
	}
	x[size - 1] = eta[size];
	for (size_t i = 1; i < size; i++)
		x[size - i - 1] = xi[size - i] * x[size - i] + eta[size - i];
	return x;
}

/**
 * \brief Обобщенная невязка 
 * \param alpha_ значение параметра регуляризации
 * \return значение невязки
 */
double IterativeProcess::Residual(/*vector<double> aa,*/ double alpha_) {
	vector<double> aa(this->a.size());
	for (size_t i = 0; i < size; i++)
		aa[i] = -a[i] - alpha_;
	Solution = Marching(aa);//
	double nz = norm(Solution);
	vector<double> pz(size);
	for (size_t i = 0; i < size - 1; i++) pz[i] = p1[i] * Solution[i] + p2[i + 1] * Solution[i + 1];
	pz[size - 1] = p1[size - 1] * Solution[size - 1];
	for (size_t i = 0; i < size; i++) pz[i] -= Qtu[i];
	double npz = norm(pz);
	return step*step*npz*npz - (delta + h*nz)*(delta + h*nz);
}

/**
 * \brief Итеративный процесс для подбора оптимального 
 *        значения параметра регуляризации
 */
void IterativeProcess::IterationsRun(){
	double AlphaS = alpha;
	double AlphaN = alpha * 0.5;
	double ss = Residual(AlphaS);
	double sn = Residual(AlphaN);
	for (int i = 0; i < Iterations; i++) {
		const double alpha_ = AlphaN / (1 - (1 / AlphaS)*(AlphaS - AlphaN)*sn / (sn - ss));
		ss = sn;
		sn = Residual(alpha_);
		AlphaS = AlphaN;
		AlphaN = alpha_;
		if (abs(AlphaS - AlphaN) < eps) break;
	}
}

/**
 * \brief построение трёхдиагональной матрицы
 */
void IterativeProcess::tridiag(){
	a.resize(size);
	b.resize(size);
	c.resize(size);
	a[0] = p1[0] * p1[0];
	a[size - 1] = p1[size - 1] * p1[size - 1] + p2[size - 1] * p2[size - 1];
	b[0] = p1[0] * p2[0];
	b[size - 1] = 0;
	for (size_t i = 1; i < size - 1; i++) {
		a[i] = p2[i - 1] * p2[i - 1] + p1[i] * p1[i];
		b[i] = p1[i] * p2[i];
	}
	c[0] = 0;
	for (size_t i = 1; i < size; i++) c[i] = b[i - 1];
}


IterativeProcess::IterativeProcess() = default;


IterativeProcess::IterativeProcess(const vector<double>& p1, vector<double> p2, 
	vector<double> RightPart, vector<double> Qtu, double Alpha, double step, 
	double h, double delta, double eps, int iterations): p1(p1), 
	p2(std::move(p2)), Qtu(std::move(Qtu)), RightPart(std::move(RightPart)), 
	alpha(Alpha), Iterations(iterations), step(step), h(h), delta(delta), 
	eps(eps)
{
	size = p1.size();
	tridiag(); 
	IterationsRun();
}


IterativeProcess::~IterativeProcess() = default;
