#include "IterativeProcess.h"
#include<cmath>
#include "MatrixSystem.h"


/**
 * \brief ћетод прогонка
 * \param aa главна€ диагональ трЄхдиагональной матрицв
 * \return решение —Ћј”
 */
vector<double> iterative_process::marching(const vector<double> &aa) {
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
 * \brief ќбобщенна€ нев€зка 
 * \param alpha_ значение параметра регул€ризации
 * \return значение нев€зки
 */
double iterative_process::residual(double alpha_) {
	vector<double> aa(this->a.size());
	for (size_t i = 0; i < size; i++)
		aa[i] = -a[i] - alpha_;
	Solution = marching(aa);//
	const double nz = norm(Solution);
	vector<double> pz(p1.size());
	vector<double> pp(p2);
	p2.push_back(0);
	for (size_t i = 0; i < p1.size() - 1; i++) pz[i] = p1[i] * Solution[i] + p2[i + 1] * Solution[i + 1];
	p2.erase(p2.end() - 1);
	pz[p1.size() - 1] = p1[p1.size() - 1] * Solution[p1.size() - 1];
	for (size_t i = 0; i < pz.size(); i++) pz[i] -= Qtu[i];
	const double npz = norm(pz);
	return step*step*npz*npz - (delta + h*nz)*(delta + h*nz);
}

/**
 * \brief »теративный процесс дл€ подбора оптимального 
 *        значени€ параметра регул€ризации
 */
void iterative_process::iterations_run(){
	double alpha_s = alpha;
	double alpha_n = alpha * 0.5;
	double ss = residual(alpha_s);
	double sn = residual(alpha_n);
	for (int i = 0; i < Iterations; i++) {
		const double alpha_ = alpha_n / (1 - 1 / alpha_s * (alpha_s - alpha_n) * sn / (sn - ss));
		ss = sn;
		sn = residual(alpha_);
		alpha_s = alpha_n;
		alpha_n = alpha_;
		if (abs(alpha_s - alpha_n) < eps) break;
	}
}

/**
 * \brief построение трЄхдиагональной матрицы
 */
void iterative_process::tridiag(){
	size = this->RightPart.size();
	a.resize(size);
	b.resize(size);
	c.resize(size);
	const auto size_ = p1.size();
	a[0] = p1[0] * p1[0];
	a[size_ - 1] = p1[size_ - 1] * p1[size_ - 1]; //  +p2[size_ - 1] * p2[size_ - 1];
	if (p2.size() == p1.size())
	{
		a[size_ - 1] += p2[size_ - 1] * p2[size_ - 1];
	}
	b[0] = p1[0] * p2[0];
	b[size_ - 1] = 0;
	for (size_t i = 1; i < size_ - 1; i++) {
		a[i] = p2[i - 1] * p2[i - 1] + p1[i] * p1[i];
		b[i] = p1[i] * p2[i];
	}
	c[0] = 0;
	for (size_t i = 1; i < size_; i++) c[i] = b[i - 1];
}

/**
 * \brief конструктор
 * \param p1 диагональ двухдиагональной матрицы
 * \param p2 наддиагональ двухдиагональной матрицы
 * \param right_part права€ часть
 * \param qtu преобразованна€ права€ часть
 * \param alpha параметр регул€ризации
 * \param step шаг
 * \param h погрешность оператора 
 * \param delta погрешность правой части
 * \param eps точность определени€ параметра регул€ризации
 * \param iterations предельное количество итераций
 */
iterative_process::iterative_process(const vector<double>& p1, vector<double> p2, 
	vector<double> right_part, vector<double> qtu, double alpha, double step, 
	double h, double delta, double eps, int iterations): p1(p1), 
	p2(std::move(p2)), Qtu(std::move(qtu)), RightPart(std::move(right_part)), 
	alpha(alpha), Iterations(iterations), step(step), h(h), delta(delta), 
	eps(eps)
{
	size = right_part.size();
	tridiag(); 
	iterations_run();
}



vector<double> iterative_process::solution() const
{
	return Solution;
}
