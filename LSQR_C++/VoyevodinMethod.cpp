#include "VoyevodinMethod.h"
#include "IterativeProcess.h"

/**
 * \brief Конструктор метода Воеводина
 * \param matrix Матрица СЛАУ;
 * \param rightpart вектор правой части;
 * \param Step Длина подотрезка;
 * \param Left Краевое условие;
 * \param Right Краевое условие;
 * \param p Параметр стабилизатора;
 * \param alphaInitialValue Начальное значение параметра регуляризации;
 * \param H Погрешность оператора;
 * \param Delta Погрешность правой части;
 * \param eps Погрешность определения параметра регуляризации;
 */
VoyevodinMethod::VoyevodinMethod(const vector<vector<double>>& matrix,
	vector<double> rightpart, double Step, BoundaryCondition Left,
	BoundaryCondition Right, double p, double alphaInitialValue, double H,
	double Delta, double eps) : AlphaInitialValue(alphaInitialValue), step(Step),
	h(H), delta(Delta), RightPart(std::move(rightpart))
{
	//1. Создаём систему и приводим её к двухдиагональному виду
	// size = RightPart.size();
	//const auto rp = RightPart;
	matrix_system matrixSystem = { matrix, RightPart, step, p, Left, Right };
	//2. Запускаем итерационный процесс
	iterative_process iterativeProcess = { matrixSystem.diagonal(),
		matrixSystem.up_diagonal(), matrixSystem.rightPart(),
		matrixSystem.multiply_qtu(RightPart), AlphaInitialValue, step, h,
		delta, eps };
	//3. Получаем решение
	Solution = iterativeProcess.solution();
	//4. Возвращаемя к изначальным неизвестнымю
	matrixSystem.multiply_rtx(Solution);
	matrixSystem.multiply_Sinv(Solution);
}

vector<double> VoyevodinMethod::solution() const
{
	return Solution;
}
