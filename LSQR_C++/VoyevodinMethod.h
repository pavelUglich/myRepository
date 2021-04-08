#pragma once
#include <utility>
#include<vector>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
//#include "Stabilizer.h"
using namespace std;

class VoyevodinMethod
{
	double AlphaInitialValue;//Начальное значение параметра регуляризации
	double step; //длина отрезка разбиения
	//double eps;//Погрешность определения параметра регуляризации
	double h; //Погрешность оператора
	double delta;//Погрешность правой части
	vector<double> RightPart, p1, p2, Qtu;
	size_t size;
	vector<double> Solution;

public:
	VoyevodinMethod();

	///<summary>
	///Конструктор метода Воеводина
	///</summary>
	///<param name="matrix">Матрица СЛАУ;</param>
	///<param name="rightpart">вектор правой части;</param>
	///<param name="Step">Длина подотрезка;</param>
	///<param name="Left">Краевое условие;</param>
	///<param name="Right">Краевое условие;</param>
	///<param name="p">Параметр стабилизатора;</param>
	///<param name="alphaInitialValue">Начальное значение параметра регуляризации;</param>
	///<param name="H">Погрешность оператора;</param>
	///<param name="Delta">Погрешность правой части;</param>
	///<param name="eps">Погрешность определения параметра регуляризации;</param>
	VoyevodinMethod(const vector<vector<double>> & matrix,
        vector<double> rightpart,
		double Step,
		BoundaryCondition Left = Neumann,
		BoundaryCondition Right = Neumann,
		double p = 1.0,
		double alphaInitialValue = 0.1e-1,
		double H = 1.0e-4,
		double Delta = 1.0e-4,
		double eps = 0.1e-3) :
		AlphaInitialValue(alphaInitialValue), step(Step), /*eps(eps),*/ h(H), delta(Delta), RightPart(
			std::move(rightpart)) {
		//1. Создаём систему и приводим её к двухдиагональному виду
		size = RightPart.size();
		MatrixSystem * matrixSystem = new MatrixSystem(matrix, RightPart, step, p, Left, Right);
		//2. Запускаем итерационный процесс
		IterativeProcess *iterativeProcess = new IterativeProcess(
			matrixSystem->Diagonal(),
			matrixSystem->UpDiagonal(),
			matrixSystem->rightPart(),
			matrixSystem->MultiplyQtu(RightPart),
			AlphaInitialValue,
			step,
			h,
			delta,
			eps);
		//3. Получаем решение
		Solution = iterativeProcess->solution();
		delete iterativeProcess;
		//4. Возвращаемя к изначальным неизвестнымю
		matrixSystem->multiply_Rtx(Solution);
		matrixSystem->multiply_Sinv(Solution);
		delete matrixSystem;
	};


	~VoyevodinMethod();
	vector<double> solution() const {
		return Solution;
	};

};

