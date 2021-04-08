#pragma once
#include <utility>
#include<vector>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
//#include "Stabilizer.h"
using namespace std;

class VoyevodinMethod
{
	double AlphaInitialValue;//��������� �������� ��������� �������������
	double step; //����� ������� ���������
	//double eps;//����������� ����������� ��������� �������������
	double h; //����������� ���������
	double delta;//����������� ������ �����
	vector<double> RightPart, p1, p2, Qtu;
	size_t size;
	vector<double> Solution;

public:
	VoyevodinMethod();

	///<summary>
	///����������� ������ ���������
	///</summary>
	///<param name="matrix">������� ����;</param>
	///<param name="rightpart">������ ������ �����;</param>
	///<param name="Step">����� ����������;</param>
	///<param name="Left">������� �������;</param>
	///<param name="Right">������� �������;</param>
	///<param name="p">�������� �������������;</param>
	///<param name="alphaInitialValue">��������� �������� ��������� �������������;</param>
	///<param name="H">����������� ���������;</param>
	///<param name="Delta">����������� ������ �����;</param>
	///<param name="eps">����������� ����������� ��������� �������������;</param>
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
		//1. ������ ������� � �������� � � ����������������� ����
		size = RightPart.size();
		MatrixSystem * matrixSystem = new MatrixSystem(matrix, RightPart, step, p, Left, Right);
		//2. ��������� ������������ �������
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
		//3. �������� �������
		Solution = iterativeProcess->solution();
		delete iterativeProcess;
		//4. ����������� � ����������� ������������
		matrixSystem->multiply_Rtx(Solution);
		matrixSystem->multiply_Sinv(Solution);
		delete matrixSystem;
	};


	~VoyevodinMethod();
	vector<double> solution() const {
		return Solution;
	};

};

