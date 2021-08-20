#include "VoyevodinMethod.h"
#include "IterativeProcess.h"

/**
 * \brief ����������� ������ ���������
 * \param matrix ������� ����;
 * \param rightpart ������ ������ �����;
 * \param Step ����� ����������;
 * \param Left ������� �������;
 * \param Right ������� �������;
 * \param p �������� �������������;
 * \param alphaInitialValue ��������� �������� ��������� �������������;
 * \param H ����������� ���������;
 * \param Delta ����������� ������ �����;
 * \param eps ����������� ����������� ��������� �������������;
 */
VoyevodinMethod::VoyevodinMethod(const vector<vector<double>>& matrix,
	vector<double> rightpart, double Step, BoundaryCondition Left,
	BoundaryCondition Right, double p, double alphaInitialValue, double H,
	double Delta, double eps) : AlphaInitialValue(alphaInitialValue), step(Step),
	h(H), delta(Delta), RightPart(std::move(rightpart))
{
	//1. ������ ������� � �������� � � ����������������� ����
	// size = RightPart.size();
	//const auto rp = RightPart;
	matrix_system matrixSystem = { matrix, RightPart, step, p, Left, Right };
	//2. ��������� ������������ �������
	iterative_process iterativeProcess = { matrixSystem.diagonal(),
		matrixSystem.up_diagonal(), matrixSystem.rightPart(),
		matrixSystem.multiply_qtu(RightPart), AlphaInitialValue, step, h,
		delta, eps };
	//3. �������� �������
	Solution = iterativeProcess.solution();
	//4. ����������� � ����������� ������������
	matrixSystem.multiply_rtx(Solution);
	matrixSystem.multiply_Sinv(Solution);
}

vector<double> VoyevodinMethod::solution() const
{
	return Solution;
}
