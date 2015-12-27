#include "Vector.h"


class GlobalMatrix
{
	public:
        std::vector<double> di; //������������ ��������
        std::vector<double> ggl; // �������� ������ ���������
        std::vector<double> ggu; // �������� ������� ���������
        std::vector<int> jg; // ������ � ��������� ������ �����(��������) � ������� ggl(ggu) 
        std::vector<int> ig; // ������ ��������(�����) ��������� � ������� ggl(ggu)
		Vector F; // ������ ������ �����
		std::vector<dot> list; //������ ������� ���������, �����, ������ ��� ����
        int n; // ����������� �������
        double normF; // ����� ������� F
		GlobalMatrix(void);
		~GlobalMatrix(void);
		void addElem(std::vector<int> el, std::vector<std::vector<double>> &G, std::vector<std::vector<double>> &M, std::vector<double> &F);
        void make(std::vector<pt> &point, std::vector<std::vector<int>> &finElem);
        void multMV(Vector &X, Vector &Y);
		void addToList(std::vector<int> &el);
		int findInList(dot x);
		void sortList();
		void addFirstBound(std::vector<int> &el,std::vector<pt> &point);
		void addSecondBound(std::vector<int> &el,std::vector<pt> &point);
		void addThirdBound(std::vector<int> &el,std::vector<pt> &point);
};

