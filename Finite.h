#include "LOS.h"
class Finite
{
public:
	int h;
	LOS global;
	std::vector<std::vector<double>> G, M;
	std::vector<std::vector<int>> finElem, bound;
	std::vector<pt> point;
	std::vector<double> F;
	void makeGlobal();
	void makeGMF(int numElem);
	void recCond();
	void divElem();
	Finite(int x);
	~Finite(void);
};

