#include "Finite.h"



// считывание и выделелние памяти
Finite::Finite(int x)
{
	int n = 0;
	h = x;
	std::ifstream coord("coord.txt");
	std::ifstream elem("elem.txt");
	std::ifstream condition("cond.txt");
	coord >> n;
	point.resize(n);
	for (int i = 0; i < n; i++)
	{
		coord >> point[i].x >> point[i].y;
	}
	elem >> n;
	finElem.resize(n);
	for (int i = 0; i < n; i++)
	{
		finElem[i].resize(4);
		for (int j = 0; j < 4; j++)
		{
			elem >> finElem[i][j];
		}
	}
	condition >> n;
	bound.resize(n);
	for (int i = 0; i < n; i++)
	{
		bound[i].resize(4);
		for (int j = 0; j < 4; j++)
		{
			condition >> bound[i][j];
		}
	}
	global.make(point, finElem);
	M.resize(4);
	G.resize(4);
	F.resize(4);
	for (int i = 0; i < 4; i++)
	{
		M[i].resize(4);
		G[i].resize(4);
	}
}

// сборка глобальной матрицы
void Finite::makeGlobal()
{
	for (int i = 0; i < finElem.size(); i++)
	{
		makeGMF(i);
		global.add(finElem[i], G, M, F);
	}
}

// вычисление матриц масс и жесткости и вектора для элемента
void Finite::makeGMF(int numElem)
{
	double hx, hy;
	std::vector<double> gam(4), f(4);
	pt max, min;
	// вот за это внатуре стыдно, но все для пользователей
	max = point[finElem[numElem][0]];
	min = point[finElem[numElem][0]];

	for (int i = 1; i < 4; i++)
	{
		if (point[finElem[numElem][i]].x > max.x && point[finElem[numElem][i]].y > max.y)
		{
			max = point[finElem[numElem][i]];
		}
		if (point[finElem[numElem][i]].x < min.x && point[finElem[numElem][i]].y < min.y)
		{
			min = point[finElem[numElem][i]];
		}
	}
	hx = max.x - min.x;
	hy = max.y - min.y;
	for (int i = 0; i < 4; i++)
	{
		gam[i] = gamma(point[finElem[numElem][i]]);
		f[i] = func(point[finElem[numElem][i]]);
	}
	// строрим матрицу жесткости
	G[0][0] = 1.0/6.0 * (2.0 * hy/hx + 2.0 * hx/hy); 
	G[0][1] = 1.0/6.0 * (-2.0 * hy/hx + hx/hy);
	G[0][2] = 1.0/6.0 * (hy/hx - 2.0 * hx/hy);
	G[0][3] = 1.0/6.0 * (-hy/hx - hx/hy);

	G[1][0] = G[0][1];
	G[1][1] = G[0][0];
	G[1][2] = G[0][3];
	G[1][3] = G[0][2];

	G[2][0] = G[0][1];
	G[2][1] = G[0][3];
	G[2][2] = G[0][0];
	G[2][3] = G[0][1];

	G[3][0] = G[0][3];
	G[3][1] = G[0][2];
	G[3][2] = G[0][1];
	G[3][3] = G[0][0];
	// строим матрицу массы
	M[0][0] = hx * hy * (gam[0]/16.0 * gam[1]/48.0 * gam[2]/48.0 * gam[3]/144.0); 
	M[0][1] = hx * hy * (gam[0]/48.0 * gam[1]/48.0 * gam[2]/144.0 * gam[3]/144.0);
	M[0][2] = hx * hy * (gam[0]/48.0 * gam[1]/144.0 * gam[2]/48.0 * gam[3]/144.0);
	M[0][3] = hx * hy * (gam[0]/144.0 * gam[1]/144.0 * gam[2]/144.0 * gam[3]/144.0);

	M[1][0] = hx * hy * (gam[0]/48.0 * gam[1]/48.0 * gam[2]/144.0 * gam[3]/144.0);
	M[1][1] = hx * hy * (gam[0]/48.0 * gam[1]/16.0 * gam[2]/144.0 * gam[3]/48.0);
	M[1][2] = hx * hy * (gam[0]/144.0 * gam[1]/144.0 * gam[2]/144.0 * gam[3]/144.0);
	M[1][3] = hx * hy * (gam[0]/144.0 * gam[1]/48.0 * gam[2]/144.0 * gam[3]/48.0);

	M[2][0] = hx * hy * (gam[0]/48.0 * gam[1]/144.0 * gam[2]/48.0 * gam[3]/144.0);
	M[2][1] = hx * hy * (gam[0]/144.0 * gam[1]/144.0 * gam[2]/144.0 * gam[3]/144.0);
	M[2][2] = hx * hy * (gam[0]/48.0 * gam[1]/144.0 * gam[2]/16.0 * gam[3]/48.0);
	M[2][3] = hx * hy * (gam[0]/144.0 * gam[1]/144.0 * gam[2]/48.0 * gam[3]/48.0);

	M[3][0] = hx * hy * (gam[0]/144.0 * gam[1]/144.0 * gam[2]/144.0 * gam[3]/144.0);
	M[3][1] = hx * hy * (gam[0]/144.0 * gam[1]/48.0 * gam[2]/144.0 * gam[3]/48.0);
	M[3][2] = hx * hy * (gam[0]/144.0 * gam[1]/144.0 * gam[2]/48.0 * gam[3]/48.0);
	M[3][3] = hx * hy * (gam[0]/144.0 * gam[1]/48.0 * gam[2]/48.0 * gam[3]/16.0);

	F[0] = f[0]/36 * (hx * hy * 4) + f[1]/36 * (hx * hy * 2) + f[2]/36 * (hx * hy * 2) + f[3]/36 * (hx * hy);
	F[1] = f[0]/36 * (hx * hy * 2) + f[1]/36 * (hx * hy * 4) + f[2]/36 * (hx * hy) + f[3]/36 * (hx * hy * 2);
	F[2] = f[0]/36 * (hx * hy * 2) + f[1]/36 * (hx * hy) + f[2]/36 * (hx * hy * 4) + f[3]/36 * (hx * hy * 2);
	F[3] = f[0]/36 * (hx * hy) + f[1]/36 * (hx * hy * 2) + f[2]/36 * (hx * hy * 2) + f[3]/36 * (hx * hy * 4);


}

// чушь
void Finite::divElem()
{
	int ind = 0, i = 0;
	pt buf;
	std::vector<pt> newpoint;
	newpoint.reserve(point.size());
	h /= 2;
	while(h != 0)
	{
		for (ind = 0; ind < point.size(); )
		{
			while(point[i].y == point[i+1].y)
			{
				newpoint.push_back(point[i]);
				buf.x = point[i + 1].x / 2.0;
				buf.y = point[i + 1].y;
				newpoint.push_back(buf);
				i++;
			}
		}
	}
}

void Finite::recCond()
{
	for (int i = 0; i < bound.size(); i++)
	{
		if(bound[i][2] == 2)
		{
			global.secondBound(bound[i], point);
		}
		if(bound[i][2] == 3)
		{
			global.thirdBound(bound[i], point);
		}
	}
	for (int i = 0; i < bound.size(); i++)
	{
		if(bound[i][2] == 1)
		{
			global.firstBound(bound[i], point);
		}
	}
}

Finite::~Finite(void)
{
}
