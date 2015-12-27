#include "Vector.h"


Vector::Vector(void)
{
}

//���� �������� � ��� ��������� ������� ����������� � ���������� ������
Vector::Vector(int n)
{
        V.resize(n);
}

Vector::~Vector(void)
{
}

void Vector::make(int n)
{
        V.resize(n);
}


void Vector::read(std::ifstream &vect)
{
        for (int i = 0; i < V.size(); i++)
        {
                vect >> V[i];
        }
}

//���������� ��������� ������������
Vector& Vector::operator=(const Vector& newVect)
        {
                if (this != &newVect)
                        this->V = newVect.V;
                return *this;
        }

//���������� ������������ �������� ��� ��������� ������������
double Vector::operator* (const Vector& newVect)
{
        double sum = 0;
        if (this->V.size() != newVect.V.size())
                return sum;
        else
        {
                for(int i = 0; i < this->V.size(); i++)
                        sum += this->V[i] * newVect.V[i];
        }
        return sum;
}

// ����� �������. � �� ����� ����� ������ �� ������������ ����, �� � ���� ������ �����������
double Vector::norm()
{
        double norma = 0;
        for(int i = 0; i < this->V.size(); i++)
        {
			norma += pow(this->V[i],2.0);
        }
        norma = sqrt(norma);
        return norma;
}

//���������� �������� ��������. ��������� �������� � ������ �������
Vector& Vector::operator-(Vector& newVect)
{
        if (this->V.size() != newVect.V.size())
        {
                return *this;
        }
        else
        {
                for (int i = 0; i < this->V.size(); i++)
                {
                        newVect.V[i] = this->V[i] - newVect.V[i]; 
                }
                return newVect;
        }
}

//���������� ��������� ��������, ��������� �������� � ����� �������, � ���� � ��� ���� ����� �� ���������
Vector& Vector::operator+(Vector& newVect)
{
        if (this->V.size() != newVect.V.size())
        {
                return *this;
        }
        else
        {
                for (int i = 0; i < this->V.size(); i++)
                {
                        this->V[i] = this->V[i] + newVect.V[i]; 
                }
                return *this;
        }

}