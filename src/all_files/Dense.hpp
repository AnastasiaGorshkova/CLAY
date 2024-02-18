#pragma once

#include <iostream>
#include <vector>

class Dense 
{
	private:
	    int rows, cols;
	    std::vector<double> data;
	public:
	    // конструктор класса 
	    Dense(int rows, int cols);

	    // доступ к элементу матрицы по индексам для изменения
	    double& operator()(int i, int j);
	    
	    // доступ к элементу матрицы по индексам для чтения
	    const double& operator()(int i, int j) const;
	    
	    // сложение двух матриц
	    Dense operator+(const Dense& other) const;
	    
	    // умножение матрицы на число
	    Dense operator*(double alpha) const;
	    
	    // умножение матрицы на вектор
	    std::vector<double> operator*(const std::vector<double>& vec) const;
	    
	    // умножение двух матриц
	    Dense operator*(const Dense& other) const;
	    
	    // вывод матрицы
	    void print() const;
};
