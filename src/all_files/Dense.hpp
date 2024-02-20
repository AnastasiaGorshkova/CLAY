#pragma once

#include <iostream>
#include <vector>
#include "vector.hpp"

class Dense 
{
	private:
	    std::size_t rows, cols;
	    std::vector<double> data;
	public:
	    // конструктор класса 
	    Dense(const std::vector<double>& data, std::size_t rows, std::size_t cols);
	    
	    // доступ к элементу матрицы по индексам для чтения
	    double operator()(std::size_t i, std::size_t j) const;
	    
	    // количество столбцов
	    int get_cols() const;
	    
	    // количество строк
	    int get_rows() const;
	    
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
