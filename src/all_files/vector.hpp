#pragma once
#include <vector>
#include <iostream>

// Умножение векторв на вектор
double operator*(const std::vector<double>& a, const std::vector<double>& b);

// Умножение вектора на число
std::vector<double> operator*(const std::vector<double>& a, const double& alpha);

// Сложение векторов
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);

// Умножение числа на вектор
std::vector<double> operator*(const double& alpha, const std::vector<double>& a);

// Вычитание векторов
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);

void print(const std::vector<double>& result) ;
