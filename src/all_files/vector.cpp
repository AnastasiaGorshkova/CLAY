#include "vector.hpp"

double operator*(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0;
    for (std::size_t i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

std::vector<double> operator*(const std::vector<double>& a, const double& alpha) {
    std::vector<double> c;
    for (std::size_t i = 0; i < a.size(); i++) {
    	c.push_back(alpha * a[i]);
    }
    return c;
}

std::vector<double> operator*(const double& alpha, const std::vector<double>& a) {
    std::vector<double> c;
    for (std::size_t i = 0; i < a.size(); i++) {
    	c.push_back(a[i] * alpha );
    }
    return c;
}

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> sum;
    for (std::size_t i = 0; i < a.size(); i++) {
    	sum.push_back(a[i] + b[i]);
    }
    return sum;
}

std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> res;
    for (std::size_t i = 0; i < a.size(); i++) {
    	res.push_back(a[i] - b[i]);
    }
    return res;
}

void print(const std::vector<double>& result) {
    for (std::size_t i = 0; i < result.size(); ++i) {
        std::cout << result[i] << " ";
        }
    std::cout << std::endl;    
}
