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

std::vector<double> operator/(const std::vector<double>& a, const double& alpha) {
    std::vector<double> c(a.size());
    for (std::size_t i = 0; i < a.size(); i++) {
    	c[i] = a[i] / alpha ;
    }
    return c;
}



double modul(const std::vector<double>& x){
    double result = 0;
    for(std::size_t i = 0; i < x.size(); i++){
        result+= std::pow(x[i],2);
    }
    return std::pow(result,0.5);
}

double sign(const std::vector<double>& x){
    if(x[0] < 0) return -1;
    return 1;
}

std::vector<double> transposed(const std::vector<double>& x) {
    std::vector<double> result;
    for (std::size_t i = 0; i < x.size(); i++) {
        for (std::size_t j = 0; j < x.size(); j++) {
            result.push_back(x[i] * x[j]);
        }
    }
    return result;
}
