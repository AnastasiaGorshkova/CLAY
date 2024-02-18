#include "Dense.hpp"


Dense::Dense(int rows, int cols) : rows(rows), cols(cols), data(rows*cols) {}

double& Dense::operator()(int i, int j) {
    return data[i*cols + j];
}

const double& Dense::operator()(int i, int j) const {
    return data[i*cols + j];
}

Dense Dense::operator+(const Dense& other) const {
    Dense result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result(i, j) = (*this)(i, j) + other(i, j);
        }
    }
    return result;
}

Dense Dense::operator*(double alpha) const {
    Dense result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result(i, j) = (*this)(i, j) * alpha;
        }
    }
    return result;
}

std::vector<double> Dense::operator*(const std::vector<double>& vec) const {
    std::vector<double> result(rows);
    for (int i = 0; i < rows; ++i) {
        result[i] = 0;
        for (int j = 0; j < cols; ++j) {
            result[i] += (*this)(i, j) * vec[j];
        }
    }
    return result;
}

Dense Dense::operator*(const Dense& other) const {
    Dense result(rows, other.cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            result(i, j) = 0;
            for (int k = 0; k < cols; ++k) {
                result(i, j) += (*this)(i, k) * other(k, j);
            }
        }
    }
    return result;
}

void Dense::print() const {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << (*this)(i, j) << " ";
        }
        std::cout << std::endl;
    }
}


