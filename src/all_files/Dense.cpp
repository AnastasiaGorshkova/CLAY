#include "Dense.hpp"

Dense::Dense(const std::vector<double>& data, std::size_t rows, std::size_t cols) : rows(rows), cols(cols), data(data) {}

double Dense::operator()(std::size_t i, std::size_t j) const {
    return data[i*cols + j];
}

int Dense::get_cols() const {
	return cols;
}
	    
int Dense::get_rows() const {
	return rows;
}


Dense Dense::operator+(const Dense& other) const {
    std::vector<double> result; 
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            result.push_back(data[i * cols + j] + other(i,j));
        }
    }
    Dense ans(result, cols, rows);
    return ans;
}

Dense Dense::operator*(double alpha) const {
    std::vector<double> result; 
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            result.push_back(data[i * cols + j] * alpha);
        }
    }
    Dense ans(result, cols, rows);
    return ans;
}

std::vector<double> Dense::operator*(const std::vector<double>& vec) const {
    std::vector<double> result;
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            result.push_back(data[i * cols + j] * vec[j]);
        }
    } 
    return result;
}

Dense Dense::operator*(const Dense& other) const {
    std::vector<double> result;
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < other.cols; ++j) {
            for (std::size_t k = 0; k < cols; ++k) {
                result.push_back(data[i * cols + k] * other(k, j));
            }
        }
    }
    Dense ans(result, cols, rows);
    return ans;
}

void Dense::print() const {
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            std::cout << data[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }
}
