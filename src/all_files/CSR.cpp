#include "CSR.hpp"

CSR::CSR(const std::vector<double>& values, const std::vector<unsigned int>& cols, const std::vector<unsigned int>& rows) : values(values), cols(cols), rows(rows) {}


CSR::CSR(const Dense& D) {
	std::size_t temp = 0;
	rows.push_back(0);
	for (std::size_t i = 0; i < D.get_rows(); i++) {
		for (std::size_t j = 0; j < D.get_cols(); j++) {
			if (D(i, j) != 0) {
				values.push_back(D(i, j));
				cols.push_back(j);
				temp++;
			}
		}
		rows.push_back(temp);
	}

}

double CSR::operator()(std::size_t i, std::size_t j) const {
    for (std::size_t k = rows[i]; k < rows[i + 1]; ++k) {
        if (cols[k] == j) {
            return values[k];
        }
    }
    return 0;
}

CSR CSR::operator*(double alpha) const {
    std::vector<double> result;
    for (std::size_t i = 0; i < values.size(); ++i) {
        result.push_back(values[i]* alpha);
    }
    CSR res(result, cols, rows);
    return res;
}

std::vector<double> CSR::operator*(const std::vector<double>& vec) const {
    std::vector<double> result;
    double temp;
    for (std::size_t i = 0; i < rows.size() -1; ++i) {
    	temp = 0;
        for (std::size_t k = rows[i]; k < rows[i + 1]; ++k) {
            temp += values[k] * vec[cols[k]];
        }
        result.push_back(temp);
    }
    return result;
}

void CSR::print() const {
        std::cout << "values: ";
        for(std::size_t i = 0; i < values.size(); i++) std::cout << values[i] << " ";
        std::cout << std::endl;
        std::cout << "cols: ";
        for(std::size_t i = 0; i < cols.size(); i++) std::cout << cols[i] << " ";
        std::cout << std::endl;
        std::cout << "rows: ";
        for(std::size_t i = 0; i < rows.size(); i++) std::cout << rows[i] << " ";
        std::cout << std::endl;    
}
