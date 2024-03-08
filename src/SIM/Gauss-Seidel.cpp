#include "SIM.hpp"

std::vector<double> Gauss_Seidel(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, double tolerance) {
    	std::vector<double> x = x0; // x0 - начальное приближение
    	std::vector<double> res = A * x - b;
    	double diag; // диагональный элемент
    	
    	std::vector<unsigned int> c = A.get_cols();
    	std::vector<unsigned int> r = A.get_rows();
    	std::vector<double> val = A.get_values();
    	
    	while (modul(res) > tolerance){
	for (std::size_t i = 0; i < x.size(); ++i) {
	    x[i] = b[i];
	    for(std::size_t k = r[i]; k < r[i+1]; k++){
               if (c[k] == i){
                    diag = A(i,i);
                    continue;
               }
               x[i] -= val[k] * x[c[k]];
            }
            x[i] /= diag;
        }
    	res = A * x - b;
    	}	
        return x;
}
