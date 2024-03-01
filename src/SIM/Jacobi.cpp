#include <Matrix & Vector/CSR.hpp>
#include <Matrix & Vector/vector.hpp>

std::vector<double> Jacobi( CSR& A, const std::vector<double>& b, const std::vector<double>& x0, double tolerance) {
    	std::vector<double> x = x0; // x0 - начальное приближение
    	std::vector<double> res = A * x - b;
    	std::vector<double> LU;
    	while (modul(res) > tolerance){
	    	
	    	LU = A.Jakobi_non_diagonal(x);
	    	x = mul(A.Jakobi_reverse_diag(), (b - LU));
	    	res = A * x - b;
    	}   	
        return x;
}
