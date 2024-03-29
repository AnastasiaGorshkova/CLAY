#include <iostream>
#include <fstream>
#include <vector>

#include "../src/SIM/SIM.hpp"
#include <chrono>

using high_resolution_clock = std::chrono::high_resolution_clock;

void get_time_Jacobi(){
	std::ofstream file;
        file.open("res_Jacobi.csv"); //файл для записи данных
	std::vector<double> a = {250, 2.75, 1.5, 3.5, 8.6, 2.75, 360, 0.63, 6.75, 2.3, 1.5, 0.63, 165, 5.4, 1.12, 2.5, 6.75, 5.4, 95, 0.25, 8.6, 2.3, 1.12, 0.25, 175};
	std::vector<double> x = {3.5, 1.75, 2.1, 0.7, 4.9};
	std::vector<double> x0 = {174, 33.7, 9.11, 27.3, 45.5};
	
	CSR A(a, 5, 5);

	std::vector<double> b = A * x;  // вектор свободных членов
	double tolerance = pow(10,-12);
		
	std::vector<double> y = x0; // x0 - начальное приближение
	std::vector<double> res = A * y - b;
	std::vector<double> LU;
		
	std::size_t n = 0;
	auto begin = high_resolution_clock::now(); // начало записи времени
    	while (modul(res) > tolerance){
		LU = A.Jakobi_non_diagonal(y);
		y = mul(A.Jakobi_reverse_diag(), (b - LU));
		res = A * y - b;
		    	
		auto now = high_resolution_clock::now();
		auto t = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
		file << modul(res) << " " << n << " "<< t.count() << std::endl;
		    	
		n++;
	} 	
	
	auto end = high_resolution_clock::now(); // конец записи времени
	auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	file << time.count() << std::endl;
	file.close();
}


void get_time_Gauss_Seidel(){
	std::ofstream file;
        file.open("res_Gauss_Seidel.csv"); //файл для записи данных
	std::vector<double> a = {250, 2.75, 1.5, 3.5, 8.6, 2.75, 360, 0.63, 6.75, 2.3, 1.5, 0.63, 165, 5.4, 1.12, 2.5, 6.75, 5.4, 95, 0.25, 8.6, 2.3, 1.12, 0.25, 175};
	std::vector<double> x = {3.5, 1.75, 2.1, 0.7, 4.9};
	std::vector<double> x0 = {174, 33.7, 9.11, 27.3, 45.5};
	
	CSR A(a, 5, 5);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
	    	double diag; // диагональный элемент
		
		std::size_t n = 0;
		auto begin = high_resolution_clock::now(); // начало записи времени
    		while(modul(res) > tolerance){
		     for (std::size_t c = 0; c < y.size(); ++c) {
			 y[c] = b[c];
			 for(std::size_t k = A.get_rows()[c]; k < A.get_rows()[c+1]; k++){
			    if (A.get_cols()[k] == c){
			       diag = A(c,c);
			       continue;
			    }
			    
			    y[c] -= A.get_values()[k] * y[A.get_cols()[k]];
			 }
			 
			 y[c] /= diag;
		    }
		    
		    res = A * y - b;
		    
		    auto now = high_resolution_clock::now();
		    auto t = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
		    file << modul(res) << " " << n << " "<< t.count() << std::endl;
		    	n++;
		}	
	auto end = high_resolution_clock::now(); // конец записи времени
	auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	file << time.count() << std::endl;
	
	file.close();
}



void get_time_Simple_iteration(){
	std::ofstream file;
        file.open("res_Simple.csv"); //файл для записи данных
	std::vector<double> a = {250, 2.75, 1.5, 3.5, 8.6, 2.75, 360, 0.63, 6.75, 2.3, 1.5, 0.63, 165, 5.4, 1.12, 2.5, 6.75, 5.4, 95, 0.25, 8.6, 2.3, 1.12, 0.25, 175};
	std::vector<double> x = {3.5, 1.75, 2.1, 0.7, 4.9};
	std::vector<double> x0 = {174, 33.7, 9.11, 27.3, 45.5};
	
	CSR A(a, 5, 5);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		double tau = 0.00439898;
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
		
		std::size_t n = 0;
    		auto begin = high_resolution_clock::now(); // начало записи времени
    		
    		while(modul(res) > tolerance){
			y = y - tau * res;
	    		res = A * y - b;
		        
		        auto now = high_resolution_clock::now();
		    	auto t = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
		    	file << modul(res) << " " << n << " "<< t.count() << std::endl;
		    	n++;
		        
         	}
		auto end = high_resolution_clock::now(); // конец записи времени
		auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
		file << time.count() << std::endl;
	file.close();
}



void get_time_Chebyshev(){
	std::ofstream file;
        file.open("res_Chebyshev.csv"); //файл для записи данных
	std::vector<double> a = {250, 2.75, 1.5, 3.5, 8.6, 2.75, 360, 0.63, 6.75, 2.3, 1.5, 0.63, 165, 5.4, 1.12, 2.5, 6.75, 5.4, 95, 0.25, 8.6, 2.3, 1.12, 0.25, 175};
	std::vector<double> x = {3.5, 1.75, 2.1, 0.7, 4.9};
	std::vector<double> x0 = {174, 33.7, 9.11, 27.3, 45.5};
	
	CSR A(a, 5, 5);
	std::vector<double> b = A * x;  // вектор свободных членов
	double tolerance = pow(10,-12);
	double lambda_max = 360.282;
	double lambda_min = 94.369;
	std::size_t res = 6;
	
		double R = static_cast<double>(res);
		double n = pow(2,R); //количество прогонок
		
		std::size_t N = static_cast<std::size_t>(n);
		
		std::vector<std::size_t > res_num(N,0);
		std::size_t P = static_cast<std::size_t>(pow(2,R-1));
	    	res_num[P] = 1;
	    	
	    	for(std::size_t  i = 2; i <= res; i++){
	    		int I = static_cast<int>(i);
	    		std::size_t P1 = static_cast<std::size_t>(pow(2,R-I));
			std::size_t k = P1;
			for(std::size_t  j = 0; j < N; j+= 2*k){
				std::size_t P2 = static_cast<std::size_t>(pow(2,I));
		    		res_num[j+k] = P2 - res_num[j] - 1;
			}
	    	}
	    
	    	double cos_a = cos(M_PI/n);
	    	double sin_a = sin(M_PI/n);
	    	double sin_b = sin(M_PI/(2*n));
	    
	    	std::vector<double> z(N,0);
	    	z[0] = sin(M_PI/(2*n));

	    	for(std::size_t  i = 1; i < N; i++){
			z[i] = z[i-1]*cos_a - sin_b*sin_a;
			sin_b = z[i-1]*cos_a + sin_b*sin_a;
	       	 	z[i-1] = z[i-1]*(lambda_max-lambda_min)/2 + (lambda_max+lambda_min)/2;
	    	}
	    
	    	z[N-1] = z[N-1]*(lambda_max-lambda_min)/2 + (lambda_max+lambda_min)/2;



    	std::vector<double> y = x0; // x0 - начальное приближение
    	std::vector<double> r = A * y - b;
   
   	std::size_t T = 0;
    	auto begin = high_resolution_clock::now(); // начало записи времени
   
    	while(modul(r) > tolerance){
        	for(std::size_t i = 0; i < res_num.size(); i++) {
            		y = y - r * (1 / z[res_num[i]]);
            		r = A * y - b;
            		auto now = high_resolution_clock::now();
		    	auto t = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
		    	file << modul(r) << " " << T << " "<< t.count() << std::endl;
		    	T++;
        	}
    	}
    	
    	auto end = high_resolution_clock::now(); // конец записи времени
	auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	file << time.count() << std::endl;
	file.close();
}

int main(){
    get_time_Jacobi();
    get_time_Gauss_Seidel();
    get_time_Simple_iteration();
    get_time_Chebyshev();
}
