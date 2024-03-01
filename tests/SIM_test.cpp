#include<gtest/gtest.h>
#include <Matrix & Vector/Dense.hpp>
#include <SIM/Simple_iteration.cpp>
#include <SIM/Jacobi.cpp>
#include <SIM/Gauss-Seidel.cpp>

TEST(Simple_iteration, TEST_1)
{
	std::vector<double> a = {3, 1, 0, 2, 5, 1, 6, 0, 7};
	std::vector<double> x = {2.5, 1, 3};
	CSR A(a, 3, 3);
	double tolerance = pow(10,-13);
	double tau = 0.111;
	std::vector<double> x0 = {-15, 100, 567};
	std::vector<double> b = A * x;
	
	std::vector<double> result = Simple(A, b, x0, tau, tolerance);	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

TEST(Gauss_Seidel, TEST_1)
{
	std::vector<double> a = {33.7, 2, 11, 20, 42, 3, 1, 0, 7};
	std::vector<double> x = {0.47, 5.24, 3.18};
	CSR A(a, 3, 3);
	
	double tolerance = pow(10,-12);
	std::vector<double> x0 = {-24, 5.678, 517};
	
	std::vector<double> b = A * x;
	std::vector<double> result = Gauss_Seidel(A, b, x0, tolerance);
	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

TEST(Jacobi, TEST_1)
{
	std::vector<double> a = {4, 1, 1, 0, 5, 2, 1, 4, 7};
	std::vector<double> x = {1, 1, 1};
	CSR A(a, 3, 3);
	
	double tolerance = pow(10,-12);
	std::vector<double> x0 = {2, 2, 2};
	
	std::vector<double> b = A * x;
	std::vector<double> result = Jacobi(A, b, x0, tolerance);
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

