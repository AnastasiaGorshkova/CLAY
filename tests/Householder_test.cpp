#include<gtest/gtest.h>
#include <Matrix & Vector/Dense.hpp>
#include <Solvers/Householder.cpp>

TEST(Householder, TEST_1)
{
	std::vector<double> a = {1, 3, 2, 4, 7, 6, 2, 1, 0};
	Dense A(a, 3, 3);
	Dense Q = Householder(A).first;
    	Dense R = Householder(A).second;
    	
	std::vector<double> q = {0.218, 0.535, -0.816, 0.873, 0.267, 0.408, 0.436, -0.802, -0.408};
	std::vector<double> r = {4.583, 7.201, 5.674, 0, 2.673, 2.673, 0, 0, 0.816};
	Dense B(q, 3, 3);
	Dense C(r, 3, 3);
	for (std::size_t i = 0; i < 3; ++i) {
        	for (std::size_t j = 0; j < 3; ++j) {
       			ASSERT_NEAR(Q(i,j), B(i,j) , 0.001);
       			ASSERT_NEAR(R(i,j), C(i,j) , 0.001);
        	}
    	}
}

TEST(Householder, TEST_2)
{
	std::vector<double> a = {5.2, 0.15, 3.1, 7.8, 1.4, 2.3, 3.67, 1.57, 0.24};
	Dense A(a, 3, 3);
	Dense Q = Householder(A).first;
    	Dense R = Householder(A).second;
    	
	std::vector<double> q = {0.517, -0.622, 0.589, 0.775, 0.047, -0.630, 0.365, 0.782, 0.506};
	std::vector<double> r = {10.067, 1.735, 3.471, 0, 1.200, -1.633, 0, 0, 0.496};
	Dense B(q, 3, 3);
	Dense C(r, 3, 3);
	for (std::size_t i = 0; i < 3; ++i) {
        	for (std::size_t j = 0; j < 3; ++j) {
       			ASSERT_NEAR(Q(i,j), B(i,j) , 0.001);
       			ASSERT_NEAR(R(i,j), C(i,j) , 0.001);
        	}
    	}
}
