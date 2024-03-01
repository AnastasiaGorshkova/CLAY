#include <utility>
#include <cmath>
#include <Matrix & Vector/Dense.hpp>
#include <Matrix & Vector/vector.hpp>

std::pair<Dense, Dense> Householder(const Dense& A) {
    Dense R = A;
    std::size_t size = R.get_rows();
    std::vector<double> e(size, 0); // vector e1
    e[0] = 1;
    std::vector<double> q(R.get_rows() * R.get_cols(), 0);
    std::vector<double> v = R.get_column(0) + sign(R.get_column(0)) * modul(R.get_column(0)) * e;
    for (std::size_t i = 0; i < R.get_rows() * R.get_cols(); i++){
        for(std::size_t j = 0; j < R.get_cols(); j++){
            if(i == j)
                q[i * R.get_cols()+ j] = 1;
        }
    }
    q = q - 2 * transposed(v)/ sq(v);
    Dense Q(q, R.get_rows(), R.get_cols()); //calculating Q[0]
    
    std::vector<double> prod = R.get_column(0) - 2 * v * ((v * R.get_column(0)) / pow(modul(v), 2));
    
    R.swap_elements(0, 0, prod[0]);
    for (std::size_t i = 1; i < R.get_cols(); i++) {
        R.swap_elements(i, 0, 0); //transforming the first column of R
    }
    for (std::size_t j = 1; j < R.get_cols(); j++) {
        prod = R.get_column(j) - 2 * v * ((v * R.get_column(j)) / pow(modul(v), 2));
        for (std::size_t i = 0; i < size; i++) {
            R.swap_elements(i, j, prod[i]); //transforming rest of the columns
        }
    }
    for (std::size_t n = 1; n < R.get_cols(); n++) { //repeating for other columns
        size--;
        std::vector<double> e1(size, 0); //new vector e1
        e1[0] = 1;
        std::vector<double> a = R.get_column(n);
        int t = static_cast<int>(n);
        std::vector<double> x = std::vector<double>(a.begin() + t, a.end());
        v = x + sign(x) * modul(x) * e1;
        prod = x - 2 * v * ((v * x) / pow(modul(v), 2));
        R.swap_elements(n, n, prod[0]); //transforming n-th column
        
        for (std::size_t i = 1 + n; i < size + 1; i++) {
            R.swap_elements(i, n, 0);
        }
        for (std::size_t j = 1 + n; j < R.get_cols(); j++) { //transforming columns after n
            a = R.get_column(j);
            x = std::vector<double>(a.begin() + t, a.end());
            prod = x - 2 * v * ((v * x) / pow(modul(v), 2));
            for (std::size_t i = 0; i < size; i++) {
                R.swap_elements(i + n, j, prod[i]);
            }
        }
        for (std::size_t i = 0; i < R.get_cols() ; i++) {
            a = Q.get_row(i);
            x = std::vector<double>(a.begin() + t, a.end());
            for (std::size_t j = 0; j < size; j++) {
                prod = x - 2 * v * ((v * x) / pow(modul(v), 2));
                Q.swap_elements(i, j + n, prod[j]);
            }
        }
    }
    return std::make_pair(Q,R);
}
