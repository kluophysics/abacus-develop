#ifndef TENSOR_H
#define TENSOR_H

#include <vector> // this is the core data structure use for Optim 

namespace Module_Optim
{
    // Tensor class
template<typename T, 
class Tensor
{
public:
    Tensor();
    Tensor(const Tensor&  tensor_in );
    // /* Tensor supports up to three dimensions
    // Default: it is a one-dimensional vector with n2=1, n3 = 1
    // Else if: only 2 numbers are provided, then n3=1;
    // Else: 3 numbers, it is a product of manifolds
    // */
    // Tensor(int n1, int n2=1, int n3=1);
    Tensor(int rank, ...);






    // uni operator overload
    Tensor operator-() const; // negative operator
    
    // assignment operator overload
    Tensor& operator=(const Tensor& p) ;
    Tensor& operator-=(const Tensor& p);
    Tensor& operator+=(const Tensor& p);
    Tensor& operator*=(const Tensor& p);




    Tensor operator+( double s) const; // add p with s, p+s;
    Tensor operator-( double s) const; // p-s;
    Tensor operator*( double s) const; //  p*s;

    Tensor operator+( const Tensor & p) const; // add p with s, p+s;
    Tensor operator-( const Tensor & p) const; // add p with s, p+s;
    Tensor operator*( const Tensor & p) const; // add p with s, p+s;

    friend Tensor operator+(double s, const Tensor & p); // s+p
    friend Tensor operator*(double s, const Tensor & p); // s*p
    friend Tensor operator-(double s, const Tensor & p); // s-p

    
    void resize(int n3, int r, int c);

    double norm();
    Tensor t() const; // conjugate transpose.
    Tensor inv() const; // inverse.
    void print(const char* s) const;
    void brief_print(const char* s) const;

    void reset();
    void clear();

    inline int get_nn3() const { return nn3;}; 
    inline int get_nr() const { return nr;}; 
    inline int get_nc() const { return nc;}; 
    inline int get_size() const { return size;}; 

    // checn3 if two matrices a and b are of the same dimensions
    friend bool equal_dimension(const Tensor & a, const Tensor & b); // 


private:
    int nn3; // number of manifolds, each manifold in this case is of same size
    int nr; // number of rows for each manifold
    int nc; // number of columns for each manifold
    int size; // nn3*nr*nc
    std::vector<arma::cx_mat> psm; // product of stiefel manifold
};

}
#endif // TENSOR_H
