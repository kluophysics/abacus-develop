#ifndef TENSOR3_H
#define TENSOR3_H

#include <armadillo>

namespace Module_Optimizer
{
    template <typename T>
    class Tensor3
    {
    public:
        Tensor3();
        Tensor3(const Tensor3<T> &);
        Tensor3(unsigned int nrows, unsigned int ncols, unsigned int nslices);
        
        // uni operator overload
        Tensor3<T> operator-() const; // negative operator
        
        // assignment operator overload
        Tensor3<T> & operator=(const Tensor3<T>& p) ;
        Tensor3<T> & operator-=(const Tensor3<T>& p);
        Tensor3<T> & operator+=(const Tensor3<T>& p);
        Tensor3<T> & operator*=(const Tensor3<T>& p);

        Tensor3 operator+( T s) const; // add p with s, p+s;
        Tensor3 operator-( T s) const; // p-s;
        Tensor3 operator*( T s) const; //  p*s;

        Tensor3 operator+( const Tensor3<T> & p) const; // add p with s, p+s;
        Tensor3 operator-( const Tensor3<T> & p) const; // add p with s, p+s;
        Tensor3 operator*( const Tensor3<T> & p) const; // add p with s, p+s;

        friend Tensor3 operator+(T s, const Tensor3<T> & p); // s+p
        friend Tensor3 operator*(T s, const Tensor3<T> & p); // s*p
        friend Tensor3 operator-(T s, const Tensor3<T> & p); // s-p

        void set_zeros();
        void clear();
        void print(const char* s) const;
        void brief_print(const char* s) const;
        void resize(int k, int r, int c);

        double norm();
        Tensor3 <T> transpose() const; // conjugate transpose.
        Tensor3 <T> inv() const; // inverse.

        inline unsigned int get_nrows() const { return nrows;}; 
        inline unsigned int get_ncols() const { return ncols;}; 
        inline unsigned int get_nslices() const { return nslices; }; 
        inline unsigned int get_size() const { return size;}; 

        // check if two matrices a and b are of the same dimensions
        friend bool equal_dimension(const Tensor3<T> & a, const Tensor3<T>  & b); // 

        private:
        unsigned int nrows; 
        unsigned int ncols;
        unsigned int nslices;
        unsigned int size; // size = nrows * ncols * nslices

        arma::Cube<T> data;
    };
}

#endif //TENSOR3_H