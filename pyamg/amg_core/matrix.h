	// Copyright (c) 2015-2017, RAPtor Developer Team
// License: Simplified BSD, http://opensource.org/licenses/BSD-2-Clause
#ifndef RAPTOR_CORE_MATRIX_HPP
#define RAPTOR_CORE_MATRIX_HPP

#include "types.h"
#include "vector.h"

/**************************************************************
 *****   Matrix Base Class
 **************************************************************
 ***** This class constructs a sparse matrix, supporting simple linear
 ***** algebra operations.
 *****
 ***** Attributes
 ***** -------------
 ***** n_rows : int
 *****    Number of rows
 ***** n_cols : int
 *****    Number of columns
 ***** nnz : int
 *****    Number of nonzeros
 ***** idx1 : std::vector<int>
 *****    List of position indices, specific to type of matrix
 ***** idx2 : std::vector<int>
 *****    List of position indices, specific to type of matrix
 ***** vals : std::vector<double>
 *****    List of values in matrix
 *****
 ***** Methods
 ***** -------
 ***** resize(int n_rows, int n_cols)
 *****    Resizes dimension of matrix to passed parameters
 ***** mult(Vector* x, Vector* b)
 *****    Sparse matrix-vector multiplication b = A * x
 ***** residual(Vector* x, Vector* b, Vector* r)
 *****    Calculates the residual r = b - A * x
 *****
 ***** Virtual Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (COO, CSR, CSC)
 ***** sort()
 *****    Sorts the matrix by position.  Whether row-wise or 
 *****    column-wise depends on matrix format.
 ***** add_value(int row, int col, double val)
 *****     Adds val to position (row, col)
 *****     TODO -- make sure this is working for CSR/CSC
 **************************************************************/
namespace raptor
{
  // Forward Declaration of classes so objects can be used
  class COOMatrix;
  class CSRMatrix;
  class CSCMatrix;
  class BSRMatrix;

  class Matrix
  {

  public:

    /**************************************************************
    *****   Matrix Base Class Constructor
    **************************************************************
    ***** Sets matrix dimensions, and sets nnz to 0
    *****
    ***** Parameters
    ***** -------------
    ***** _nrows : int
    *****    Number of rows in matrix
    ***** _ncols : int
    *****    Number of cols in matrix
    **************************************************************/
    Matrix(int _nrows, int _ncols)
    {
        n_rows = _nrows;
        n_cols = _ncols;
        nnz = 0;
        sorted = false;
        diag_first = false;
    }

    /**************************************************************
    *****   Matrix Base Class Constructor
    **************************************************************
    ***** Sets matrix dimensions and nnz based on Matrix* A
    *****
    ***** Parameters
    ***** -------------
    ***** A : Matrix*
    *****    Matrix to be copied
    **************************************************************/
    Matrix()
    {
        n_rows = 0;
        n_cols = 0;
        nnz = 0;
        sorted = false;
        diag_first = false;
    }

/*
    virtual ~Matrix(){}

    virtual format_t format() = 0;
    virtual void sort() = 0;
    virtual void move_diag() = 0;
    virtual void remove_duplicates() = 0;
    virtual void add_value(int row, int col, double val) = 0;

    virtual void print() = 0;

    virtual void copy(const COOMatrix* A) = 0;
    virtual void copy(const CSRMatrix* A) = 0;
    virtual void copy(const CSCMatrix* A) = 0;
    virtual void copy(const BSRMatrix* A) = 0;
*/
    void jacobi(Vector& x, Vector& b, Vector& tmp, double omega = .667);
    void gauss_seidel(Vector& x, Vector& b);
    void SOR(Vector& x, Vector& b, double omega = .667);

    Matrix* strength(double theta = 0.0);
    Matrix* aggregate();

    void mult(Vector& x, Vector& b)
    {
        mult(x.values, b.values);
    }
    void mult(std::vector<double>& x, Vector& b)
    {
        mult(x, b.values);
    }
    void mult(Vector& x, std::vector<double>& b)
    {
        mult(x.values, b);
    }
//    virtual void mult(std::vector<double>& x, std::vector<double>& b) = 0;
     void mult(std::vector<double>& x, std::vector<double>& b);

    void mult_T(Vector& x, Vector& b)
    {
        mult_T(x.values, b.values);
    }
    void mult_T(std::vector<double>& x, Vector& b)
    {
        mult_T(x, b.values);
    }
    void mult_T(Vector& x, std::vector<double>& b)
    {
        mult_T(x.values, b);
    }
//    virtual void mult_T(std::vector<double>& x, std::vector<double>& b) = 0;
    void mult_T(std::vector<double>& x, std::vector<double>& b);

    void mult_append(Vector& x, Vector& b)
    {
        mult_append(x.values, b.values);
    }
    void mult_append(std::vector<double>& x, Vector& b)
    {
        mult_append(x, b.values);
    }
    void mult_append(Vector& x, std::vector<double>& b)
    {
        mult_append(x.values, b);
    }
//    virtual void mult_append(std::vector<double>& x, std::vector<double>& b) = 0;
    void mult_append(std::vector<double>& x, std::vector<double>& b);

    void mult_append_T(Vector& x, Vector& b)
    {
        mult_append_T(x.values, b.values);
    }
    void mult_append_T(std::vector<double>& x, Vector& b)
    {
        mult_append_T(x, b.values);
    }
    void mult_append_T(Vector& x, std::vector<double>& b)
    {
        mult_append_T(x.values, b);
    }
//    virtual void mult_append_T(std::vector<double>& x, std::vector<double>& b) = 0;
    void mult_append_T(std::vector<double>& x, std::vector<double>& b);

    void mult_append_neg(Vector& x, Vector& b)
    {
        mult_append_neg(x.values, b.values);
    }
    void mult_append_neg(std::vector<double>& x, Vector& b)
    {
        mult_append_neg(x, b.values);
    }
    void mult_append_neg(Vector& x, std::vector<double>& b)
    {
        mult_append_neg(x.values, b);
    }
//    virtual void mult_append_neg(std::vector<double>& x, std::vector<double>& b) = 0;
    void mult_append_neg(std::vector<double>& x, std::vector<double>& b) ;

    void mult_append_neg_T(Vector& x, Vector& b)
    {
        mult_append_neg_T(x.values, b.values);
    }
    void mult_append_neg_T(std::vector<double>& x, Vector& b)
    {
        mult_append_neg_T(x, b.values);
    }
    void mult_append_neg_T(Vector& x, std::vector<double>& b)
    {
        mult_append_neg_T(x.values, b);
    }
//    virtual void mult_append_neg_T(std::vector<double>& x, std::vector<double>& b) = 0;
    void mult_append_neg_T(std::vector<double>& x, std::vector<double>& b);

    void residual(const Vector& x, const Vector& b, Vector& r)
    {
        residual(x.values, b.values, r.values);
    }
    void residual(const std::vector<double>& x, const Vector& b, Vector& r)
    {
        residual(x, b.values, r.values);
    }
//    virtual void residual(const std::vector<double>& x, const std::vector<double>& b,
//            std::vector<double>& r) = 0;
    void residual(const std::vector<double>& x, const std::vector<double>& b,
            std::vector<double>& r) ;

    CSRMatrix* mult(const CSRMatrix* B){ return NULL; }
    CSRMatrix* mult(const CSCMatrix* B){ return NULL; }
    CSRMatrix* mult(const COOMatrix* B){ return NULL; }
    CSRMatrix* mult_T(const CSRMatrix* A){ return NULL; }
    CSRMatrix* mult_T(const CSCMatrix* A){ return NULL; }
    CSRMatrix* mult_T(const COOMatrix* A){ return NULL; }

    void RAP(const CSCMatrix& P, CSCMatrix* Ac);
    void RAP(const CSCMatrix& P, CSRMatrix* Ac);

//	virtual Matrix* ilu_k(int lof) = 0;
	Matrix* ilu_k(int lof);
	//virtual Matrix* ilu_levels() = 0;
	//virtual Matrix* ilu_sparsity(Matrix* levls, int lof) = 0;
	//virtual Matrix* ilu_symbolic(int lof) = 0;
	//virtual std::vector<double> ilu_numeric(Matrix* levls) = 0;

    Matrix* subtract(Matrix* B);

//    virtual void add_block(int row, int col, std::vector<double>& values) = 0;
    void add_block(int row, int col, std::vector<double>& values);

    void resize(int _n_rows, int _n_cols);

//    virtual Matrix* transpose() = 0;
    Matrix* transpose();

    std::vector<int>& index1()
    {
        return idx1;
    }

    std::vector<int>& index2()
    {
        return idx2;
    }
    
    std::vector<double>& values()
    {
        return vals;
    }

    std::vector<int> idx1;
    std::vector<int> idx2;
    std::vector<double> vals;

    int n_rows;
    int n_cols;
    int nnz;

    bool sorted;
    bool diag_first;

  };


/**************************************************************
 *****   COOMatrix Class (Inherits from Matrix Base Class)
 **************************************************************
 ***** This class constructs a sparse matrix in COO format.
 *****
 ***** Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (COO)
 ***** sort()
 *****    Sorts the matrix by row, and by column within each row.
 ***** add_value(int row, int col, double val)
 *****     Adds val to position (row, col)
 ***** rows()
 *****     Returns std::vector<int>& containing the rows corresponding
 *****     to each nonzero
 ***** cols()
 *****     Returns std::vector<int>& containing the cols corresponding
 *****     to each nonzero
 ***** data()
 *****     Returns std::vector<double>& containing the nonzero values
 **************************************************************/
  class COOMatrix : public Matrix
  {

  public:

    /**************************************************************
    *****   COOMatrix Class Constructor
    **************************************************************
    ***** Initializes an empty COOMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** _nrows : int
    *****    Number of rows in Matrix
    ***** _ncols : int
    *****    Number of columns in Matrix
    ***** nnz_per_row : int
    *****    Prediction of (approximately) number of nonzeros 
    *****    per row, used in reserving space
    **************************************************************/
    COOMatrix(int _nrows, int _ncols, int nnz_per_row = 1) : Matrix(_nrows, _ncols)
    {
        if (nnz_per_row)
        {
            int _nnz = nnz_per_row * _nrows;
            if (_nnz)
            {
                idx1.reserve(_nnz);
                idx2.reserve(_nnz);
                vals.reserve(_nnz);
            }
        }        
    }

    COOMatrix(int _nrows, int _ncols, double* _data) : Matrix(_nrows, _ncols)
    {
        nnz = 0;
        int nnz_dense = n_rows*n_cols;

        if (nnz_dense)
        {
            idx1.reserve(nnz_dense);
            idx2.reserve(nnz_dense);
            vals.reserve(nnz_dense);
        }

        for (int i = 0; i < n_rows; i++)
        {
            for (int j = 0; j < n_cols; j++)
            {
                double val = _data[i*n_cols + j];
                if (fabs(val) > zero_tol)
                {
                    idx1.push_back(i);
                    idx2.push_back(j);
                    vals.push_back(val);
                    nnz++;
                }
            }
        }
    }

    COOMatrix(int _nrows, int _ncols, std::vector<int>& rows, std::vector<int>& cols, 
            std::vector<double>& data) : Matrix(_nrows, _ncols)
    {
        nnz = idx1.size();
        idx1.resize(nnz);
        idx2.resize(nnz);
        vals.resize(nnz);

        std::copy(rows.begin(), rows.end(), idx1.begin());
        std::copy(cols.begin(), cols.end(), idx2.begin());
        std::copy(data.begin(), data.end(), vals.begin());
    }

    COOMatrix()
    {
    }


    /**************************************************************
    *****   COOMatrix Class Constructor
    **************************************************************
    ***** Constructs a COOMatrix from a CSRMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const CSRMatrix*
    *****    CSRMatrix A, from which to copy data
    **************************************************************/
    explicit COOMatrix(const CSRMatrix* A)
    {
        copy(A);
    }

    /**************************************************************
    *****   COOMatrix Class Constructor
    **************************************************************
    ***** Copies matrix, constructing new COOMatrix from 
    ***** another COOMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const COOMatrix*
    *****    COOMatrix A, from which to copy data
    **************************************************************/
    explicit COOMatrix(const COOMatrix* A)
    {
        copy(A);
    }

    /**************************************************************
    *****   COOMatrix Class Constructor
    **************************************************************
    ***** Constructs a COOMatrix from a CSCMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const CSCMatrix*
    *****    CSCMatrix A, from which to copy data
    **************************************************************/
    explicit COOMatrix(const CSCMatrix* A)
    {
        copy(A);
    }

    /**************************************************************
    *****   COOMatrix Class Constructor
    **************************************************************
    ***** Constructs a COOMatrix from a BSRMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const BSRMatrix*
    *****    BSRMatrix A, from which to copy data
    **************************************************************/
    explicit COOMatrix(const BSRMatrix* A)
    {
        copy(A);
    }

    ~COOMatrix()
    {

    }

    Matrix* transpose();

    void print();

    void copy(const COOMatrix* A);
    void copy(const CSRMatrix* A);
    void copy(const CSCMatrix* A);
    void copy(const BSRMatrix* A);
    void block_copy(const BSRMatrix* A, int row, int num_blocks_prev, int col);

    std::vector<double> to_dense() const;

    void add_value(int row, int col, double value);
    void sort();
    void move_diag();
    void remove_duplicates();

    template <typename T, typename U> void mult(T& x, U& b)
    { 
        Matrix::mult(x, b);
    }
    template <typename T, typename U> void mult_T(T& x, U& b)
    { 
        Matrix::mult_T(x, b);
    }
    template <typename T, typename U> void mult_append(T& x, U& b)
    { 
        Matrix::mult_append(x, b);
    }
    template <typename T, typename U> void mult_append_T(T& x, U& b)
    { 
        Matrix::mult_append_T(x, b);
    }
    template <typename T, typename U> void mult_append_neg(T& x, U& b)
    { 
        Matrix::mult_append_neg(x, b);
    }
    template <typename T, typename U> void mult_append_neg_T(T& x, U& b)
    { 
        Matrix::mult_append_neg_T(x, b);
    }
    template <typename T, typename U, typename V> 
    void residual(const T& x, const U& b, V& r)
    { 
        Matrix::residual(x, b, r);
    }

    void mult(std::vector<double>& x, std::vector<double>& b)
    {
        for (int i = 0; i < n_rows; i++)
            b[i] = 0.0;
        mult_append(x, b);
    }
    void mult_T(std::vector<double>& x, std::vector<double>& b)
    {
        for (int i = 0; i < n_cols; i++)
            b[i] = 0.0;

        mult_append_T(x, b);
    }
    void mult_append(std::vector<double>& x, std::vector<double>& b)
    { 
        for (int i = 0; i < nnz; i++)
        {
            b[idx1[i]] += vals[i] * x[idx2[i]];
        }
    }
    void mult_append_T(std::vector<double>& x, std::vector<double>& b)
    {
        for (int i = 0; i < nnz; i++)
        {
            b[idx2[i]] += vals[i] * x[idx1[i]];
        }
    }
    void mult_append_neg(std::vector<double>& x, std::vector<double>& b)
    {
        for (int i = 0; i < nnz; i++)
        {
            b[idx1[i]] -= vals[i] * x[idx2[i]];
        }
    }
    void mult_append_neg_T(std::vector<double>& x, std::vector<double>& b)
    {
        for (int i = 0; i < nnz; i++)
        {
            b[idx2[i]] -= vals[i] * x[idx1[i]];
        }
    }
    void residual(const std::vector<double>& x, const std::vector<double>& b,
            std::vector<double>& r)
    {
        for (int i = 0; i < n_rows; i++)
            r[i] = b[i];
     
        for (int i = 0; i < nnz; i++)
        {
            r[idx1[i]] -= vals[i] * x[idx2[i]];
        }
    }

    CSRMatrix* mult(const CSRMatrix* B);
    CSRMatrix* mult(const CSCMatrix* B);
    CSRMatrix* mult(const COOMatrix* B);
    CSRMatrix* mult_T(const CSRMatrix* A);
    CSRMatrix* mult_T(const CSCMatrix* A);
    CSRMatrix* mult_T(const COOMatrix* A);

    void mult_append(Vector& x, Vector& b);
    void mult_append_neg(Vector& x, Vector& b);
    void mult_append_T(Vector& x, Vector& b);
    void mult_append_neg_T(Vector& x, Vector& b);

	Matrix* ilu_k(int lof);
	//Matrix* ilu_levels();
	//Matrix* ilu_sparsity(Matrix* levls, int lof);
	//Matrix* ilu_symbolic(int lof);
	//std::vector<double> ilu_numeric(Matrix* levls);



    void add_block(int row, int col, std::vector<double>& values);

    format_t format()
    {
        return COO;
    }

    std::vector<int>& rows()
    {
        return idx1;
    }

    std::vector<int>& cols()
    {
        return idx2;
    }

    std::vector<double>& data()
    {
        return vals;

    }
};

/**************************************************************
 *****   CSRMatrix Class (Inherits from Matrix Base Class)
 **************************************************************
 ***** This class constructs a sparse matrix in CSR format.
 *****
 ***** Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (CSR)
 ***** sort()
 *****    Sorts the matrix.  Already in row-wise order, but sorts
 *****    the columns in each row.
 ***** add_value(int row, int col, double val)
 *****     TODO -- add this functionality
 ***** indptr()
 *****     Returns std::vector<int>& row pointer.  The ith element points to
 *****     the index of indices() corresponding to the first column to lie on 
 *****     row i.
 ***** indices()
 *****     Returns std::vector<int>& containing the cols corresponding
 *****     to each nonzero
 ***** data()
 *****     Returns std::vector<double>& containing the nonzero values
 **************************************************************/
  class CSRMatrix : public Matrix
  {

  public:

    /**************************************************************
    *****   CSRMatrix Class Constructor
    **************************************************************
    ***** Initializes an empty CSRMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** _nrows : int
    *****    Number of rows in Matrix
    ***** _ncols : int
    *****    Number of columns in Matrix
    ***** nnz_per_row : int
    *****    Prediction of (approximately) number of nonzeros 
    *****    per row, used in reserving space
    **************************************************************/
    CSRMatrix(int _nrows, int _ncols, int _nnz = 0): Matrix(_nrows, _ncols)
    {
        idx1.resize(_nrows + 1);
        if (_nnz)
        {
            idx2.reserve(_nnz);
            vals.reserve(_nnz);
        }
    }

    CSRMatrix(int _nrows, int _ncols, double* _data) : Matrix(_nrows, _ncols)
    {
        n_rows = _nrows;
        n_cols = _ncols;
        nnz = 0;

        int nnz_dense = n_rows*n_cols;

        idx1.resize(n_rows + 1);
        if (nnz_dense)
        {
            idx2.reserve(nnz_dense);
            vals.reserve(nnz_dense);
        }

        idx1[0] = 0;
        for (int i = 0; i < n_rows; i++)
        {
            for (int j = 0; j < n_cols; j++)
            {
                double val = _data[i*n_cols + j];
                if (fabs(val) > zero_tol)
                {
                    idx2.push_back(j);
                    vals.push_back(val);
                    nnz++;
                }
            }
            idx1[i+1] = nnz;
        }
    }

    CSRMatrix(int _nrows, int _ncols, std::vector<int>& rowptr, 
            std::vector<int>& cols, std::vector<double>& data) : Matrix(_nrows, _ncols)
    {
        nnz = cols.size();
        idx1.resize(n_rows+1);
        idx2.resize(nnz);
        vals.resize(nnz);

        std::copy(rowptr.begin(), rowptr.end(), idx1.begin());
        std::copy(cols.begin(), cols.end(), idx2.begin());
        std::copy(data.begin(), data.end(), vals.begin());
    }

    /**************************************************************
    *****   CSRMatrix Class Constructor
    **************************************************************
    ***** Constructs a CSRMatrix from a COOMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const COOMatrix*
    *****    COOMatrix A, from which to copy data
    **************************************************************/
    explicit CSRMatrix(const COOMatrix* A) 
    {
        copy(A);
    }

    /**************************************************************
    *****   CSRMatrix Class Constructor
    **************************************************************
    ***** Constructs a CSRMatrix from a CSCMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const CSCMatrix*
    *****    CSCMatrix A, from which to copy data
    **************************************************************/
    explicit CSRMatrix(const CSCMatrix* A)
    {
        copy(A);
    }

    /**************************************************************
    *****   CSRMatrix Class Constructor
    **************************************************************
    ***** Constructs a CSRMatrix from a CSRMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const CSRMatrix*
    *****    CSRMatrix A, from which to copy data
    **************************************************************/
    explicit CSRMatrix(const CSRMatrix* A) 
    {
        copy(A);
    }

    CSRMatrix()
    {
    }

    ~CSRMatrix()
    {

    }

    Matrix* transpose();

    void print();

    void copy(const COOMatrix* A);
    void copy(const CSRMatrix* A);
    void copy(const CSCMatrix* A);
    void copy(const BSRMatrix* A);

    // Converts matrix to a dense flattened vector
    std::vector<double> to_dense() const;

    void add_value(int row, int col, double value);
    void sort();
    void move_diag();
    void remove_duplicates();

    template <typename T, typename U> void mult(T& x, U& b)
    { 
        Matrix::mult(x, b);
    }
    template <typename T, typename U> void mult_T(T& x, U& b)
    { 
        Matrix::mult_T(x, b);
    }
    template <typename T, typename U> void mult_append(T& x, U& b)
    { 
        Matrix::mult_append(x, b);
    }
    template <typename T, typename U> void mult_append_T(T& x, U& b)
    { 
        Matrix::mult_append_T(x, b);
    }
    template <typename T, typename U> void mult_append_neg(T& x, U& b)
    { 
        Matrix::mult_append_neg(x, b);
    }
    template <typename T, typename U> void mult_append_neg_T(T& x, U& b)
    { 
        Matrix::mult_append_neg_T(x, b);
    }
    template <typename T, typename U, typename V> 
    void residual(const T& x, const U& b, V& r)
    { 
        Matrix::residual(x, b, r);
    }

    void mult(std::vector<double>& x, std::vector<double>& b)
    {
        for (int i = 0; i < n_rows; i++)
            b[i] = 0.0;
        mult_append(x, b);
    }
    void mult_T(std::vector<double>& x, std::vector<double>& b)

    {
        for (int i = 0; i < n_cols; i++)
            b[i] = 0.0;

        mult_append_T(x, b);    
    }
    void mult_append(std::vector<double>& x, std::vector<double>& b)
    { 
        int start, end;
        for (int i = 0; i < n_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                b[i] += vals[j] * x[idx2[j]];
            }
        }
    }
    void mult_append_T(std::vector<double>& x, std::vector<double>& b)
    {
        int start, end;
        for (int i = 0; i < n_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                b[idx2[j]] += vals[j] * x[i];
            }
        }
    }
    void mult_append_neg(std::vector<double>& x, std::vector<double>& b)
    {
        int start, end;
        for (int i = 0; i < n_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                b[i] -= vals[j] * x[idx2[j]];
            }
        }
    }
    void mult_append_neg_T(std::vector<double>& x, std::vector<double>& b)
    {
        int start, end;
        for (int i = 0; i < n_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                b[idx2[j]] -= vals[j] * x[i];
            }
        }
    }
    void residual(const std::vector<double>& x, const std::vector<double>& b, 
            std::vector<double>& r)
    {
        for (int i = 0; i < n_rows; i++)
            r[i] = b[i];
     
        int start, end;
        for (int i = 0; i < n_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                r[i] -= vals[j] * x[idx2[j]];
            }
        }
    }


    CSRMatrix* mult(const CSRMatrix* B);
    CSRMatrix* mult(const CSCMatrix* B);
    CSRMatrix* mult(const COOMatrix* B);
    CSRMatrix* mult_T(const CSCMatrix* A);
    CSRMatrix* mult_T(const CSRMatrix* A);
    CSRMatrix* mult_T(const COOMatrix* A);

    CSRMatrix* subtract(CSRMatrix* B);

    CSRMatrix* strength(double theta = 0.0);
    CSRMatrix* aggregate();
    CSRMatrix* fit_candidates(data_t* B, data_t* R, int num_candidates, 
            double tol = 1e-10);

	Matrix* ilu_k(int lof);
	CSRMatrix* ilu_levels();
	CSRMatrix* ilu_sparsity(CSRMatrix* levls, int lof);
	CSRMatrix* ilu_symbolic(int lof);
	std::vector<double> ilu_numeric(CSRMatrix* levls);
    
	void add_block(int row, int col, std::vector<double>& values);

    format_t format()
    {
        return CSR;
    }

    std::vector<int>& row_ptr()
    {
        return idx1;
    }

    std::vector<int>& cols()
    {
        return idx2;
    }

    std::vector<double>& data()
    {
        return vals;
    }
};

/**************************************************************
 *****   CSCMatrix Class (Inherits from Matrix Base Class)
 **************************************************************
 ***** This class constructs a sparse matrix in CSC format.
 *****
 ***** Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (CSC)
 ***** sort()
 *****    Sorts the matrix.  Already in col-wise order, but sorts
 *****    the rows in each column.
 ***** add_value(int row, int col, double val)
 *****     TODO -- add this functionality
 ***** indptr()
 *****     Returns std::vector<int>& column pointer.  The ith element points to
 *****     the index of indices() corresponding to the first row to lie on 
 *****     column i.
 ***** indices()
 *****     Returns std::vector<int>& containing the rows corresponding
 *****     to each nonzero
 ***** data()
 *****     Returns std::vector<double>& containing the nonzero values
 **************************************************************/
  class CSCMatrix : public Matrix
  {

  public:

    CSCMatrix(int _nrows, int _ncols, int _nnz = 0): Matrix(_nrows, _ncols)
    {
        idx1.resize(_ncols + 1);
        if (_nnz)
        {
            idx2.reserve(_nnz);
            vals.reserve(_nnz);
        }
        nnz = _nnz;
    }

    CSCMatrix(int _nrows, int _ncols, double* _data) : Matrix(_nrows, _ncols)
    {
        int nnz_dense = n_rows*n_cols;

        idx1.resize(n_cols + 1);
        if (nnz_dense)
        {
            idx2.reserve(nnz_dense);
            vals.reserve(nnz_dense);
        }

        idx1[0] = 0;
        for (int i = 0; i < n_cols; i++)
        {
            for (int j = 0; j < n_rows; j++)
            {
                double val = _data[i*n_cols + j];
                if (fabs(val) > zero_tol)
                {
                    idx2.push_back(j);
                    vals.push_back(val);
                    nnz++;
                }
            }
            idx1[i+1] = nnz;
        }
    }

    CSCMatrix(int _nrows, int _ncols, std::vector<int>& colptr, 
            std::vector<int>& rows, std::vector<double>& data) : Matrix(_nrows, _ncols)
    {
        nnz = rows.size();
        idx1.resize(n_cols+1);
        idx2.resize(nnz);
        vals.resize(nnz);

        std::copy(colptr.begin(), colptr.end(), idx1.begin());
        std::copy(rows.begin(), rows.end(), idx2.begin());
        std::copy(data.begin(), data.end(), vals.begin());
    }

    /**************************************************************
    *****   CSCMatrix Class Constructor
    **************************************************************
    ***** Constructs a CSCMatrix from a COOMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const COOMatrix*
    *****    COOMatrix A, from which to copy data
    **************************************************************/
    explicit CSCMatrix(const COOMatrix* A) 
    {
        copy(A);
    }

    /**************************************************************
    *****   CSCMatrix Class Constructor
    **************************************************************
    ***** Constructs a CSCMatrix from a CSRMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const CSRMatrix*
    *****    CSRMatrix A, from which to copy data
    **************************************************************/
    explicit CSCMatrix(const CSRMatrix* A) 
    {
        copy(A);
    }

    /**************************************************************
    *****   CSCMatrix Class Constructor
    **************************************************************
    ***** Constructs a CSCMatrix from a CSCMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const CSCMatrix*
    *****    CSCMatrix A, from which to copy data
    **************************************************************/
    explicit CSCMatrix(const CSCMatrix* A) 
    {
        copy(A);
    }

    CSCMatrix()
    {
    }

    ~CSCMatrix()
    {

    }


    Matrix* transpose();
    void print();

    void copy(const COOMatrix* A);
    void copy(const CSRMatrix* A);
    void copy(const CSCMatrix* A);
    void copy(const BSRMatrix* A);

    void sort();
    void move_diag();
    void remove_duplicates();
    void add_value(int row, int col, double value);

    template <typename T, typename U> void mult(T& x, U& b)
    { 
        Matrix::mult(x, b);
    }
    template <typename T, typename U> void mult_T(T& x, U& b)
    { 
        Matrix::mult_T(x, b);
    }
    template <typename T, typename U> void mult_append(T& x, U& b)
    { 
        Matrix::mult_append(x, b);
    }
    template <typename T, typename U> void mult_append_T(T& x, U& b)
    { 
        Matrix::mult_append_T(x, b);
    }
    template <typename T, typename U> void mult_append_neg(T& x, U& b)
    { 
        Matrix::mult_append_neg(x, b);
    }
    template <typename T, typename U> void mult_append_neg_T(T& x, U& b)
    { 
        Matrix::mult_append_neg_T(x, b);
    }
    template <typename T, typename U, typename V> 
    void residual(const T& x, const U& b, V& r)
    { 
        Matrix::residual(x, b, r);
    }

    void mult(std::vector<double>& x, std::vector<double>& b)
    {
        for (int i = 0; i < n_rows; i++)
            b[i] = 0.0;
        mult_append(x, b);
    }
    void mult_T(std::vector<double>& x, std::vector<double>& b)
    {
        for (int i = 0; i < n_cols; i++)
            b[i] = 0.0;

        mult_append_T(x, b);
    }
    void mult_append(std::vector<double>& x, std::vector<double>& b)
    { 
        int start, end;
        for (int i = 0; i < n_cols; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                b[idx2[j]] += vals[j] * x[i];
            }
        }
    }
    void mult_append_T(std::vector<double>& x, std::vector<double>& b)
    {
        int start, end;
        for (int i = 0; i < n_cols; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                b[i] += vals[j] * x[idx2[j]];
            }
        }
    }
    void mult_append_neg(std::vector<double>& x, std::vector<double>& b)
    {
        int start, end;
        for (int i = 0; i < n_cols; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                b[idx2[j]] -= vals[j] * x[i];
            }
        }
    }
    void mult_append_neg_T(std::vector<double>& x, std::vector<double>& b)
    {
        int start, end;
        for (int i = 0; i < n_cols; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                b[i] -= vals[j] * x[idx2[j]];
            }
        }
    }
    void residual(const std::vector<double>& x, const std::vector<double>& b, 
            std::vector<double>& r)
    {
        for (int i = 0; i < n_rows; i++)
            r[i] = b[i];

        int start, end;
        for (int i = 0; i < n_cols; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                r[idx2[j]] -= vals[j] * x[i];
            }
        }
    }

    CSRMatrix* mult(const CSRMatrix* B);
    CSRMatrix* mult(const CSCMatrix* B);
    CSRMatrix* mult(const COOMatrix* B);
    CSRMatrix* mult_T(const CSRMatrix* A);
    CSRMatrix* mult_T(const CSCMatrix* A);
    CSRMatrix* mult_T(const COOMatrix* A);

    void add_block(int row, int col, std::vector<double>& values);

    void jacobi(Vector& x, Vector& b, Vector& tmp, double omega = .667);    

    format_t format()
    {
        return CSC;
    }

    std::vector<int>& col_ptr()
    {
        return idx1;
    }

    std::vector<int>& rows()
    {
        return idx2;
    }

    std::vector<double>& data()
    {
        return vals;
    }
	
	Matrix* ilu_k(int lof);
	//Matrix* ilu_levels();
	//Matrix* ilu_sparsity(Matrix* levls, int lof);
	//Matrix* ilu_symbolic(int lof);
	//std::vector<double> ilu_numeric(Matrix* levls);

  };


/**************************************************************
 *****   BSRMatrix Class (Inherits from Matrix Base Class)
 **************************************************************
 ***** This class constructs a sparse matrix in BSR format.
 *****
 ***** Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (BSR)
 ***** add_value(int row, int col, double val)
 *****     TODO -- add this functionality
 ***** add_block(int row, int col, std::vector<double>& data)
 *****     Adds the row-wise flattened block 'data' to the matrix
 *****     at block location (row, col) in the coarse matrix defined 
 *****     by blocks - NOT global row and column indices
 ***** row_ptr()
 *****     Returns std::vector<int>& row pointer.  The ith element points to
 *****     the index of indices() corresponding to the first block column to 
 *****     lie on block row i.
 ***** cols()
 *****     Returns std::vector<int>& containing the cols corresponding
 *****     to each nonzero dense block
 ***** data()
 *****     Returns std::vector<double>& containing the nonzero values
 *****     - flattened array of block values 
 ***** block_rows()
 *****     Returns b_rows - number of rows per block
 ***** block_cols()
 *****     Returns b_cols - number of columns per block
 ***** block_size()
 *****     Returns nnz in dense block
 ***** num_blocks()
 *****     Returns number of dense blocks in sparse matrix
 **************************************************************/
  class BSRMatrix : public Matrix
  {
  public:

    /**************************************************************
    *****   BSRMatrix Class Constructor
    **************************************************************
    ***** Initializes an empty BSRMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** _nrows : int
    *****    Number of rows in Matrix
    ***** _ncols : int
    *****    Number of columns in Matrix
    ***** _brows : int
    *****    Number of rows in block
    ***** _bcols : int
    *****    Number of columns in block
    ***** _nblocks : int
    *****    Number of nonzero blocks in matrix
    ***** nnz_per_row : int
    *****    Prediction of (approximately) number of nonzeros 
    *****    per row, used in reserving space
    *****
    ***** idx2 : columns of each block in matrix (row-ordered)
    ***** idx1 : block row pointer - first index in idx1 of block in row i
    ***** b_rows : rows per block
    ***** b_cols : columns per block
    ***** n_blocks : number of dense blocks in matrix
    ***** b_size : number of non-zeros in a block
    **************************************************************/
    BSRMatrix(int _nrows, int _ncols, int _brows, int _bcols, 
            int _nblocks=0, int _nnz = 0): Matrix(_nrows, _ncols)
    {
        if (_nrows % _brows != 0 || _ncols % _bcols != 0)
	{
            printf("Matrix dimensions must be divisible by block dimensions.\n");
	    exit(-1);
	}

	n_rows = _nrows;
	n_cols = _ncols;
	b_rows = _brows;
	b_cols = _bcols;
	b_size = b_rows * b_cols;

	if (_nblocks)
	{
	    n_blocks = _nblocks;
	}
	else if (_brows != 0 && _bcols != 0)
	{
	    // Assume dense number of blocks
            n_blocks = _nrows/_brows * _ncols/_bcols;
	}

        idx1.resize(n_rows/b_rows + 1);
        idx2.reserve(n_blocks);
        vals.reserve(b_size * n_blocks);
    }

    // Constructs BSRMatrix from flattened _data array of entire matrix 
    // - dropping blocks that are entirely zero
    // Assumes data array is flattened array of matrix in 'block' format
    BSRMatrix(int _nrows, int _ncols, int _brows, int _bcols, double* _data) : Matrix(_nrows, _ncols)
    {
        if (_nrows % _brows != 0 || _ncols % _bcols != 0)
	{
            printf("Matrix dimensions must be divisible by block dimensions.\n");
	    exit(-1);
	}

	// Assumes dense data array given
	n_rows = _nrows;
	n_cols = _ncols;
	b_rows = _brows;
	b_cols = _bcols;
	b_size = b_rows * b_cols;
	n_blocks = 0;
	nnz = 0;

        int nnz_dense = n_rows*n_cols;

        idx1.resize(n_rows/b_rows + 1);
        //idx2.reserve(n_blocks);
        //vals.reserve(nnz_dense);

	std::vector<double> test;
	double val;
	int data_offset = 0;
	idx1[0] = 0;
	for (int i=0; i<n_rows/b_rows; i++)
	{
            for (int j=0; j<n_cols/b_cols; j++)
	    {
		// 1. Push block data to test vector & check if it's a 0 block
                for (int k=data_offset; k<data_offset+b_size; k++){
		    val = _data[k];
		    if (fabs(val) > zero_tol){
		         test.push_back(val);
	            }
		}

		// 2. If not all 0 then add block
		if (test.size() > 0)
		{
		    for (int k=data_offset; k<data_offset+b_size; k++){
			val = _data[k];
		        vals.push_back(val);
			nnz++;
		    }
		    n_blocks++;
		    idx2.push_back(j);
		}
                data_offset += b_size;
		test.clear();
	    }
	    idx1[i+1] = idx2.size();
	}

    }

    // Constructs BSRMatrix of size _nrows * _ncols with blocks of size _brows * _bcols
    // and rowptr, cols, and data vectors given
    BSRMatrix(int _nrows, int _ncols, int _brows, int _bcols, 
            std::vector<int>& rowptr, std::vector<int>& cols, 
	    std::vector<double>& data) : Matrix(_nrows, _ncols)
    {
        if (_nrows % _brows != 0 || _ncols % _bcols != 0)
	{
            printf("Matrix dimensions must be divisible by block dimensions.\n");
	    exit(-1);
	}

        nnz = data.size();
	n_rows = _nrows;
	n_cols = _ncols;
	b_rows = _brows;
	b_cols = _bcols;
	n_blocks = cols.size();
	b_size = nnz/n_blocks;
        idx1.resize(n_rows/b_rows + 1);
        idx2.resize(n_blocks);
        vals.resize(nnz);

        std::copy(rowptr.begin(), rowptr.end(), idx1.begin());
        std::copy(cols.begin(), cols.end(), idx2.begin());
        std::copy(data.begin(), data.end(), vals.begin());
    }

    /**************************************************************
    *****   BSRMatrix Class Constructor
    **************************************************************
    ***** Constructs a BSRMatrix from a COOMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const COOMatrix*
    *****    COOMatrix A, from which to copy data
    **************************************************************/
    explicit BSRMatrix(const COOMatrix* A, int _brows, int _bcols) 
    {
	b_rows = _brows;
	b_cols = _bcols;
	b_size = b_rows * b_cols;

	copy(A);
    }

    /**************************************************************
    *****   BSRMatrix Class Constructor
    **************************************************************
    ***** Constructs a BSRMatrix from a CSCMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const CSCMatrix*
    *****    CSCMatrix A, from which to copy data
    **************************************************************/
    /*explicit BSRMatrix(const CSCMatrix* A)
    {
        copy(A);
    }*/

    /**************************************************************
    *****   BSRMatrix Class Constructor
    **************************************************************
    ***** Constructs a BSRMatrix from a CSRMatrix
    *****
    ***** Parameters
    ***** -------------
    ***** A : const CSRMatrix*
    *****    CSRMatrix A, from which to copy data
    **************************************************************/
    explicit BSRMatrix(const CSRMatrix* A, int _brows, int _bcols) 
    {
        b_rows = _brows;
    	b_cols = _bcols;
	    b_size = b_rows * b_cols;

        copy(A);
    }

    BSRMatrix()
    {
    }

    ~BSRMatrix()
    {

    }

    Matrix* transpose();

    void print();
    void block_print(int row, int num_blocks_prev, int col);
    std::vector<double> to_dense();

    void copy(const COOMatrix* A);
    void copy(const CSRMatrix* A);
    void copy(const CSCMatrix* A);
    void copy(const BSRMatrix* A);

    void add_value(int row, int col, double value);
    void add_block(int row, int col, std::vector<double>& values);
    void sort();
    void move_diag();
    void remove_duplicates();

    template <typename T, typename U> void mult(T& x, U& b)
    { 
        Matrix::mult(x, b);
    }
    template <typename T, typename U> void mult_T(T& x, U& b)
    { 
        Matrix::mult_T(x, b);
    }
    template <typename T, typename U> void mult_append(T& x, U& b)
    { 
        Matrix::mult_append(x, b);
    }
    template <typename T, typename U> void mult_append_T(T& x, U& b)
    { 
        Matrix::mult_append_T(x, b);
    }
    template <typename T, typename U> void mult_append_neg(T& x, U& b)
    { 
        Matrix::mult_append_neg(x, b);
    }
    template <typename T, typename U> void mult_append_neg_T(T& x, U& b)
    { 
        Matrix::mult_append_neg_T(x, b);
    }
    template <typename T, typename U, typename V> 
    void residual(const T& x, const U& b, V& r)
    { 
        Matrix::residual(x, b, r);
    }

    // STANDARD MULTIPLICATION
    void mult(std::vector<double>& x, std::vector<double>& b)
    {
        for (int i = 0; i < n_rows; i++)
            b[i] = 0.0;
        mult_append(x, b);
    }

    // TRANSPOSE MULTIPLICATION
    void mult_T(std::vector<double>& x, std::vector<double>& b)

    {
        for (int i = 0; i < n_cols; i++)
            b[i] = 0.0;

        mult_append_T(x, b);    
    }

    // STANDARD MULTIPLICATION HELPER
    void mult_append(std::vector<double>& x, std::vector<double>& b)
    { 
        int start, end;
	int rowsOfBlocks = n_rows/b_rows;
        for (int i = 0; i < rowsOfBlocks; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
		// Dense multiplication on block
                block_mult(i, j, idx2[j], x, b);
            }
        }
    }

    void block_mult(int row, int num_blocks_prev, int col,
		    std::vector<double>& x, std::vector<double>& b)
    {
        int upper_i = row * b_rows;
	int upper_j = col * b_cols;
	int data_offset = num_blocks_prev * b_size;

	int glob_i, glob_j, ind;
	for(int i=0; i<b_rows; i++){
            for(int j=0; j<b_cols; j++){
		glob_i = upper_i + i;
		glob_j = upper_j + j;
		ind = i * b_cols + j + data_offset;
		b[glob_i] += vals[ind] * x[glob_j];
	    }
	}
    }

    void mult_append_T(std::vector<double>& x, std::vector<double>& b)
    {
        int start, end;
        for (int i = 0; i < n_rows/b_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                // Dense transpose multiplication on block
                block_mult_T(i, j, idx2[j], x, b);
            }
        }
    }

    void block_mult_T(int row, int num_blocks_prev, int col,
		    std::vector<double>& x, std::vector<double>& b)
    {
        int upper_i = row * b_rows;
	int upper_j = col * b_cols;
	int data_offset = num_blocks_prev * b_size;

	int glob_i, glob_j, ind;
	for(int i=0; i<b_rows; i++){
            for(int j=0; j<b_cols; j++){
		glob_i = upper_i + i;
		glob_j = upper_j + j;
		ind = i * b_cols + j + data_offset;
		b[glob_j] += vals[ind] * x[glob_i];
	    }
	}
    }

    void mult_append_neg(std::vector<double>& x, std::vector<double>& b)
    {
        int start, end;
        for (int i = 0; i < n_rows/b_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                // Dense negative multiplication on block
                block_mult_neg(i, j, idx2[j], x, b);
            }
        }
    }

    void block_mult_neg(int row, int num_blocks_prev, int col,
		    std::vector<double>& x, std::vector<double>& b)
    {
        int upper_i = row * b_rows;
	int upper_j = col * b_cols;
	int data_offset = num_blocks_prev * b_size;

	int glob_i, glob_j, ind;
	for(int i=0; i<b_rows; i++){
            for(int j=0; j<b_cols; j++){
		glob_i = upper_i + i;
		glob_j = upper_j + j;
		ind = i * b_cols + j + data_offset;
		b[glob_i] -= vals[ind] * x[glob_j];
	    }
	}
    }

    void mult_append_neg_T(std::vector<double>& x, std::vector<double>& b)
    {
        int start, end;
        for (int i = 0; i < n_rows/b_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                block_mult_neg_T(i, j, idx2[j], x, b);
            }
        }
    }

    void block_mult_neg_T(int row, int num_blocks_prev, int col,
		    std::vector<double>& x, std::vector<double>& b)
    {
        int upper_i = row * b_rows;
	int upper_j = col * b_cols;
	int data_offset = num_blocks_prev * b_size;

	int glob_i, glob_j, ind;
	for(int i=0; i<b_rows; i++){
            for(int j=0; j<b_cols; j++){
		glob_i = upper_i + i;
		glob_j = upper_j + j;
		ind = i * b_cols + j + data_offset;
		b[glob_j] -= vals[ind] * x[glob_i];
	    }
	}
    }

    void residual(const std::vector<double>& x, const std::vector<double>& b, 
            std::vector<double>& r)
    {
        for (int i = 0; i < n_rows; i++)
            r[i] = b[i];

        int start, end;
        for (int i = 0; i < n_rows/b_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                block_res(i, j, idx2[j], x, r);
            }
        }
    }

    void block_res(int row, int num_blocks_prev, int col,
		    const std::vector<double>& x, std::vector<double>& r)
    {
        int upper_i = row * b_rows;
	int upper_j = col * b_cols;
	int data_offset = num_blocks_prev * b_size;

	int glob_i, glob_j, ind;
	for(int i=0; i<b_rows; i++){
            for(int j=0; j<b_cols; j++){
		glob_i = upper_i + i;
		glob_j = upper_j + j;
		ind = i * b_cols + j + data_offset;
		r[glob_i] -= vals[ind] * x[glob_j];
	    }
	}
    }


    //CSRMatrix* mult(const CSRMatrix* B);
    //CSRMatrix* mult(const CSCMatrix* B);
    //CSRMatrix* mult(const COOMatrix* B);
    //CSRMatrix* mult_T(const CSCMatrix* A);
    //CSRMatrix* mult_T(const CSRMatrix* A);
    //CSRMatrix* mult_T(const COOMatrix* A);

    //CSRMatrix* subtract(CSRMatrix* B);

    //CSRMatrix* strength(double theta = 0.0);
    //CSRMatrix* aggregate();
    //CSRMatrix* fit_candidates(data_t* B, data_t* R, int num_candidates, 
      //      double tol = 1e-10);

    format_t format()
    {
        return BSR;
    }

    std::vector<int>& row_ptr()
    {
        return idx1;
    }

    std::vector<int>& cols()
    {
        return idx2;
    }

    std::vector<double>& data()
    {
        return vals;
    }

    int block_rows()
    {
        return b_rows;
    }

    int block_cols()
    {
        return b_cols;
    }

    int block_size()
    {
        return b_size;
    }

    int num_blocks()
    {
        return n_blocks;
    }

    int b_rows;
    int b_cols;
    int n_blocks;
    int b_size;

//	Matrix* ilu_k(int lof);
	BSRMatrix* ilu_k(int lof);
	BSRMatrix* ilu_levels();
	BSRMatrix* ilu_sparsity(BSRMatrix* levls, int lof);
	BSRMatrix* ilu_symbolic(int lof);
	std::vector<double> ilu_numeric(BSRMatrix* levls);

	std::vector<double> get_diag(BSRMatrix* levls, std::vector<double> data);
	std::vector<double> fill_factors(BSRMatrix* levls);
	std::vector<double> inv_diag_block(std::vector<double> diag_vec, int k);
	std::vector<double> mult_b(std::vector<double> block_a, std::vector<double> block_b);
};

}


// Copyright (c) 2015-2017, RAPtor Developer Team
// License: Simplified BSD, http://opensource.org/licenses/BSD-2-Clause


using namespace raptor;

/**************************************************************
*****  Matrix Print
**************************************************************
***** Print the nonzeros in the matrix, as well as the row
***** and column according to each nonzero
**************************************************************/
void COOMatrix::print()
{
    int row, col;
    double val;

    for (int i = 0; i < nnz; i++)
    {
        row = idx1[i];
        col = idx2[i];
        val = vals[i];

        printf("A[%d][%d] = %e\n", row, col, val);
    }
}
void CSRMatrix::print()
{
    int col, start, end;
    double val;

    for (int row = 0; row < n_rows; row++)
    {
        start = idx1[row];
        end = idx1[row+1];
        for (int j = start; j < end; j++)
        {
            col = idx2[j];
            val = vals[j];

            printf("A[%d][%d] = %e\n", row, col, val);
        }
    }
}
void CSCMatrix::print()
{
    int row, start, end;
    double val;

    for (int col = 0; col < n_cols; col++)
    {
        start = idx1[col];
        end = idx1[col+1];
        for (int j = start; j < end; j++)
        {
            row = idx2[j];
            val = vals[j];

            printf("A[%d][%d] = %e\n", row, col, val);
        }
    }
}

void BSRMatrix::print()
{
    printf("void BSRMatrix::print()\n");
    int col, start, end;

//    print("n_rows=%d\tb_rows=%d\n", n_rows, b_rows);
    for (int i = 0; i < n_rows/b_rows; i++)
    {
        start = idx1[i];
        end = idx1[i+1];
        for (int j = start; j < end; j++)
        {
            // Call block print function
	    block_print(i, j, idx2[j]);
	    printf("----------\n");
        }
    }
}

void BSRMatrix::block_print(int row, int num_blocks_prev, int col)
{
    // The upper left corner indices of the matrix
    int upper_i = row * b_rows;
    int upper_j = col * b_cols;
    // The offset to determine where this block starts in the data array
    int data_offset = num_blocks_prev * b_size;

    int glob_i, glob_j, ind;
    double val;
    for (int i=0; i<b_rows; i++)
    {
        for (int j=0; j<b_cols; j++)
	{
            glob_i = upper_i + i;
	    glob_j = upper_j + j;
	    ind = i * b_cols + j + data_offset;
	    val = vals[ind];
	    printf("A[%d][%d] = %e\n", glob_i, glob_j, val);
	}
    }
}


Matrix* COOMatrix::transpose()
{
    Matrix* T = new COOMatrix(n_rows, n_cols, idx2, idx1, vals);

    return T;
}

Matrix* CSRMatrix::transpose()
{
    // Create CSC Matrix... rowptr is now colptr
    CSCMatrix* T_csc = new CSCMatrix(n_rows, n_cols, idx1, idx2, vals); 

    // Convert back to CSR to tranpose
    Matrix* T = new CSRMatrix(T_csc);

    delete T_csc;

    return T;
}

Matrix* CSCMatrix::transpose()
{
    // Create CSR Matrix... colptr is now rowptr
    CSRMatrix* T_csr = new CSRMatrix(n_rows, n_cols, idx1, idx2, vals); 

    // Convert back to CSC to tranpose
    Matrix* T = new CSCMatrix(T_csr);

    delete T_csr;

    return T;
}

Matrix* BSRMatrix::transpose()
{
    printf("Currently not implemented.\n");	
    return NULL;
}

/**************************************************************
*****   Matrix Resize
**************************************************************
***** Set the matrix dimensions to those passed as parameters
*****
***** Parameters
***** -------------
***** _nrows : int
*****    Number of rows in matrix
***** _ncols : int
*****    Number of cols in matrix
**************************************************************/
void Matrix::resize(int _n_rows, int _n_cols)
{
    n_rows = _n_rows;
    n_cols = _n_cols;
}


/**************************************************************
*****  COOMatrix Add Value
**************************************************************
***** Inserts value into the position (row, col) of the matrix
*****
***** Parameters
***** -------------
***** row : int
*****    Row in which to insert value 
***** col : int
*****    Column in which to insert value
***** value : double
*****    Nonzero value to be inserted into the matrix
**************************************************************/
void COOMatrix::add_value(int row, int col, double value)
{
    idx1.push_back(row);
    idx2.push_back(col);
    vals.push_back(value);
    nnz++;
}

void COOMatrix::add_block(int row, int col, std::vector<double>& values){
    printf("Not implemented.\n");
}

void COOMatrix::copy(const COOMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = A->nnz;

    idx1.clear();
    idx2.clear();
    vals.clear();

    idx1.reserve(A->nnz);
    idx2.reserve(A->nnz);
    vals.reserve(A->nnz);
    for (int i = 0; i < A->nnz; i++)
    {
        idx1.push_back(A->idx1[i]);
        idx2.push_back(A->idx2[i]);
        vals.push_back(A->vals[i]);
    }
}
void COOMatrix::copy(const CSRMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = A->nnz;

    idx1.clear();
    idx2.clear();
    vals.clear();

    idx1.reserve(A->nnz);
    idx2.reserve(A->nnz);
    vals.reserve(A->nnz);
    for (int i = 0; i < A->n_rows; i++)
    {
        int row_start = A->idx1[i];
        int row_end = A->idx1[i+1];
        for (int j = row_start; j < row_end; j++)
        {
            idx1.push_back(i);
            idx2.push_back(A->idx2[j]);
            vals.push_back(A->vals[j]);
        }
    }
}
void COOMatrix::copy(const CSCMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = A->nnz;

    idx1.clear();
    idx2.clear();
    vals.clear();

    idx1.reserve(A->nnz);
    idx2.reserve(A->nnz);
    vals.reserve(A->nnz);
    for (int i = 0; i < A->n_cols; i++)
    {
        int col_start = A->idx1[i];
        int col_end = A->idx1[i+1];
        for (int j = col_start; j < col_end; j++)
        {
            idx1.push_back(A->idx2[j]);
            idx2.push_back(i);
            vals.push_back(A->vals[j]);
        }
    }
}

void COOMatrix::copy(const BSRMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = 0;

    idx1.clear();
    idx2.clear();
    vals.clear();

    idx1.reserve(A->nnz);
    idx2.reserve(A->nnz);
    vals.reserve(A->nnz);

    for (int i = 0; i < n_rows/A->b_rows; i++)
    {
        int row_start = A->idx1[i];
        int row_end = A->idx1[i+1];
        for (int j = row_start; j < row_end; j++)
        {
            // Call block copy function
	    block_copy(A, i, j, A->idx2[j]);
        }
    }
}

void COOMatrix::block_copy(const BSRMatrix* A, int row, int num_blocks_prev, int col)
{
    int upper_i = row * A->b_rows;
    int upper_j = col * A->b_cols;
    int data_offset = num_blocks_prev * A->b_size;

    int glob_i, glob_j, ind;
    double val;
    for (int i = 0; i < A->b_rows; i++)
    {
        for (int j = 0; j < A->b_cols; j++)
	{
            glob_i = upper_i + i;
	    glob_j = upper_j + j;
	    ind = i * A->b_cols + j + data_offset;
	    val = A->vals[ind];

	    if (fabs(val) > zero_tol)
	    {
	        idx1.push_back(glob_i);
	        idx2.push_back(glob_j);
                vals.push_back(val);
		nnz++;
	    }
	}
    }
}

/**************************************************************
*****   COOMatrix to_dense
**************************************************************
***** Converts the COOMatrix into a dense matrix
***** in the form of a flattened vector ordered row-wise
**************************************************************/
std::vector<double> COOMatrix::to_dense() const
{
    std::vector<double> dense(n_rows * n_cols);
    std::fill(dense.begin(), dense.end(), 0.0);

    for (int i = 0; i < nnz; i++)
    {
        dense[idx1[i]*n_cols + idx2[i]] = vals[i];
    }

    return dense;
}

/**************************************************************
*****   COOMatrix Sort
**************************************************************
***** Sorts the sparse matrix by row, and by column within 
***** each row.  Removes duplicates, summing their values 
***** together.
**************************************************************/
void COOMatrix::sort()
{
    if (sorted || nnz == 0)
    {
        sorted = true;
        return;
    }

    int k, prev_k;

    std::vector<int> permutation(nnz);
    std::vector<bool> done(nnz, false);

    // Create permutation vector p
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(permutation.begin(), permutation.end(),
        [&](int i, int j){ 
            if (idx1[i] == idx1[j])
                return idx2[i] < idx2[j];
            else
                return idx1[i] < idx1[j];
        });

    // Permute all vectors (rows, cols, data) 
    // according to p
    for (int i = 0; i < nnz; i++)
    {
        if (done[i]) continue;

        done[i] = true;
        prev_k = i;
        k = permutation[i];
        while (i != k)
        {
            std::swap(idx1[prev_k], idx1[k]);
            std::swap(idx2[prev_k], idx2[k]);
            std::swap(vals[prev_k], vals[k]);
            done[k] = true;
            prev_k = k;
            k = permutation[k];
        }
    }

    sorted = true;
    diag_first = false;
}

void COOMatrix::move_diag()
{
    if (diag_first || nnz == 0)
    {
        return;
    }

    if (!sorted)
    {
        sort();
    }

    int row_start, prev_row;
    int row, col;
    double tmp;

    // Move diagonal entry to first in row
    row_start = 0;
    prev_row = 0;
    for (int i = 0; i < nnz; i++)
    {
        row = idx1[i];
        col = idx2[i];
        if (row != prev_row)
        {
            prev_row = row;
            row_start = i;
        }
        else if (row == col)
        {
            tmp = vals[i];
            for (int j = i; j > row_start; j--)
            {
                idx2[j] = idx2[j-1];
                vals[j] = vals[j-1];
            }
            idx2[row_start] = row;
            vals[row_start] = tmp;
        }
    }

    diag_first = true;
}

void COOMatrix::remove_duplicates()
{
    if (!sorted)
    {
        sort();
        diag_first = false;
    }

    int prev_row, prev_col, ctr;
    int row, col;
    double val;

    // Remove duplicates (sum together)
    prev_row = idx1[0];
    prev_col = idx2[0];
    ctr = 1;
    for (int i = 1; i < nnz; i++)
    {
        row = idx1[i];
        col = idx2[i];
        val = vals[i];
        if (row == prev_row && col == prev_col)
        {
            vals[ctr-1] += val;
        }
        else
        { 
            if (ctr != i)
            {
                idx1[ctr] = row;
                idx2[ctr] = col;
                vals[ctr] = val;
            }
            ctr++;

            prev_row = row;
            prev_col = col;
        }
    }

    nnz = ctr;
}

/**************************************************************
*****  CSRMatrix Add Value
**************************************************************
***** Inserts value into the position (row, col) of the matrix.
***** Values must be inserted in row-wise order, so if the row
***** is not equal to the row of the previously inserted value,
***** indptr is edited, and it is assumed that row is complete.
***** TODO -- this method needs further testing
*****
***** Parameters
***** -------------
***** row : int
*****    Row in which to insert value 
***** col : int
*****    Column in which to insert value
***** value : double
*****    Nonzero value to be inserted into the matrix
**************************************************************/
void CSRMatrix::add_value(int row, int col, double value)
{
    // Assumes idx1 is created separately
    idx2.push_back(col);
    vals.push_back(value);
    nnz++;
}

void CSRMatrix::add_block(int row, int col, std::vector<double>& values){
    printf("Not implemented.\n");
}

void CSRMatrix::copy(const COOMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = A->nnz;

    idx1.resize(n_rows + 1);
    std::fill(idx1.begin(), idx1.end(), 0);
    if (nnz)
    {
        idx2.resize(nnz);
        vals.resize(nnz);
    }

    // Calculate indptr
    for (int i = 0; i < nnz; i++)
    {
        int row = A->idx1[i];
        idx1[row+1]++;
    }
    for (int i = 0; i < n_rows; i++)
    {
        idx1[i+1] += idx1[i];
    }

    // Add indices and data
    std::vector<int> ctr;
    if (n_rows)
    {
    	ctr.resize(n_rows, 0);
    }
    for (int i = 0; i < nnz; i++)
    {
        int row = A->idx1[i];
        int col = A->idx2[i];
        double val = A->vals[i];
        int index = idx1[row] + ctr[row]++;
        idx2[index] = col;
        vals[index] = val;
    }
}
void CSRMatrix::copy(const CSRMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = A->nnz;

    idx1.resize(A->n_rows + 1);
    idx2.resize(A->nnz);
    vals.resize(A->nnz);

    idx1[0] = 0;
    for (int i = 0; i < A->n_rows; i++)
    {
        idx1[i+1] = A->idx1[i+1];
        int row_start = idx1[i];
        int row_end = idx1[i+1];
        for (int j = row_start; j < row_end; j++)
        {
            idx2[j] = A->idx2[j];
            vals[j] = A->vals[j];
        }
    }
}
void CSRMatrix::copy(const CSCMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = A->nnz;

    idx1.clear();
    idx2.clear();
    vals.clear();

    // Resize vectors to appropriate dimensions
    idx1.resize(A->n_rows + 1);
    idx2.resize(A->nnz);
    vals.resize(A->nnz);

    // Create indptr, summing number times row appears in CSC
    for (int i = 0; i <= A->n_rows; i++) idx1[i] = 0;
    for (int i = 0; i < A->nnz; i++)
    {
        idx1[A->idx2[i] + 1]++;
    }
    for (int i = 1; i <= A->n_rows; i++)
    {
        idx1[i] += idx1[i-1];
    }

    // Add values to indices and data
    std::vector<int> ctr(n_rows, 0);
    for (int i = 0; i < A->n_cols; i++)
    {
        int col_start = A->idx1[i];
        int col_end = A->idx1[i+1];
        for (int j = col_start; j < col_end; j++)
        {
            int row = A->idx2[j];
            int idx = idx1[row] + ctr[row]++;
            idx2[idx] = i;
            vals[idx] = A->vals[j];
        }
    }
}
void CSRMatrix::copy(const BSRMatrix* A)
{
    printf("Currently not implemented\n");
}

/**************************************************************
*****   CSRMatrix to_dense
**************************************************************
***** Converts the CSRMatrix into a dense matrix
***** in the form of a flattened vector ordered row-wise
**************************************************************/
std::vector<double> CSRMatrix::to_dense() const
{
    std::vector<double> dense(n_rows * n_cols);
    std::fill(dense.begin(), dense.end(), 0.0);

    for (int i = 0; i < n_rows; i++)
    {
        for (int j = idx1[i]; j < idx1[i+1]; j++)
	{
            dense[i*n_cols + idx2[j]] = vals[j];
	}
    }
    return dense;
}

/**************************************************************
*****   CSRMatrix Sort
**************************************************************
***** Sorts the sparse matrix by columns within each row.  
***** Removes duplicates, summing their values 
***** together.
**************************************************************/
void CSRMatrix::sort()
{
    int start, end, row_size;
    int k, prev_k;

    if (sorted || nnz == 0)
    {
        sorted = true;
        return;
    }

    std::vector<int> permutation;
    std::vector<bool> done;

    // Sort the columns of each row (and data accordingly) and remove
    // duplicates (summing values together)
    for (int row = 0; row < n_rows; row++)
    {
        start = idx1[row];
        end = idx1[row+1];
        row_size = end - start;
        if (row_size == 0) 
        {
            continue;
        }

        // Create permutation vector p for row
        permutation.resize(row_size);
        std::iota(permutation.begin(), permutation.end(), 0);
        std::sort(permutation.begin(), permutation.end(),
                [&](int i, int j)
                { 
                    return idx2[i+start] < idx2[j+start];
                });


        // Permute columns and data according to p
        done.resize(row_size);
        for (int i = 0; i < row_size; i++)
        {
            done[i] = false;
        }
        if (vals.size())
        {
            for (int i = 0; i < row_size; i++)
            {
                if (done[i]) continue;

                done[i] = true;
                prev_k = i;
                k = permutation[i];
                while (i != k)
                {
                    std::swap(idx2[prev_k + start], idx2[k + start]);
                    std::swap(vals[prev_k + start], vals[k + start]);
                    done[k] = true;
                    prev_k = k;
                    k = permutation[k];
                }
            }
        }
        else
        {
            for (int i = 0; i < row_size; i++)
            {
                if (done[i]) continue;

                done[i] = true;
                prev_k = i;
                k = permutation[i];
                while (i != k)
                {
                    std::swap(idx2[prev_k + start], idx2[k + start]);
                    done[k] = true;
                    prev_k = k;
                    k = permutation[k];
                }
            }
        }
    }

    sorted = true;
    diag_first = false;
}

void CSRMatrix::move_diag()
{
    int start, end;
    int col;
    double tmp;

    if (diag_first || nnz == 0)
    {
        return;
    }

    // Move diagonal values to beginning of each row
    if (vals.size())
    {
        for (int i = 0; i < n_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                col = idx2[j];
                if (col == i)
                {
                    tmp = vals[j];
                    for (int k = j; k > start; k--)
                    {
                        idx2[k] = idx2[k-1];
                        vals[k] = vals[k-1];
                    }
                    idx2[start] = i;
                    vals[start] = tmp;
                    break;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < n_rows; i++)
        {
            start = idx1[i];
            end = idx1[i+1];
            for (int j = start; j < end; j++)
            {
                col = idx2[j];
                if (col == i)
                {
                    for (int k = j; k > start; k--)
                    {
                        idx2[k] = idx2[k-1];
                    }
                    idx2[start] = i;
                    break;
                }
            }
        }
    }
    diag_first = true;
}

void CSRMatrix::remove_duplicates()
{
    int orig_start, orig_end;
    int new_start;
    int col, prev_col;
    int ctr, row_size;
    double val;

    if (!sorted)
    {
        sort();
        diag_first = false;
    }

    orig_start = idx1[0];
    for (int row = 0; row < n_rows; row++)
    {
        new_start = idx1[row];
        orig_end = idx1[row+1];
        row_size = orig_end - orig_start;
        if (row_size == 0) 
        {
            orig_start = orig_end;
            idx1[row+1] = idx1[row];
            continue;
        }

        // Remove Duplicates
        col = idx2[orig_start];
        val = vals[orig_start];
        idx2[new_start] = col;
        vals[new_start] = val;
        prev_col = col;
        ctr = 1;
        for (int j = orig_start + 1; j < orig_end; j++)
        {
            col = idx2[j];
            val = vals[j];
            if (col == prev_col)
            {
                vals[ctr - 1 + new_start] += val;
            }
            else
            {
                if (fabs(vals[ctr - 1 + new_start]) < zero_tol)
                {
                    ctr--;
                }

                idx2[ctr + new_start] = col;
                vals[ctr + new_start] = val;
                ctr++;
                prev_col = col;
            }
        }
        if (fabs(vals[ctr - 1 + new_start]) < zero_tol)
        {
            ctr--;
        }

        orig_start = orig_end;
        idx1[row+1] = idx1[row] + ctr;
    }
    nnz = idx1[n_rows];
    idx2.resize(nnz);
    vals.resize(nnz);
}

/**************************************************************
*****  BSRMatrix Copy
**************************************************************
**************************************************************/
void BSRMatrix::copy(const COOMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;

    std::vector<double> A_dense = A->to_dense();

    double block_dense[n_rows*n_cols];

    int block_ind, glob_i, glob_j;
    for (int k = 0; k < n_rows/b_rows; k++)
    {
        for (int w = 0; w < n_cols/b_cols; w++)
	{
            for (int i = 0; i < b_rows; i++)
	    {
                for (int j = 0; j < b_cols; j++)
		{
                    block_ind = k*n_cols/b_cols + w;
		    glob_i = k * b_rows + i;
		    glob_j = w * b_cols + j;
		    block_dense[block_ind * b_size + i * b_cols + j] = A_dense[glob_i * n_cols + glob_j];
		}
	    }
	}
    }    

    const BSRMatrix* B = new BSRMatrix(n_rows, n_cols, b_rows, b_cols, block_dense);
    copy(B);

    delete B;
}

void BSRMatrix::copy(const CSRMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;

    std::vector<double> A_dense = A->to_dense();
    double block_dense[n_rows*n_cols];

    int block_ind, glob_i, glob_j;
    for (int k = 0; k < n_rows/b_rows; k++)
    {
        for (int w = 0; w < n_cols/b_cols; w++)
	{
            for (int i = 0; i < b_rows; i++)
	    {
                for (int j = 0; j < b_cols; j++)
		{
                    block_ind = k*n_cols/b_cols + w;
		    glob_i = k * b_rows + i;
		    glob_j = w * b_cols + j;
		    block_dense[block_ind * b_size + i * b_cols + j] = A_dense[glob_i * n_cols + glob_j];
		}
	    }
	}
    }

    const BSRMatrix* B = new BSRMatrix(n_rows, n_cols, b_rows, b_cols, block_dense);
    copy(B);

    delete B;
}

void BSRMatrix::copy(const CSCMatrix* A)
{
    printf("Currently not implemented\n");
}

void BSRMatrix::copy(const BSRMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    b_rows = A->b_rows;
    b_cols = A->b_cols;
    b_size = A->b_size;
    n_blocks = A->n_blocks;
    nnz = A->nnz;

    idx1.resize(A->idx1.size());
    idx2.resize(A->idx2.size());
    vals.resize(A->vals.size());

    std::copy(A->idx1.begin(), A->idx1.end(), idx1.begin());
    std::copy(A->idx2.begin(), A->idx2.end(), idx2.begin());
    std::copy(A->vals.begin(), A->vals.end(), vals.begin());
}

/**************************************************************
*****   BSRMatrix to_dense
**************************************************************
***** Converts the BSRMatrix into a dense matrix
***** in the form of a flattened vector ordered row-wise
**************************************************************/
std::vector<double> BSRMatrix::to_dense()
{
    std::vector<double> dense(n_rows * n_cols);
    std::fill(dense.begin(), dense.end(), 0.0);

    int start, end;
    int upper_i, upper_j, data_offset;
    double val;
    int glob_i, glob_j, ind;
    for (int i=0; i<n_rows/b_rows; i++)
    {
	start = idx1[i];
	end = idx1[i+1];
        for (int j=start; j<end; j++)
	{
            upper_i = i * b_rows;
	    upper_j = idx2[j] * b_cols;
	    data_offset = j * b_size;
	    for (int block_i = 0; block_i < b_rows; block_i++)
	    {
                for (int block_j = 0; block_j < b_cols; block_j++)
		{
                    glob_i = upper_i + block_i;
		    glob_j = upper_j + block_j;
		    ind = block_i * b_cols + block_j + data_offset;
		    val = vals[ind];
                    dense[glob_i*n_cols + glob_j] = val;
		}
	    }
	}
    }

    return dense;
}

/**************************************************************
*****  BSRMatrix Add Value
**************************************************************
**************************************************************/
void BSRMatrix::add_value(int row, int col, double value)
{
    printf("Currently not implemented\n");
}

/**************************************************************
*****  BSRMatrix Add Block
***************************************************************
***** Inserts nonzero block of values into position (row, col)
***** of the block-matrix - NOT the global indices.
***** Values of the non-zero block being added should be row
***** ordered within the values array.
*****
***** Input:
*****   row:
*****     row index of non-zero block in matrix - NOT global row
*****   col:
*****     col index of non-zero block in matrix - NOT global col
*****   values:
*****     vector of block values - ordered row-wise within the block
**************************************************************/
void BSRMatrix::add_block(int row, int col, std::vector<double>& values)
{
    // Only add correct number of elements for block if values is longer than
    // block size
    if (values.size() > b_size) values.erase(values.begin()+b_size, values.end());

    // Add zeros to end of values vector if smaller than block size
    if (values.size() < b_size)
    {
        for(int k=values.size(); k<b_size; k++)
	{
            values.push_back(0.0);
	}
    }

    int start, end, j, data_offset;
    start = idx1[row];
    end = idx1[row+1];
    data_offset = idx1[row] * b_size;

    // Update cols vector and data offset
    if(idx2.size() < 1)
    {
        // First block added
        idx2.push_back(col);
	data_offset = 0;
    }
    else if(start == end || col > idx2[end-1])
    {
        idx2.insert(idx2.begin()+end, col);
	data_offset += b_size * (end-start);
    }
    else if(col < idx2[start]){
        idx2.insert(idx2.begin()+start, col);
    }
    else
    {
        while(j < end)
        {
            if(col < idx2[j])
            {
                idx2.insert(idx2.begin()+j, col);
		data_offset += b_size * (j-start);
            }
	    else j++;
        }
    }

    // Update rowptr
    for(int i=row+1; i<idx1.size(); i++){
        idx1[i]++;
    }

    // Update vals array
    vals.insert(vals.begin()+data_offset, values.begin(), values.end());

    // Update matrix variables
    nnz += b_size;
    n_blocks++;

}

/**************************************************************
*****   BSRMatrix Sort
**************************************************************
**************************************************************/
void BSRMatrix::sort()
{
    //printf("Currently not implemented\n");
    return;
}

void BSRMatrix::move_diag()
{
    printf("Currently not implemented\n");
}

void BSRMatrix::remove_duplicates()
{
    printf("Currently not implemented\n");
}

/**************************************************************
*****  CSCMatrix Add Value
**************************************************************
***** Inserts value into the position (row, col) of the matrix.
***** Values must be inserted in column-wise order, so if the col
***** is not equal to the col of the previously inserted value,
***** indptr is edited, and it is assumed that col is complete.
***** TODO -- this method needs further testing
*****
***** Parameters
***** -------------
***** row : int
*****    Row in which to insert value 
***** col : int
*****    Column in which to insert value
***** value : double
*****    Nonzero value to be inserted into the matrix
**************************************************************/
void CSCMatrix::add_value(int row, int col, double value)
{
    // Assumes idx1 is created separately
    idx2.push_back(row);
    vals.push_back(value);
    nnz++;
}

void CSCMatrix::add_block(int row, int col, std::vector<double>& values){
    printf("Not implemented.\n");
}

void CSCMatrix::copy(const COOMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = A->nnz;

    idx1.resize(n_rows + 1);
    std::fill(idx1.begin(), idx1.end(), 0);
    if (nnz)
    {
        idx2.resize(nnz);
        vals.resize(nnz);
    }

    // Calculate indptr
    for (int i = 0; i < n_cols + 1; i++)
    {
        idx1[i] = 0;
    }
    for (int i = 0; i < A->nnz; i++)
    {
        idx1[A->idx2[i]+1]++;
    }
    for (int i = 0; i < A->n_cols; i++)
    {
        idx1[i+1] += idx1[i];
    }

    // Add indices and data
    std::vector<int> ctr;
    if (n_cols)
    {
        ctr.resize(n_cols, 0);
    }
    for (int i = 0; i < A->nnz; i++)
    {
        int row = A->idx1[i];
        int col = A->idx2[i];
        double val = A->vals[i];
        int index = idx1[col] + ctr[col]++;
        idx2[index] = row;
        vals[index] = val;
    }        
}
void CSCMatrix::copy(const CSRMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = A->nnz;

    // Resize vectors to appropriate dimensions
    idx1.resize(A->n_cols + 1);
    if (A->nnz)
    {
        idx2.resize(A->nnz);
        vals.resize(A->nnz);
    }

    // Create indptr, summing number times col appears in CSR
    for (int i = 0; i <= A->n_cols; i++) 
    {
        idx1[i] = 0;
    }
    for (int i = 0; i < A->nnz; i++)
    {
        idx1[A->idx2[i] + 1]++;
    }
    for (int i = 0; i < A->n_cols; i++)
    {
        idx1[i+1] += idx1[i];
    }

    // Add values to indices and data
    if (A->n_cols)
    {
        std::vector<int> ctr(A->n_cols, 0);
        for (int i = 0; i < A->n_rows; i++)
        {
            int row_start = A->idx1[i];
            int row_end = A->idx1[i+1];
            for (int j = row_start; j < row_end; j++)
            {
                int col = A->idx2[j];
                int idx = idx1[col] + ctr[col]++;
                idx2[idx] = i;
                vals[idx] = A->vals[j];
            }
        }
    }
}

void CSCMatrix::copy(const CSCMatrix* A)
{
    n_rows = A->n_rows;
    n_cols = A->n_cols;
    nnz = A->nnz;

    idx1.resize(A->n_cols + 1);
    idx2.resize(A->nnz);
    vals.resize(A->nnz);

    idx1[0] = 0;
    for (int i = 0; i < A->n_cols; i++)
    {
        int col_start = A->idx1[i];
        int col_end = A->idx1[i+1];
        idx1[i+1] = col_end;
        for (int j = col_start; j < col_end; j++)
        {
            idx2[j] = A->idx2[j];
            vals[j] = A->vals[j];
        }
    }
}

void CSCMatrix::copy(const BSRMatrix* A)
{
    printf("Currently not implemented");
    return;
}

/**************************************************************
*****   CSCMatrix Sort
**************************************************************
***** Sorts the sparse matrix by rows within each column.  
***** Removes duplicates, summing their values 
***** together.
**************************************************************/
void CSCMatrix::sort()
{
    int start, end, col_size;
    int prev_k, k;

    std::vector<int> permutation;
    std::vector<bool> done;

    if (sorted || nnz == 0)
    {
        return;
    }

    // Sort the columns of each row (and data accordingly) and remove
    // duplicates (summing values together)
    for (int col = 0; col < n_cols; col++)
    {
        start  = idx1[col];
        end = idx1[col+1];
        col_size = end - start;
        if (col_size == 0)
        {
            continue;
        }

        // Create permutation vector p for row
        permutation.resize(col_size);
        std::iota(permutation.begin(), permutation.end(), 0);
        std::sort(permutation.begin(), permutation.end(),
                [&](int i, int j)
                { 
                    return idx2[i + start] < idx2[j + start];
                });

        // Permute columns and data according to p
        done.resize(col_size);
        for (int i = 0; i < col_size; i++)
        {
            done[i] = false;
        }
        for (int i = 0; i < col_size; i++)
        {
            if (done[i]) continue;

            done[i] = true;
            prev_k = i;
            k = permutation[i];
            while (i != k)
            {
                std::swap(idx2[prev_k + start], idx2[k + start]);
                std::swap(vals[prev_k + start], vals[k + start]);
                done[k] = true;
                prev_k = k;
                k = permutation[k];
            }
        }
    }

    sorted = true;
    diag_first = false;
}

void CSCMatrix::remove_duplicates()
{
    if (!sorted)
    {
        sort();
        diag_first = false;
    }
    
    int orig_start, orig_end, new_start;
    int col_size;
    int row, prev_row, ctr;
    double val;

    // Sort the columns of each row (and data accordingly) and remove
    // duplicates (summing values together)
    orig_start = idx1[0];
    for (int col = 0; col < n_cols; col++)
    {
        new_start = idx1[col];
        orig_end = idx1[col+1];
        col_size = orig_end - orig_start;
        if (col_size == 0) 
        {
            orig_start = orig_end;
            idx1[col+1] = idx1[col];
            continue;
        }

        // Remove Duplicates
        row = idx2[orig_start];
        val = vals[orig_start];
        idx2[new_start] = row;
        vals[new_start] = val;
        prev_row = row;
        ctr = 1;
        for (int j = orig_start + 1; j < orig_end; j++)
        {
            row = idx2[j];
            val = vals[j];
            if (row == prev_row)
            {
                vals[ctr - 1 + new_start] += val;
            }
            else
            {
                idx2[ctr + new_start] = row;
                vals[ctr + new_start] = val;
                ctr++;
                prev_row = row;
            }
        }
        orig_start = orig_end;
        idx1[col+1] = idx1[col] + ctr;
    }
}

void CSCMatrix::move_diag()
{
    if (diag_first || nnz == 0)
    {
       return;
    }

    int start, end, row;
    double tmp;

    // Move diagonal values to beginning of each column
    for (int i = 0; i < n_cols; i++)
    {
        start = idx1[i];
        end = idx1[i+1];
        for (int j = start; j < end; j++)
        {
            row = idx2[j];
            if (row == i)
            {
                tmp = vals[j];
                for (int k = j; k > start; k--)
                {
                    idx2[k] = idx2[k-1];
                    vals[k] = vals[k-1];
                }
                idx2[start] = i;
                vals[start] = tmp;
                break;
            }
        }
    }

    diag_first = true;
}

#endif
