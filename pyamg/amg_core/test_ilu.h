//#include "ILU/types.hpp"
//#include "ILU/matrix.hpp"
//#include "ILU/vector.hpp"
#include "types.h"
#include "matrix.h"
#include "vector.h"
#include "ilu.h"
//#include "ILU/Num.h"

//using namespace raptor; 

struct ILU_data
{
    bool init = false;
    int *idx1; int idx1_size;
    int *idx2; int idx2_size;
    double *vals; int vals_size;
};

ILU_data ilu_tmp;

template<class I, class F>
void ilu_bsr_fact(  
//input
                const I n,          const I lof,
                const I rows_in_block,   const I cols_in_block,
                const I row_ptr_arr[],   const I row_ptr_arr_size,
                const I cols_arr[],      const I cols_arr_size,
                const F val_arrIn[],     const I val_arrIn_size, //)
//output
                I row_ptr_out[],   I row_ptr_out_size,
                I cols_out[],      I cols_out_size,
                F val_out[],     I val_out_size,
                I out_len[],   I out_len_size)
{
    printf("Input data array sizes\n (%d,%d,%d)\n", row_ptr_arr_size, cols_arr_size, val_arrIn_size);
    std::vector<I>  row_ptr(row_ptr_arr,    row_ptr_arr +   row_ptr_arr_size);
	std::vector<I>  cols(   cols_arr,       cols_arr    +   cols_arr_size);
	std::vector<F>  vals(   val_arrIn,      val_arrIn   +   val_arrIn_size);

    raptor::BSRMatrix* A_bsr = new raptor::BSRMatrix(n, n, rows_in_block, cols_in_block, 
                                     row_ptr, cols, 
                                     vals);

    raptor::BSRMatrix* ilu_bsr = A_bsr->ilu_k(lof);

//    factors->print();

    
    raptor::COOMatrix* factors =  new raptor::COOMatrix(); 
    factors->copy(ilu_bsr);
    
    
    printf("sizes provided: (%d,%d,%d)\nsizes needed:  (%d,%d,%d)\n", 
                row_ptr_out_size, cols_out_size, val_out_size,
                factors->idx1.size(), factors->idx2.size(), factors->vals.size());
    if (row_ptr_out_size < factors->idx1.size() || cols_out_size < factors->idx2.size() || val_out_size < factors->vals.size())
    {
        printf("\n\n\nrow_ptr_out_size is too small!!!!\n\n\n");
        printf("\n\n\ncol_ptr_out_size is too small!!!!\n\n\n");
        printf("\n\n\nrow_ptr_out is too small!!!!\n\n\n");
        return;
    }
    out_len[0] = factors->idx1.size();
    out_len[1] = factors->idx2.size();
    out_len[2] = factors->vals.size();
    for (unsigned long i = 0; i < factors->idx1.size(); i++)
        row_ptr_out[i] = (factors->idx1.data())[i];
    for (unsigned long i = 0; i < factors->idx2.size(); i++)
        cols_out[i] = (factors->idx2.data())[i];
    for (unsigned long i = 0; i < factors->vals.size(); i++)
        val_out[i] = (factors->vals.data())[i];
    
/*
    I *idx1_arr; // I idx1_arr_size;
    I *idx2_arr; // I idx2_arr_size;
    F *val_arr;  // I val_arr_size;
    idx1_arr = (I *) malloc(factors->idx1.size()*sizeof(I));
    idx2_arr = (I *) malloc(factors->idx2.size()*sizeof(I));
    val_arr = (F *) malloc(factors->vals.size()*sizeof(F));
    for (unsigned long i = 0; i < factors->idx1.size(); i++i)
        idx1_arr[i] = (factors->idx1.data())[i]; 
    for (unsigned long i = 0; i < factors->idx2.size(); i++)
        idx2_arr[i] = (factors->idx2.data())[i]; 
    for (unsigned long i = 0; i < factors->vals.size(); i++)
        val_arr[i] = (factors->vals.data())[i]; 
    
    ilu_tmp.idx1 = idx1_arr; ilu_tmp.idx1_size = factors->idx1.size();
    ilu_tmp.idx2 = idx2_arr; ilu_tmp.idx2_size = factors->idx2.size();
    ilu_tmp.vals = val_arr; ilu_tmp.vals_size = factors->vals.size();
    ilu_tmp.init = true;
    */
}            
/*
template<class I>
I ilu_get_idx1_size()
{
    return (ilu_tmp.init ? ilu_tmp.idx1_size : 0);
}
template<class I>
I ilu_get_idx2_size()
{
    return (ilu_tmp.init ? ilu_tmp.idx2_size : 0);
}
template<class I>
I ilu_get_vals_size()
{
    return (ilu_tmp.init ? ilu_tmp.vals_size : 0);
}
//template<class F>
template<class F>
void ilu_get_vals(F arr[])
{
    if (ilu_tmp.init)
    {
        for (int i =0; i < ilu_tmp.vals_size; i++)
        {
            arr[i] = ilu_tmp.vals[i];
        }
    }
}
*/
/*
int main(int argc, char** argv){
	std::vector<int> row_ptr = {0,2,3,5};
	std::vector<int> cols = {0,1,1,1,2};
	std::vector<double> vals = {1.0,0.0,2.0,1.0,6.0,7.0,8.0,2.0,1.0,4.0,5.0,1.0,4.0,3.0,0.0,0.0,7.0,2.0,0.0,0.0};

	int rows_in_block = 2;
	int cols_in_block = 2;
	int n = 6;

	BSRMatrix* A_bsr = new BSRMatrix(n, n, rows_in_block, cols_in_block, row_ptr, cols, vals);
	A_bsr->print();

	int lof = 0;

	Matrix* factors = A_bsr->ilu_k(lof);
	factors->print();
	printf("After ilu k funtion returns\n");	


    int *idx1_arr; int *idx2_arr;
    double *val_arr;
    int idx1_arr_len, idx2_arr_len, val_arr_len;
    ilu_bsr_fact<int, double>(n, lof,
                 rows_in_block, cols_in_block,
                 row_ptr.data(),   row_ptr.size(),
                 cols.data(),      cols.size(),
                 vals.data(),      vals.size(),
 // return values
                idx1_arr,      idx1_arr_len,
                idx2_arr,      idx2_arr_len,
                val_arr,    val_arr_len);

    printf("\n\n\n--------------------");
    printf("\nidx1[%d]",factors->idx1.size());
    for (int i = 0; i < factors->idx1.size(); i++)
        printf("%d\t", (factors->idx1.data())[i]);
    printf("\nidx2[%d]",factors->idx2.size());
    for (int i = 0; i < factors->idx2.size(); i++)
        printf("%d\t", (factors->idx2.data())[i]);
    printf("\nvals[%d]=", factors->vals.size());
    for (int i = 0; i < factors->vals.size(); i++)
        printf("%f\t", (factors->vals.data())[i]);
    printf("\n");

    printf("--------------------");
    printf("\nidx1[%d]", idx1_arr_len);
    for (int i = 0; i < idx1_arr_len; i++)
        printf("%d\t", idx1_arr[i]);
    printf("\nidx2[%d]",idx2_arr_len);
    for (int i = 0; i < idx2_arr_len; i++)
        printf("%d\t", idx2_arr[i]);
    printf("\nvals[%d]=", val_arr_len);
    for (int i = 0; i < val_arr_len; i++)
        printf("%f\t", val_arr[i]);
    printf("\n");
	return 0;
}
*/
