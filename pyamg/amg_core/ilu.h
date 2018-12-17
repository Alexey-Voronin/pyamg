//#include "core/types.hpp"
//#include "core/matrix.hpp"
//#include "core/vector.hpp"
#include "types.h"
#include "matrix.h"
#include "vector.h"
#include <vector>
#include <math.h>

#include</usr/local/Cellar/eigen/3.3.5/include/eigen3/Eigen/Core>
#include</usr/local/Cellar/eigen/3.3.5/include/eigen3/Eigen/SVD>


using namespace raptor;


Matrix* COOMatrix::ilu_k(int lof)
{
	printf("Function not implemented for COO type\n");
	return NULL;
}

/*Matrix* COOMatrix::ilu_levels()
{
	printf("Function not implemented for COO type\n");
	return NULL;
}

Matrix* COOMatrix::ilu_sparsity(Matrix* levls, int lof)
{
	printf("Function not implemented for COO type\n");
	return NULL;
}

Matrix* COOMatrix::ilu_symbolic(int lof)
{
	printf("Function not implemented for COO type\n");
	return NULL;
}

std::vector<double> COOMatrix::ilu_numeric(Matrix* levls)
{
	printf("Function not implemented for COO type\n");
	std::vector<double> emp;
	return emp;
}*/


Matrix* CSRMatrix::ilu_k(int lof)
{
    printf("@CSRMatrix::ilu_k()\n");
	printf("Begin ilu symbolic phase \n");
	CSRMatrix* sparsity = ilu_symbolic(lof);

	printf("\n");
	
	printf("Begin ilu numeric phase \n");
	std::vector<double> factors_data = this->ilu_numeric(sparsity);
	printf("After numeric ilu phase completed\n");
	
	CSRMatrix* factors = new CSRMatrix(sparsity->n_rows,sparsity->n_cols);
	factors->nnz = sparsity->nnz;
	int factors_nnz_dense = factors->n_rows*factors->n_cols;

	factors->idx1.resize(sparsity->n_rows+1);
	if(factors_nnz_dense){
		factors->idx2.resize(sparsity->nnz);
		factors->vals.resize(sparsity->nnz);
	}

	for(int i = 0; i < sparsity->n_rows+1;i++){
		factors->idx1[i] = sparsity->idx1[i];
	}

	for(int i =0; i < sparsity->nnz;i++){
		factors->idx2[i] = sparsity->idx2[i];
		factors->vals[i] = factors_data[i];
	}
	
	return factors;
}

CSRMatrix* CSRMatrix::ilu_symbolic(int lof)
{
    printf("@CSRMatrix::ilu_symbolic()\n");
	printf("Begin levls phase \n");
	CSRMatrix* levels = this->ilu_levels();
	printf("Begin sparsity phase \n");
	CSRMatrix* sparsity = this->ilu_sparsity(levels, lof);
	return sparsity;
}

CSRMatrix* CSRMatrix::ilu_sparsity(CSRMatrix* levls, int lof)
{
    printf("@CSRMatrix::ilu_sparsity()\n");
	printf("Begin ilu sparsity \n");
	int adj_lof = lof + 1;

	CSRMatrix* sparsity = new CSRMatrix(levls->n_rows,levls->n_cols);

	sparsity->n_rows = levls->n_rows;
	sparsity->n_cols = levls->n_cols;
	int sparsity_nnz = 0;
	sparsity->nnz = 0;
	int sparsity_nnz_dense = sparsity->n_rows*sparsity->n_cols;

	sparsity->idx1.resize(levls->n_rows+1);
	if(sparsity_nnz_dense){
		sparsity->idx2.reserve(sparsity_nnz_dense);
		sparsity->vals.reserve(sparsity_nnz_dense);
	}
	
	for(int row_i = 0; row_i<sparsity->n_rows;row_i++){
		sparsity->idx1[row_i] = sparsity->nnz;
		int start_ri = levls->idx1[row_i];
		int end_ri = levls->idx1[row_i+1];
		for(int j = start_ri; j<end_ri;j++){
			int col_j = levls->idx2[j];
			int levl_j = levls->vals[j];
			if(levl_j <= adj_lof){
				//sparsity->idx2[sparsity_nnz] = (col_j);
				sparsity->idx2.push_back(col_j);
				//sparsity->vals[sparsity_nnz] = (levl_j);
				sparsity->vals.push_back(levl_j);
				sparsity_nnz++;
				sparsity->nnz = sparsity_nnz;
			}
		}
	}
	sparsity->idx1[sparsity->n_rows] = sparsity->nnz;
/*	
	printf("sparsity idx1 = ");
    for (auto i: sparsity->idx1)
		printf("%d ",i);
	printf("\n");
  
	printf("sparsity idx2 = ");
    for (auto i: sparsity->idx2)
		printf("%d ",i);
	printf("\n");
  	
	printf("sparsity vals = ");
    for (auto i: sparsity->vals)
		printf("%.2f ",i);
	printf("\n");
*/	

	return sparsity; 
}

std::vector<double> CSRMatrix::ilu_numeric(CSRMatrix* levls){
    printf("@CSRMatrix::ilu_numerics()\n");

	int m = n_rows;
	int n = n_cols;

	std::vector<double> factors_data(levls->nnz);
	std::fill(factors_data.begin(), factors_data.end(), 0.0);
/*	
	printf("factors data = ");
	for(auto i: factors)
		printf("%.2f ",i);
	printf("\n");
*/
	//initialize entries in factors to entries in A
	int index_i = 0;
	for(int i = 0; i<levls->nnz;i++){
		if(levls->vals[i] == 1){
			factors_data[i] = vals[index_i];
			index_i++;
		}
	}
/*
	printf("initialize factors data = ");
	for(auto i: factors)
		printf("%.2f ",i);
	printf("\n");
*/
	//get vector of diagonal elements
	std::vector<double> diag_vec(n_rows);
	std::fill(diag_vec.begin(),diag_vec.end(),0.0);
	
	for(int row = 0; row < n_rows; row++){
		int start_i = levls->idx1[row];
		int end_i = levls->idx1[row+1];
		for(int j = start_i; j<end_i; j++){
			int col = levls->idx2[j];
			double val = factors_data[j];
			if(row == col)
				diag_vec[row] = val;
		}
	}
/*
	printf("diag vec = ");
	for(auto i: diag_vec)
		printf("%.2f ",i);
	printf("\n");
*/	
	//ikj Gaussian Elimination
	for(int row_i = 1; row_i<n_rows; row_i++){
		//printf("row i = %d\n", row_i);
		std::vector<int> current_row_levls(n_rows);
		std::fill(current_row_levls.begin(),current_row_levls.end(),0);

		std::vector<double> current_row_factors(n_rows);
		std::fill(current_row_factors.begin(),current_row_factors.end(),0.0);
		
		int start_ri = levls->idx1[row_i];
		int end_ri = levls->idx1[row_i+1];

		for(int jj=start_ri;jj<end_ri;jj++){
			int col = levls->idx2[jj];
			current_row_levls[col] = levls->vals[jj];
			current_row_factors[col] = factors_data[jj];
		}

		for(int k = 0; k<row_i;k++){
			if(diag_vec[k] == 0){
				printf("Matrix is singular!\n");
				//throw error
			}

			//compute multiplier
			if(current_row_levls[k] != 0){
				current_row_factors[k]=current_row_factors[k]/diag_vec[k];
				if(k==row_i)
					diag_vec[k] = current_row_factors[k];
			}

			//get row k
			int start_rk = levls->idx1[k];
			int end_rk = levls->idx1[k+1];
			std::vector<double> row_k_factors(n_rows);
			std::fill(row_k_factors.begin(),row_k_factors.end(),0.0);

			for(int jj=start_rk;jj<end_rk;jj++){
				int col = levls->idx2[jj];
				row_k_factors[col]=factors_data[jj];
			}

			for(int j=k+1;j<n_rows;j++){
				if(current_row_levls[j] != 0){
					current_row_factors[j] = current_row_factors[j]-current_row_factors[k]*row_k_factors[j];
					if(j==row_i)
						diag_vec[j] = current_row_factors[j];
				}
			}
		}

		for(int jj=start_ri;jj<end_ri;jj++){
			int col = levls->idx2[jj];
			factors_data[jj] = current_row_factors[col];
		}
	}
	
	/*
	printf("factors = ");
    for (auto i: factors)
		printf("%.2f ",i);
	printf("\n");
 
 	printf("ilu numeric phase done\n");
	printf("\n");
	*/
	return factors_data;
}

CSRMatrix* CSRMatrix::ilu_levels()
{
    printf("@CSRMatrix*  CSRMatrix::ilu_k()\n");
	printf("Begin ilu levels \n");
	CSRMatrix * levls = new CSRMatrix(n_rows,n_cols);

	//initialize vectors for final levels matrix
	levls->n_rows = n_rows;
	levls->n_cols = n_cols;
	int levls_nnz = 0;
	levls->nnz = 0;
	int levls_nnz_dense = n_rows*n_cols;

	
	levls->idx1.resize(n_rows+1);
	if(levls_nnz_dense){
		levls->idx2.reserve(levls_nnz_dense);
		levls->vals.reserve(levls_nnz_dense);
	}

	levls->idx1[0] = 0;

	//copy first row indices into final levls matrix cz that doesn't change
	int start_r = idx1[0];
	int end_r = idx1[1];
	for(int i = start_r; i<end_r;i++){
		//printf("levels nnz = %d\n",levls_nnz);
		levls->idx2.push_back(idx2[i]);
		levls->vals.push_back(1);
		levls_nnz++;
		levls->nnz = levls_nnz;
		//printf("First row: Added %d at 0,%d\n",1,idx2[i]);
	}
	levls->idx1[1] = levls_nnz;
	
	/*	
	printf("Levls idx1 = ");
    for (auto i: levls->idx1)
		printf("%d ",i);
	printf("\n");
  
	printf("Levls idx2 = ");
    for (auto i: levls->idx2)
		printf("%d ",i);
	printf("\n");
  	
	printf("Levls vals = ");
    for (auto i: levls->vals)
		printf("%e ",i);
	printf("\n");
	*/
	
	//begin ILU process
	for(int row_i = 1; row_i < n_rows;row_i++){
		std::vector<int> current_row_levls(n_rows);
		std::fill(current_row_levls.begin(), current_row_levls.end(), 100);
		//get row i

		int start_ri = idx1[row_i];
		int end_ri = idx1[row_i+1];

		//initialize temporary levls and idx2 vectors for row i
		for(int jj = start_ri; jj < end_ri; jj++){
			current_row_levls[idx2[jj]] = 1;
		}

		//printf("Initialize temp levls and idx2 for row i\n");
		
		/*printf("current row Levls = ");
    	for (auto i: current_row_levls)
			printf("%d ",i);
		printf("\n");
  	
		printf("\n");*/
		
		
		for(int k = 0; k <row_i; k++){
			//get row k
			int start_rk = levls->idx1[k];
			int end_rk = levls->idx1[k+1];
			 
			std::vector<int> row_k_levls(n_rows);
			std::fill(row_k_levls.begin(), row_k_levls.end(), 100);
	
			for(int t =start_rk;t<end_rk;t++){
				row_k_levls[levls->idx2[t]] =  levls->vals[t];
			}

			for(int j = k+1; j<n_rows; j++)
				current_row_levls[j] = min(current_row_levls[j], current_row_levls[k]+row_k_levls[j]);
		}

		for(int jj=0; jj<n_rows;jj++){
			if(current_row_levls[jj]<100){
				levls->idx2.push_back(jj);
				levls->vals.push_back(current_row_levls[jj]);
				levls_nnz++;
				levls->nnz = levls->nnz + 1;
			}
		}
		levls->idx1[row_i+1] = levls->nnz;
	}
		/*printf("Levls idx1 = ");
    	for (auto i: levls->idx1)
			printf("%d ",i);
		printf("\n");
  
		printf("Levls idx2 = ");
    	for (auto i: levls->idx2)
			printf("%d ",i);
		printf("\n");
  	
		printf("Levls vals = ");
    	for (auto i: levls->vals)
			printf("%e ",i);
		printf("\n");*/
	

	return levls;    
}


Matrix* CSCMatrix::ilu_k(int lof)
{
	printf("Function not implemented \n");
	return NULL;
}

///////////////////////////BSR////////////////////////// 
BSRMatrix* BSRMatrix::ilu_k(int lof)
{
    printf("@Matrix* BSRMatrix::ilu_k(int lof)\n");
	printf("Begin ilu symbolic phase \n");
	BSRMatrix* sparsity = ilu_symbolic(lof);

	printf("\n");
	
	printf("Begin ilu numeric phase \n");
	std::vector<double> factors_data = this->ilu_numeric(sparsity);
	printf("After numeric ilu phase completed\n");

	//for testing without numeric phase	
	//std::vector<double> factors_data(b_size* sparsity->n_blocks);
	//std::fill(factors_data.begin(), factors_data.end(), 1.0);


	BSRMatrix* factors = new BSRMatrix(sparsity->n_rows,sparsity->n_cols,sparsity->b_rows,sparsity->b_cols);
	factors->nnz = sparsity->nnz;
	
	int factors_nb_dense = factors->n_rows/sparsity->b_rows * factors->n_cols/sparsity->b_cols;

	factors->idx1.resize(sparsity->n_rows/sparsity->b_rows+1);
	if(factors_nb_dense){
		factors->idx2.resize(sparsity->n_blocks);
		factors->vals.resize(b_size * sparsity->n_blocks);
	}

	for(int i = 0; i < sparsity->n_rows/sparsity->b_rows+1;i++)
		factors->idx1[i] = sparsity->idx1[i];

	for(int i =0; i < sparsity->n_blocks;i++)
		factors->idx2[i] = sparsity->idx2[i];

	for(int i =0; i<sparsity->n_blocks*b_size;i++)
		factors->vals[i] = factors_data[i];
	
	return factors;
}

BSRMatrix* BSRMatrix::ilu_levels()
{
    printf("@BSRMatrix* BSRMatrix::ilu_levels()\n");
	printf("Begin ilu levels \n");
	BSRMatrix * levls = new BSRMatrix(n_rows,n_cols,b_rows,b_cols);

	//initialize vectors for final levels matrix
	levls->n_rows = n_rows;
	levls->n_cols = n_cols;
	levls->b_rows = b_rows;
	levls->b_cols = b_cols;
	levls->b_size = b_size;

	int levls_nnz = 0;
	levls->nnz = 0;
	int levls_n_blocks = 0;
	levls->n_blocks = 0;
	
	int levls_nb_dense = n_rows/b_rows * n_cols/b_cols;

	
	levls->idx1.resize(n_rows/b_rows+1);
	if(levls_nb_dense){
		levls->idx2.reserve(levls_nb_dense);
		levls->vals.reserve(b_size*levls_nb_dense);
	}

	levls->idx1[0] = 0;

	//copy first row indices into final levls matrix cz that doesn't change
	int start_r = idx1[0];
	int end_r = idx1[1];
	for(int i = start_r; i<end_r;i++){
		//levls->idx2[levls_n_blocks] = idx2[i];
		levls->idx2.push_back(idx2[i]);
		for(int j = 0; j <b_size;j++){
			//levls->vals[j] = 1;
			levls->vals.push_back(1);
			levls_nnz++;
			levls->nnz = levls_nnz;
		}
		levls_n_blocks++;
		levls->n_blocks = levls_n_blocks;
		//printf("First row: Added %lf at 0,%d\n",val,col);
	}
	levls->idx1[1] = levls_n_blocks;
	

	/*	
	printf("Levls idx1 = ");
    for (auto i: levls->idx1)
		printf("%d ",i);
	printf("\n");
  
	printf("Levls idx2 = ");
    for (auto i: levls->idx2)
		printf("%d ",i);
	printf("\n");
  	
	printf("Levls vals = ");
    for (auto i: levls->vals)
		printf("%e ",i);
	printf("\n");
	*/
	
	//begin ILU process
	for(int row_i = 1; row_i < n_rows/b_rows;row_i++){
		std::vector<int> current_row_levls(n_cols*b_rows);
		std::fill(current_row_levls.begin(), current_row_levls.end(), 100);
		//get row i

		int start_ri = idx1[row_i];
		int end_ri = idx1[row_i+1];

		//initialize temporary levls and idx2 vectors for row i
		for(int jj = start_ri; jj < end_ri; jj++){
			for(int ii=0; ii<b_size;ii++)
				current_row_levls[idx2[jj]*b_size+ii] = 1;
		}

		//printf("Initialize temp levls and idx2 for row i\n");
		
		/*printf("current row Levls = ");
    	for (auto i: current_row_levls)
			printf("%d ",i);
		printf("\n");
  	
		printf("\n");*/
		
		
		for(int k = 0; k <row_i; k++){
			//get row k
			int start_rk = levls->idx1[k];
			int end_rk = levls->idx1[k+1]; 

			std::vector<int> row_k_levls(n_cols*b_rows);
			std::fill(row_k_levls.begin(), row_k_levls.end(), 100);

			for(int t =start_rk;t<end_rk;t++){
				for(int ii=0; ii<b_size;ii++){
						row_k_levls[levls->idx2[t]*b_size+ii] = levls->vals[t*b_size+ii];
				}
			}

			for(int j = k+1; j<n_rows/b_rows; j++){
				for(int ii = 0; ii<b_size; ii++)
					current_row_levls[j*b_size+ii] = min(current_row_levls[j*b_size+ii], current_row_levls[k*b_size]+row_k_levls[j*b_size+ii]);
			}

		}

		for(int jj=0; jj<n_cols/b_cols;jj++){
			if(current_row_levls[jj*b_size]<100){
				levls->idx2[levls_n_blocks] = jj;
				for(int ii=0; ii<b_size;ii++){
					levls->vals[levls_n_blocks*b_size+ii] = current_row_levls[jj*b_size+ii];
					levls_nnz++;
					levls->nnz = levls->nnz + 1;
				}
				levls_n_blocks++;
				levls->n_blocks = levls_n_blocks;
			}
		}
		levls->idx1[row_i+1] = levls->n_blocks;
	}
		/*printf("Levls idx1 = ");
    	for (auto i: levls->idx1)
			printf("%d ",i);
		printf("\n");
  
		printf("Levls idx2 = ");
    	for (auto i: levls->idx2)
			printf("%d ",i);
		printf("\n");
  	
		printf("Levls vals = ");
    	for (auto i: levls->vals)
			printf("%e ",i);
		printf("\n");*/
	

	return levls;    
}


BSRMatrix* BSRMatrix::ilu_sparsity(BSRMatrix* levls, int lof)
{
    printf("@BSRMatrix* BSRMatrix::ilu_sparsity(BSRMatrix* levls, int lof)\n");
	//printf("Begin ilu sparsity \n");
	int adj_lof = lof + 1;

	BSRMatrix* sparsity = new BSRMatrix(levls->n_rows,levls->n_cols,levls->b_rows,levls->b_cols);

	sparsity->n_rows = levls->n_rows;
	sparsity->n_cols = levls->n_cols;
	sparsity->b_rows = levls->b_rows;
	sparsity->b_cols = levls->b_cols;
	sparsity->b_size = levls->b_size;

	int sparsity_nnz = 0;
	sparsity->nnz = 0;
	int sparsity_n_blocks = 0;
	sparsity->n_blocks = 0;

	int sparsity_nb_dense = sparsity->n_rows/sparsity->b_rows * sparsity->n_cols/sparsity->b_cols;

	sparsity->idx1.resize(levls->n_rows/levls->b_rows+1);
	if(sparsity_nb_dense){
		sparsity->idx2.reserve(sparsity_nb_dense);
		sparsity->vals.reserve(b_size*sparsity_nb_dense);
	}
	
	for(int row_i = 0; row_i<sparsity->n_rows/sparsity->b_rows;row_i++){
		sparsity->idx1[row_i] = sparsity->n_blocks;
		int start_ri = levls->idx1[row_i];
		int end_ri = levls->idx1[row_i+1];

		for(int j = start_ri; j<end_ri;j++){
			int col_j = levls->idx2[j];
			int levl_j = levls->vals[j*b_size];
			if(levl_j <= adj_lof){
				//sparsity->idx2[sparsity_nnz] = (col_j);
				sparsity->idx2.push_back(col_j);
				for(int ii = 0; ii <b_size; ii++){
					sparsity->vals.push_back(levl_j);
					sparsity_nnz++;
					sparsity->nnz = sparsity_nnz;
				}
				sparsity_n_blocks++;
				sparsity->n_blocks = sparsity_n_blocks;
			}
		}
	}
	sparsity->idx1[sparsity->n_rows/sparsity->b_rows] = sparsity_n_blocks;
/*	
	printf("sparsity idx1 = ");
    for (auto i: sparsity->idx1)
		printf("%d ",i);
	printf("\n");
  
	printf("sparsity idx2 = ");
    for (auto i: sparsity->idx2)
		printf("%d ",i);
	printf("\n");
  	
	printf("sparsity vals = ");
    for (auto i: sparsity->vals)
		printf("%.2f ",i);
	printf("\n");
*/	

	return sparsity; 
}

BSRMatrix* BSRMatrix::ilu_symbolic(int lof)
{
    printf("@BSRMatrix* BSRMatrix::ilu_symbolic(int lof)\n");
	printf("Begin levls phase \n");
	BSRMatrix* levels = this->ilu_levels();
	printf("Begin sparsity phase \n");
	BSRMatrix* sparsity = this->ilu_sparsity(levels, lof);
	return sparsity;
}

std::vector<double> BSRMatrix::ilu_numeric(BSRMatrix* levls)
{
    printf("@vector<double> BSRMatrix::ilu_numeric(BSRMatrix* levls)\n");
	int m = n_rows;
	int n = n_cols;
	
	//update data in levls matrix to have the matrix entries of A
    printf("goin into fill_factors\n");
	std::vector<double> factors_data = this->fill_factors(levls);
    printf("came back from fill_factors\n");

	//get the diagonals
	std::vector<double> diag_vec = this->get_diag(levls, factors_data);
	
	for(int row_i = 1; row_i < n/b_rows; row_i++){
		std::vector<int> current_row_levls(m*b_rows);
		std::fill(current_row_levls.begin(), current_row_levls.end(), 0);
		std::vector<double> current_row_factors(m*b_rows);
		std::fill(current_row_factors.begin(), current_row_factors.end(), 0.0);
		
		int start_ri = levls->idx1[row_i];
		int end_ri = levls->idx1[row_i+1];
		
		for(int jj = start_ri; jj<end_ri; jj++){
			int col = levls->idx2[jj];
			for(int ii = 0; ii<b_size; ii++){
				current_row_levls[col*b_size+ii] = levls->vals[jj*b_size+ii];
				current_row_factors[col*b_size+ii] = factors_data[jj*b_size+ii];	
			}
		}
		for(int k = 0; k < row_i; k++){
			if(current_row_levls[k*b_size] != 0){
				//compute inverse of diagonal block
				std::vector<double> inv_diag_b = this->inv_diag_block(diag_vec, k);

				//multiply inverse by block to get multiplier
				std::vector<double> current_row_fac_b(b_size);
				for(int ii = 0; ii< b_size; ii++)
					current_row_fac_b[ii] = current_row_factors[k*b_size+ii];

				std::vector<double> mult_b = this->mult_b(inv_diag_b,current_row_fac_b);

				//update in data array
				for(int ii=0; ii<b_size; ii++)
					current_row_factors[k*b_size+ii] = mult_b[ii];
			}//end if 

			//get row k
			int start_rk = levls->idx1[k];
			int end_rk = levls->idx1[k+1];
			
			std::vector<double>row_k_factors(m*b_rows);
			std::fill(row_k_factors.begin(), row_k_factors.end(), 0.0);
			
			for(int jj=start_rk; jj<end_rk; jj++){
				int col = levls->idx2[jj];
				for(int ii = 0; ii<b_size; ii++)
					row_k_factors[col*b_size+ii] = factors_data[jj*b_size+ii];
			}
			
			for(int j = k+1; j< n_rows/b_rows; j++){
				if(current_row_levls[j*b_size] != 0){
					std::vector<double> current_row_fac_k(b_size);
					std::vector<double> row_k_fac_b(b_size);
					std::fill(current_row_fac_k.begin(), current_row_fac_k.end(), 0.0);
					std::fill(row_k_fac_b.begin(), row_k_fac_b.end(), 0.0);
					
					for(int ii=0; ii<b_size; ii++){
						current_row_fac_k[ii] = current_row_factors[k*b_size+ii];
						row_k_fac_b[ii] = row_k_factors[j*b_size+ii];
					}//end for

					std::vector<double> temp = this->mult_b(current_row_fac_k, row_k_fac_b); 
					
					for(int ii=0; ii<b_rows;ii++){
						for(int jj=0; jj<b_cols; jj++)
							current_row_factors[j*b_size+ii*b_rows+jj] = current_row_factors[j*b_size+ii*b_rows+jj] - temp[ii*b_rows+jj];
					}

					//if it's a diagonal block, update in diag_vec
					if(j==row_i){
						for(int ii=0; ii<b_size; ii++)
							diag_vec[j*b_size+ii] = current_row_factors[j*b_size+ii];
					}//end if
				}//end if
			}//end j loop
		}//end k loop

		for(int jj = start_ri; jj<end_ri; jj++){
			int col = levls->idx2[jj];
			for(int ii=0; ii<b_size; ii++)
				factors_data[jj*b_size+ii] = current_row_factors[col*b_size+ii];
		}
	}//end i loop

	return factors_data;
}


std::vector<double> BSRMatrix::get_diag(BSRMatrix* levls, std::vector<double> data)
{
    printf("@vector<double> BSRMatrix::get_diag(BSRMatrix* levls, std::vector<double> data)\n");
	std::vector<double> diag_vec(n_rows/b_rows * b_size);
	std::fill(diag_vec.begin(), diag_vec.end(), 0.0);

	for(int row = 0; row<n_rows/b_rows;row++){
		int start_i = levls->idx1[row];
		int end_i = levls->idx1[row+1];
		for(int j = start_i; j < end_i; j++){
			int col = levls->idx2[j];
			if(row == col){
				for(int ii = 0; ii < b_size; ii++)
					diag_vec[row*b_size + ii] = data[j*b_size + ii];
			}
		}
	}
	return diag_vec;
}

std::vector<double> BSRMatrix::fill_factors(BSRMatrix* levls){
    printf("@vector<double> BSRMatrix::fill_factors(BSRMatrix* levls)\n");
    std::vector<double> factors_data(levls->nnz);
	std::fill(factors_data.begin(), factors_data.end(), 0.0);

	int index_i = 0;
	for(int i = 0; i < levls->nnz; i++){
		if(levls->vals[i] == 1){
			factors_data[i] = vals[index_i];
			index_i++;
		}
		else
			factors_data[i] = 0.0;
	}
    return factors_data;
}

 template<typename _Matrix_Type_>
 _Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon())
 {
     Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
     double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
     return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix(). asDiagonal() * svd.matrixU().adjoint();
 }

std::vector<double> BSRMatrix::inv_diag_block(std::vector<double> diag_vec, int k){
    printf("@vector<double> BSRMatrix::inv_diag_block(std::vector<double> diag_vec, int k)\n");
	//Compute inverse of diagonal block A[k,k]

    printf("len(diag_vec)=%d\tk=%d\tb_dim(%d,%d)\tb_size=%d\n", 
            diag_vec.size(), k, b_rows, b_cols,b_size);

	//extract diagonal block
	std::vector<double> diag_b(b_size);
	std::fill(diag_b.begin(), diag_b.end(), 0.0);
	double inv_det;

	for(int ii=0; ii< b_size; ii++)
		diag_b[ii] = diag_vec[k*b_size + ii];

	//Find inverse of diagonal block
	std::vector<double> inv_diag_b(b_size);
	std::fill(inv_diag_b.begin(), inv_diag_b.end(), 0.0);


    // Pseudo-inverse (Moorse-Penrose) using Eigen
     Eigen::MatrixXd A(b_rows,b_cols);
     printf("A:\n");
     for (int i = 0; i < b_rows; i++)
     {
         for (int j = 0; j < b_cols; j++)
         {
             A(i,j) = (diag_vec.data())[i*k+j];
             printf("%.3f\t", A(i,j) );
         }
         printf("\n");
     }
     Eigen::MatrixXd Ainv = pseudoInverse(A);
     std::vector<double> out (k*k);
     printf("Ainv:\n");
     for (int i = 0; i < b_rows; i++)
     {    
          for (int j = 0; j < b_cols; j++)
          {
              (inv_diag_b.data())[i*b_rows+j] = Ainv(i, j);
              printf("%.3f\t", Ainv(i,j) );
          }
          printf("\n");
     }
     Eigen::MatrixXd I = A*Ainv;
     printf("A*Ainv:\n");
      for (int i = 0; i < b_rows; i++)
      {
           for (int j = 0; j < b_cols; j++)
           {
               printf("%.3f\t", I(i,j) );
           }
           printf("\n");
      }
/*
//  A     = | a  b |
//          | c  d |
//  A_inv = | d  -b|
//          |-c   a| * (1/(a*d-b*c))

	if(b_rows == 2){
		inv_det = 1.0/(diag_b[0]*diag_b[3]-diag_b[1]*diag_b[2]);
		inv_diag_b[0] = inv_det * diag_b[3];
		inv_diag_b[1] = -inv_det * diag_b[1];
		inv_diag_b[2] = -inv_det * diag_b[2];
		inv_diag_b[3] = inv_det * diag_b[0];
	}
	else if(b_rows == 3){
		inv_det = 1.0/( diag_b[0]*(diag_b[4]*diag_b[8] - diag_b[5]*diag_b[7]) - diag_b[1]*(diag_b[3]*diag_b[8] - diag_b[5]*diag_b[6]) + diag_b[2]*(diag_b[3]*diag_b[7] - diag_b[4]*diag_b[6]) );  
		
		inv_diag_b[0] = inv_det*(diag_b[4]*diag_b[8] - diag_b[5]*diag_b[7]);
		inv_diag_b[1] = inv_det*(diag_b[2]*diag_b[7] - diag_b[1]*diag_b[8]);
		inv_diag_b[2] = inv_det*(diag_b[1]*diag_b[5] - diag_b[2]*diag_b[4]);

		inv_diag_b[3] = inv_det*(diag_b[5]*diag_b[6] - diag_b[3]*diag_b[8]);
		inv_diag_b[4] = inv_det*(diag_b[0]*diag_b[8] - diag_b[2]*diag_b[6]);
		inv_diag_b[5] = inv_det*(diag_b[2]*diag_b[3] - diag_b[0]*diag_b[5]);

		inv_diag_b[6] = inv_det*(diag_b[3]*diag_b[7] - diag_b[4]*diag_b[6]);
		inv_diag_b[7] = inv_det*(diag_b[1]*diag_b[6] - diag_b[0]*diag_b[7]);
		inv_diag_b[8] = inv_det*(diag_b[0]*diag_b[4] - diag_b[1]*diag_b[3]);
	}
*/
	return inv_diag_b;
}

std::vector<double> BSRMatrix::mult_b(std::vector<double> block_a, std::vector<double> block_b)
{
    printf("@vector<double> BSRMatrix::mult_b(std::vector<double> block_a, std::vector<double> block_b)\n");
	std::vector<double> mult(b_size);
	std::fill(mult.begin(), mult.end(), 0.0);

	for(int ii = 0 ; ii< b_rows; ii++){
		for(int jj = 0; jj<b_cols; jj++){
			for(int kk = 0; kk<b_cols; kk++)
				mult[ii*b_rows+jj] = mult[ii*b_rows+jj] + block_a[ii*b_rows+kk]*block_b[kk*b_rows+jj];	
		}
	}

	return mult;
}
