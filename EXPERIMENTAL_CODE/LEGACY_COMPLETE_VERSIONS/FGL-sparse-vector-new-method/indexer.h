/*
 * indexer.h
 * Copyright 2012 Martani Fayssal (UPMC University Paris 06 / INRIA)
 *
 *  Created on: 30 mai 2012
 *      Author: martani (UPMC University Paris 06)
 */
 

#ifndef __INDEXER_H
#define __INDEXER_H

#include <assert.h>

using namespace LELA;


template <typename ArrayType = uint32>
class Indexer
{
private:
	void initArrays(uint32 rowSize, uint32 colSize)
	{
		assert(colSize > 0);
		assert(rowSize > 0);

		pivot_columns_map = new ArrayType [colSize];
		memset(pivot_columns_map, MINUS_ONE_8, colSize * sizeof(ArrayType));

		non_pivot_columns_map = new ArrayType [colSize];
		memset(non_pivot_columns_map, MINUS_ONE_8, colSize * sizeof(ArrayType));

		pivot_columns_rev_map = new ArrayType [colSize];
		memset(pivot_columns_rev_map, MINUS_ONE_8, colSize * sizeof(ArrayType));

		non_pivot_columns_rev_map = new ArrayType [colSize];
		memset(non_pivot_columns_rev_map, MINUS_ONE_8, colSize * sizeof(ArrayType));

		pivot_rows_idxs_by_entry = new ArrayType [colSize];
		memset(pivot_rows_idxs_by_entry, MINUS_ONE_8, colSize * sizeof(ArrayType));

		non_pivot_rows_idxs = new ArrayType [rowSize];
		memset(non_pivot_rows_idxs, MINUS_ONE_8, rowSize * sizeof(ArrayType));
	}

public:
	uint32 coldim, rowdim;

	//static const ArrayType MINUS_ONE = 0xFFFFFFFF;	//stable for uint32!! lookout for other types
	static const ArrayType MINUS_ONE = 0 - 1 ;
	static const unsigned char MINUS_ONE_8 = 0xFF;

	//Array of size of the original matrix M.coldim() that maps a pivot (resp non pivot) column index to 
	//its index in the sub matrix A (resp B)
	//If the value at any index i is MINUS_ONE, it means that the i'th column is not a pivot
	ArrayType *pivot_columns_map, *non_pivot_columns_map;

	//the reverse map of pivot_columns_map and non_pivot_columns_map. 
	//Maps columns from submatrix A and B to the original matrix M
	ArrayType *pivot_columns_rev_map, *non_pivot_columns_rev_map;

	//the indexes of the non pivot rows
	ArrayType *non_pivot_rows_idxs;

	//the indexes of pivot rows sorted by their entry (pivot column)
	ArrayType *pivot_rows_idxs_by_entry;



	//the number of pivots
	uint32 Npiv;

	Indexer ()
	{
		pivot_columns_map = NULL;
		non_pivot_columns_map = NULL;
		pivot_columns_rev_map = NULL;
		non_pivot_columns_rev_map = NULL;
		non_pivot_rows_idxs = NULL;
		pivot_rows_idxs_by_entry = NULL;

	}

	~Indexer ()
	{
		if(pivot_columns_map != NULL)
			delete [] pivot_columns_map;

		if(non_pivot_columns_map != NULL)
			delete [] non_pivot_columns_map;

		if(pivot_columns_rev_map != NULL)
			delete [] pivot_columns_rev_map;

		if(non_pivot_columns_rev_map != NULL)
			delete [] non_pivot_columns_rev_map;

		if(non_pivot_rows_idxs != NULL)
			delete [] non_pivot_rows_idxs;

		if(pivot_rows_idxs_by_entry != NULL)
			delete [] pivot_rows_idxs_by_entry;
	}

	//constructs the maps of indexes from original matrix to submatrices
	template <class Matrix>
	void processMatrix(const Matrix& M)
	{
		this->coldim = M.coldim ();
		this->rowdim = M.rowdim ();

		typename Matrix::ConstRowIterator i_M;
		ArrayType curr_row_idx = 0, entry;

		initArrays(this->rowdim, this->coldim);

		Npiv = 0;
		
		//Sweeps the rows to identify the row pivots and column pivots
		for (i_M = M.rowBegin (); i_M != M.rowEnd (); ++i_M)
		{
			if(!i_M->empty ())
			{
				entry = i_M->front ().first;
				if(pivot_rows_idxs_by_entry[entry] == MINUS_ONE)
				{
					pivot_rows_idxs_by_entry[entry] = curr_row_idx;
					Npiv++;
				}
				else		//TODO New choose the least sparse row //ELGAB Sylvain
				{
					//check if the pivot is less dense. replace it with this one
					if(M[pivot_rows_idxs_by_entry[entry]].size () > i_M->size ())
					{
						non_pivot_rows_idxs[pivot_rows_idxs_by_entry[entry]] = pivot_rows_idxs_by_entry[entry];
						pivot_rows_idxs_by_entry[entry] = curr_row_idx;
						//Npiv++;
					}
					else
						non_pivot_rows_idxs[curr_row_idx] = curr_row_idx;
				}
			}
			else
				non_pivot_rows_idxs[curr_row_idx] = curr_row_idx;

			++curr_row_idx;
		}

		//construct the column pivots (non column pivots) map and its reverse
		ArrayType piv_col_idx = 0, non_piv_col_idx = 0;
		for(uint32 i = 0; i < this->coldim; ++i){
			if(pivot_rows_idxs_by_entry[i] != MINUS_ONE)
			{
				pivot_columns_map[i] = piv_col_idx;
				pivot_columns_rev_map[piv_col_idx] = i;
				piv_col_idx++;
			}
			else
			{
				non_pivot_columns_map[i] = non_piv_col_idx;
				non_pivot_columns_rev_map[non_piv_col_idx] = i;
				non_piv_col_idx++;
			}
		}
	}

	template <typename Matrix>
	void constructSubMatrices(Matrix& M, Matrix& A, Matrix& B, Matrix& C, Matrix& D, bool destruct_original)
	{
		//Rpiv <- the values in pivots_positions
		//Cpiv <- the indexes 0..pivots_size

		typename Matrix::ConstRow *row;
		typename Matrix::Row::const_iterator it;

		uint32 curr_piv_AB = 0;
		uint32 nb_elts_A=0, nb_elts_B=0;

		uint32 ii = 0;
		for (uint32 i = 0; i < M.coldim (); ++i) {
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)			//non pivot row
				continue;

			row = &(M[pivot_rows_idxs_by_entry[i]]);

			assert(pivot_columns_map[row->front ().first] != MINUS_ONE);

			nb_elts_A = 0;
			nb_elts_B = 0;

			//1. count the number of the elements in row Rpiv[i] with column index in Cpiv
			for(it = row->begin (); it != row->end (); ++it) {
				if(pivot_columns_map[it->first] != MINUS_ONE)
					nb_elts_A++;
			}

			nb_elts_B = row->size () - nb_elts_A;

			//2. preallocate memory and add the elements
			A[curr_piv_AB].reserve (nb_elts_A);
			B[curr_piv_AB].reserve (nb_elts_B);

			assert(pivot_columns_map[row->front ().first] == ii++);

			for(it = row->begin (); it != row->end (); ++it) {
				if(pivot_columns_map[it->first] != MINUS_ONE)
				{
					assert(A[curr_piv_AB].size () < nb_elts_A);		//make sure it doesn't exceed reserved space
					A[curr_piv_AB].push_back (typename Matrix::Row::value_type(pivot_columns_map[it->first], it->second));
				}
				else
				{
					assert(B[curr_piv_AB].size () < nb_elts_B);
					B[curr_piv_AB].push_back (typename Matrix::Row::value_type(non_pivot_columns_map[it->first], it->second));
				}
			}

			curr_piv_AB++;

			if(destruct_original)
				M[pivot_rows_idxs_by_entry[i]].free ();
		}

		uint32 curr_piv_CD = 0;
		uint32 nb_elts_C=0, nb_elts_D=0;
		for (uint32 i = 0; i < M.rowdim (); ++i) {		//non pivot rows
			if(non_pivot_rows_idxs[i] == MINUS_ONE)
				continue;

			assert(non_pivot_rows_idxs[i] == i);		//assert that non_pivot_rows_idxs is constructed correctly

			row = &(M[non_pivot_rows_idxs[i]]);
			nb_elts_C = 0;
			nb_elts_D = 0;

			assert(pivot_rows_idxs_by_entry[row->front ().first] != MINUS_ONE);		//assert that this row is actually non pivot
			assert(pivot_rows_idxs_by_entry[row->front ().first] != non_pivot_rows_idxs[i]);	//assert that this row is not the pivot in that column index

			for(it = row->begin (); it != row->end (); ++it) {
				if(pivot_columns_map[it->first] != MINUS_ONE)
					nb_elts_C++;
			}

			nb_elts_D = row->size () - nb_elts_C;

			C[curr_piv_CD].reserve (nb_elts_C);
			D[curr_piv_CD].reserve (nb_elts_D);

			for(it = row->begin (); it != row->end (); ++it) {
				if(pivot_columns_map[it->first] != MINUS_ONE)
				{
					assert(C[curr_piv_CD].size () < nb_elts_C);
					C[curr_piv_CD].push_back (typename Matrix::Row::value_type(pivot_columns_map[it->first], it->second));
				}
				else
				{
					assert(D[curr_piv_CD].size () < nb_elts_D);
					D[curr_piv_CD].push_back (typename Matrix::Row::value_type(non_pivot_columns_map[it->first], it->second));
				}
			}

			curr_piv_CD++;

			if(destruct_original)
				M[non_pivot_rows_idxs[i]].free ();
		}
	}

	template <typename Matrix>
	void constructSubMatrices(Matrix& B, Matrix& D, Matrix& B1, Matrix& B2, Matrix& D1, Matrix& D2, bool destruct_original)
	{
		assert(B.coldim () == D.coldim ());
		assert(B.coldim () == B1.coldim () + B2.coldim ());
		assert(B.rowdim () == B1.rowdim ());
		assert(B.rowdim () == B2.rowdim ());
		assert(D.coldim () == D1.coldim () + D2.coldim ());
		assert(D1.rowdim () == D2.rowdim ());

		typename Matrix::ConstRow *row;
		typename Matrix::Row::const_iterator it;

		uint32 curr_piv_D1_D2 = 0;
		uint32 nb_elts_D1=0, nb_elts_D2=0;

		// Splice D into two submatrices D1 and D2. D1 contains pivot columns, D2 contains non pivot colums.
		// Notice that in D, all the rows must be poivot or null rows
		uint32 ii = 0;
		for (uint32 i = 0; i < D.coldim (); ++i) {
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			row = &(D[pivot_rows_idxs_by_entry[i]]);

			assert(pivot_columns_map[row->front ().first] != MINUS_ONE);

			nb_elts_D1 = 0;
			nb_elts_D2 = 0;

			//1. count the number of the elements in row Rpiv[i] with column index in Cpiv
			for(it = row->begin (); it != row->end (); ++it) {
				if(pivot_columns_map[it->first] != MINUS_ONE)
					nb_elts_D1++;
			}

			nb_elts_D2 = row->size () - nb_elts_D1;

			//2. preallocate memory and add the elements
			D1[curr_piv_D1_D2].reserve (nb_elts_D1);
			D2[curr_piv_D1_D2].reserve (nb_elts_D2);

			assert(pivot_columns_map[row->front ().first] == ii++);

			for(it = row->begin (); it != row->end (); ++it) {
				if(pivot_columns_map[it->first] != MINUS_ONE)
				{
					assert(D1[curr_piv_D1_D2].size () < nb_elts_D1);		//make sure it doesn't exceed reserved space
					D1[curr_piv_D1_D2].push_back (typename Matrix::Row::value_type(pivot_columns_map[it->first], it->second));
				}
				else
				{
					assert(D2[curr_piv_D1_D2].size () < nb_elts_D2);
					D2[curr_piv_D1_D2].push_back (typename Matrix::Row::value_type(non_pivot_columns_map[it->first], it->second));
				}
			}

			++curr_piv_D1_D2;

			if(destruct_original)
				D[pivot_rows_idxs_by_entry[i]].free ();
		}

		// Splice B into two submatrices B1 and B2, B1 contains pivot columns, B2 contains non pivot colums.
		// Pivot columns are found in the matrix D

		uint32 curr_piv_B1_B2 = 0;
		uint32 nb_elts_B1=0, nb_elts_B2=0;
		typename Matrix::RowIterator i_B;

		for(i_B = B.rowBegin (); i_B != B.rowEnd (); ++i_B)
		{
			nb_elts_B1 = 0;
			nb_elts_B2 = 0;

			//1. count the number of the elements in row Rpiv[i] with column index in Cpiv
			for(it = i_B->begin (); it != i_B->end (); ++it) {
				if(pivot_columns_map[it->first] != MINUS_ONE)
					nb_elts_B1++;
			}

			nb_elts_B2 = i_B->size () - nb_elts_B1;

			//2. preallocate memory and add the elements
			B1[curr_piv_B1_B2].reserve (nb_elts_B1);
			B2[curr_piv_B1_B2].reserve (nb_elts_B2);

			for(it = i_B->begin (); it != i_B->end (); ++it) {
				if(pivot_columns_map[it->first] != MINUS_ONE)
				{
					assert(B1[curr_piv_B1_B2].size () < nb_elts_B1);		//make sure it doesn't exceed reserved space
					B1[curr_piv_B1_B2].push_back (typename Matrix::Row::value_type(pivot_columns_map[it->first], it->second));
				}
				else
				{
					assert(B2[curr_piv_B1_B2].size () < nb_elts_B2);
					B2[curr_piv_B1_B2].push_back (typename Matrix::Row::value_type(non_pivot_columns_map[it->first], it->second));
				}
			}

			++curr_piv_B1_B2;

			if(destruct_original)
				i_B->free ();
		}

	}

	template <typename Matrix>
	void reconstructMatrix(Matrix& M, const Matrix& B2, const Matrix& D2)
	{
		assert(M.coldim () == this->coldim);
		assert(B2.rowdim () == this->Npiv);
		assert(B2.coldim () == D2.coldim ());

		typename Matrix::ConstRow *row;
		typename Matrix::Row *row_M;
		typename Matrix::ConstRowIterator i_B, i_D;
		typename Matrix::Row::const_iterator it;
		uint32 nb_elts = 0, piv=0;

		for(uint32 i=0; i < this->coldim; i++)
		{
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			if(this->pivot_rows_idxs_by_entry[i] >= this->Npiv)
				row = &(D2[this->pivot_rows_idxs_by_entry[i] - this->Npiv]);
			else
				row = &(B2[this->pivot_rows_idxs_by_entry[i]]);

			assert(piv < Npiv + D2.rowdim ());

			row_M = &(M[piv]);

			nb_elts = row->size () + 1;
	
			row_M->free ();
			row_M->reserve (nb_elts);

			row_M->push_back (typename Matrix::Row::value_type(i, 1));

			for(it = row->begin (); it != row->end (); ++it) {
				assert(non_pivot_columns_rev_map[it->first] != MINUS_ONE);
				assert(non_pivot_columns_rev_map[it->first] > i);

				row_M->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[it->first], it->second));
			}

			++piv;
		}
	}

	//when RREF is not needed, takes submatrics B and D and reconstructs M
	template <typename Matrix>
	void reconstructMatrixFromBD(Matrix& M, const Matrix& B, const Matrix& D)
	{
		assert(M.coldim () == this->coldim);
		assert(B.rowdim () == this->Npiv);
		assert(B.coldim () == D.coldim ());

		typename Matrix::ConstRow *row;
		typename Matrix::Row *row_M;
		typename Matrix::ConstRowIterator i_B, i_D;
		typename Matrix::Row::const_iterator it;
		uint32 nb_elts = 0, piv=0;

		for(uint32 i=0; i < this->coldim; i++)
		{
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			if(this->pivot_rows_idxs_by_entry[i] >= this->Npiv)
			{
				row = &(D[this->pivot_rows_idxs_by_entry[i] - this->Npiv]);
				nb_elts = row->size ();
			}
			else
			{
				row = &(B[this->pivot_rows_idxs_by_entry[i]]);
				nb_elts = row->size () + 1;
			}

			row_M = &(M[piv]);
			row_M->free ();
			row_M->reserve (nb_elts);

			if(this->pivot_rows_idxs_by_entry[i] < this->Npiv)
				row_M->push_back (typename Matrix::Row::value_type(i, 1));

			for(it = row->begin (); it != row->end (); ++it) {
				assert(non_pivot_columns_rev_map[it->first] != MINUS_ONE);
				assert(non_pivot_columns_rev_map[it->first] >= piv);

				row_M->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[it->first], it->second));
			}

			++piv;
		}
	}

	template <typename Matrix>
	void reconstructMatrix(Matrix& M, const Matrix& A, const Matrix& B, const Matrix& D)
	{
		lela_check(A.rowdim() == B.rowdim());
		lela_check(A.rowdim() == A.coldim());
		lela_check(D.coldim() == B.coldim());

		typename Matrix::ConstRow *rowA, *rowB, *rowD;
		typename Matrix::Row *row_M;
		typename Matrix::Row::const_iterator it_A, it_B, it_D;

		uint32 nb_elts = 0, piv=0;
		for(uint32 i=0; i < this->coldim; i++)
		{
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			if(this->pivot_rows_idxs_by_entry[i] >= this->Npiv)				//D
			{
				rowD = &(D[this->pivot_rows_idxs_by_entry[i] - this->Npiv]);
				nb_elts = rowD->size ();

				row_M = &(M[piv]);
				row_M->free ();
				row_M->reserve (nb_elts);

				for(it_D = rowD->begin (); it_D != rowD->end (); ++it_D) {
					assert(non_pivot_columns_rev_map[it_D->first] != MINUS_ONE);
					assert(non_pivot_columns_rev_map[it_D->first] >= piv);

					row_M->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[it_D->first], it_D->second));
				}
			}
			else															//A & B
			{
				rowA = &(A[this->pivot_rows_idxs_by_entry[i]]);
				rowB = &(B[this->pivot_rows_idxs_by_entry[i]]);

				nb_elts = rowA->size () + rowB->size ();

				row_M = &(M[piv]);
				row_M->free ();
				row_M->reserve (nb_elts);

				it_A = rowA->begin ();
				it_B = rowB->begin ();
				while(it_A != rowA->end () && it_B != rowB->end ())
				{
					if(pivot_columns_rev_map[it_A->first] < non_pivot_columns_rev_map[it_B->first])
					{
						assert(pivot_columns_rev_map[it_A->first] >= piv);
						row_M->push_back (typename Matrix::Row::value_type(pivot_columns_rev_map[it_A->first], it_A->second));
						++it_A;
					}
					else if(pivot_columns_rev_map[it_A->first] > non_pivot_columns_rev_map[it_B->first])
					{
						assert(non_pivot_columns_rev_map[it_B->first] >= piv);
						row_M->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[it_B->first], it_B->second));
						++it_B;
					}
					else
						throw "ERROR IN INDEX";
				}

				while(it_A != rowA->end ())
				{
					row_M->push_back (typename Matrix::Row::value_type(pivot_columns_rev_map[it_A->first], it_A->second));
					++it_A;
				}

				while(it_B != rowB->end ())
				{
					row_M->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[it_B->first], it_B->second));
					++it_B;
				}
			}

			++piv;
		}
	}

	///param only_rows is used to remap the new pivots only, no columns are considered;
	///this is basically used when the echelon form is wanted, and not the RREF
	void combineInnerIndexer(Indexer& inner_idxr, bool only_rows = false)
	{
		uint32 i, idx_in_B, idx_in_B2, rev_idx_in_A;

		//remap pivot rows indexes
		uint32 next_piv = 0;
		for (i = 0; i < this->coldim; ++i) {
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			pivot_rows_idxs_by_entry[i] = next_piv++;
		}

		next_piv = 0;
		if(only_rows)
		{
			for (i = 0; i < inner_idxr.coldim; ++i) {
				if(inner_idxr.pivot_rows_idxs_by_entry[i] == MINUS_ONE)
					continue;

				idx_in_B = i;
				rev_idx_in_A = non_pivot_columns_rev_map[idx_in_B];

				pivot_rows_idxs_by_entry[rev_idx_in_A] = this->Npiv + inner_idxr.pivot_rows_idxs_by_entry[i];
			}

			return;
		}


		//remap the pivot columns of D (their corresponding rows indexes)
		for (i = 0; i < inner_idxr.Npiv; ++i) {
			idx_in_B = inner_idxr.pivot_columns_rev_map[i];
			rev_idx_in_A = non_pivot_columns_rev_map[idx_in_B];

			// the pivot column at index rev_idx_in_A points to a row in
			// D (Npiv rows of A + the relative index in D)
			assert(pivot_rows_idxs_by_entry[rev_idx_in_A] == MINUS_ONE);
			pivot_rows_idxs_by_entry[rev_idx_in_A] = i + this->Npiv;
		}

		//remap the non pivot columns from indexer inner_idxr to this indexer 
		//(from B2, D2 to the outer matrix M0)
		for (i = 0; i < this->coldim; ++i) {
			if(non_pivot_columns_map[i] == MINUS_ONE)
				continue;

			//points the index of the non pivot colum in the matrix B = B1|B2 to the index in the matrix B2
			idx_in_B = non_pivot_columns_map[i];

			//column is a pivot column in B, skip it
			if(inner_idxr.pivot_columns_map[idx_in_B] != MINUS_ONE)	//TODO
				continue;

			assert(inner_idxr.non_pivot_columns_map[idx_in_B] != MINUS_ONE);

			idx_in_B2 = inner_idxr.non_pivot_columns_map[idx_in_B];

			non_pivot_columns_map[i] = idx_in_B2;
			non_pivot_columns_rev_map[idx_in_B2] = i;

			//the reverse map	//these tests should pass to ensure integrity
			/*rev_idx_in_B = inner_idxr.non_pivot_columns_rev_map[idx_in_B2];
			rev_idx_in_A = non_pivot_columns_rev_map[rev_idx_in_B];

			assert(rev_idx_in_B == idx_in_B);
			assert(rev_idx_in_A != MINUS_ONE);
			assert(rev_idx_in_A == i);

			non_pivot_columns_rev_map[idx_in_B2] = rev_idx_in_A;
			//non_pivot_columns_rev_map[idx_in_B2] = i; 			//equivalent
			 */

		}
	}


};


#endif 	//__INDEXER_H
