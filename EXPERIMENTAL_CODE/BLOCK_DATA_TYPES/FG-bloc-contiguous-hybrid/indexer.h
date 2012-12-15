/*
 * indexer.h
 *
 *  Created on: 6 juil. 2012
 *      Author: martani
 */

#ifndef INDEXER_H_
#define INDEXER_H_

#include "matrix-utils.h"
#include "../only-D/indexer.h"

#include "lela/util/debug.h"
#include "lela/util/commentator.h"

using namespace LELA;

class Indexer_
{
private:
	void freeArrays()
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

	void initArrays(uint32 rowSize, uint32 colSize)
	{
		lela_check(colSize > 0);
		lela_check(rowSize > 0);

		if(colSize == 0)
			colSize++;

		if(rowSize == 0)
			rowSize++;

		freeArrays();

		pivot_columns_map = new uint32 [colSize];
		memset(pivot_columns_map, MINUS_ONE_8, colSize * sizeof(uint32));

		non_pivot_columns_map = new uint32 [colSize];
		memset(non_pivot_columns_map, MINUS_ONE_8, colSize * sizeof(uint32));

		pivot_columns_rev_map = new uint32 [colSize];
		memset(pivot_columns_rev_map, MINUS_ONE_8, colSize * sizeof(uint32));

		non_pivot_columns_rev_map = new uint32 [colSize];
		memset(non_pivot_columns_rev_map, MINUS_ONE_8, colSize * sizeof(uint32));

		pivot_rows_idxs_by_entry = new uint32 [colSize];
		memset(pivot_rows_idxs_by_entry, MINUS_ONE_8, colSize * sizeof(uint32));

		non_pivot_rows_idxs = new uint32 [rowSize];
		memset(non_pivot_rows_idxs, MINUS_ONE_8, rowSize * sizeof(uint32));
	}

	static const uint32 MINUS_ONE = (uint32)-1 ;
	static const uint8 MINUS_ONE_8 = (uint8)-1;

	/** Array of size of the original matrix M.coldim() that maps a pivot (resp non pivot) column
	 * index to its corresponding index in the sub matrix A (resp B).
	 * If the value at any index i is MINUS_ONE, it means that the i'th column is not a pivot (resp. non pivot)
	 */
	uint32 *pivot_columns_map, *non_pivot_columns_map;

	/** The reverse map of pivot_columns_map and non_pivot_columns_map.
	 * Maps columns from submatrix A and B to the original matrix M
	 */
	uint32 *pivot_columns_rev_map, *non_pivot_columns_rev_map;

	/**
	 * The indexes of the non pivot rows
	 */
	uint32 *non_pivot_rows_idxs;

	/** Array of size of the original matrix M.coldim() that maps all pivot columns to
	 * their row index.
	 * Non pivot column entries have the values MINUS_ONE
	 */
	uint32 *pivot_rows_idxs_by_entry;

	bool _index_maps_constructed;

public:

	uint32 coldim, rowdim;

	/**
	 * The number of pivots in the matrix
	 */
	uint32 Npiv;

	Indexer_ ()
	{
		pivot_columns_map = NULL;
		non_pivot_columns_map = NULL;
		pivot_columns_rev_map = NULL;
		non_pivot_columns_rev_map = NULL;
		non_pivot_rows_idxs = NULL;
		pivot_rows_idxs_by_entry = NULL;

		_index_maps_constructed = false;
	}

	~Indexer_ ()
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

	/**
	 * Constructs the maps of indexes of the original matrix.
	 * Identifies the list of pivot/non pivot rows, pivot/non pivot columns and their reverse maps
	 * (reverse map is the correspondance from indexes in sub matrices to the original matrix)
	 */
	template <class Element>
	void processMatrix(const SparseMatrix<Element>& M)
	{
		this->coldim = M.coldim ();
		this->rowdim = M.rowdim ();

		typename SparseMatrix<Element>::ConstRowIterator i_M;
		uint32 curr_row_idx = 0, entry;

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
				else	//choose the least sparse row //ELGAB Sylvain
				{
					//check if the pivot is less dense. replace it with this one
					if(M[pivot_rows_idxs_by_entry[entry]].size () > i_M->size ())
					{
						non_pivot_rows_idxs[pivot_rows_idxs_by_entry[entry]] = pivot_rows_idxs_by_entry[entry];
						pivot_rows_idxs_by_entry[entry] = curr_row_idx;
					}
					else
						non_pivot_rows_idxs[curr_row_idx] = curr_row_idx;
				}
			}
			else
				non_pivot_rows_idxs[curr_row_idx] = curr_row_idx;

			++curr_row_idx;
		}

		//construct the column pivots (non column pivots) map and its reverse map
		uint32 piv_col_idx = 0, non_piv_col_idx = 0;
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

		/*lela_check(Npiv == piv_col_idx);
		lela_check(this->coldim - Npiv == non_piv_col_idx);

		for(uint32 i = 0; i < this->coldim; ++i){
			lela_check((pivot_columns_map[i] < Npiv) || (pivot_columns_map[i] == MINUS_ONE));
		}

		for(uint32 i = 0; i < this->coldim; ++i){
			lela_check((non_pivot_columns_map[i] < (this->coldim - Npiv))  || (non_pivot_columns_map[i] == MINUS_ONE));
		}*/

		_index_maps_constructed = true;
	}



	/**
	 * Constructs the submatrices A, B, C, D from matrix M
	 * In the case where the index maps are not yet constructed, the functions constructs the maps
	 * implicitly; in that case, the matrices are re_initiliazed with the new corresponding dimentions
	 */
	template<typename Element, typename Index>
	void constructSubMatrices(SparseMatrix<Element>& M,
			SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, 
			SparseBlocMatrix<ContiguousBloc<Element, Index> >& B,
			SparseBlocMatrix<ContiguousBloc<Element, Index> >& C, 
			SparseBlocMatrix<ContiguousBloc<Element, Index> >& D,
			bool destruct_original_matrix)
	{
		//std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

		typedef SparseBlocMatrix<ContiguousBloc<Element, Index> > Matrix;
		
		if (!_index_maps_constructed)
		{
			commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "Warning: Maps are not yet constructed. Constructing maps and reinitializing submatrices"
					<< endl;

			//throw std::runtime_error ("Indexes not constructed yet! call Indexer.processMatrix() first");
			processMatrix(M);

			A = Matrix (Npiv, Npiv, Matrix::ArrangementDownTop_RightLeft,
					false, A.bloc_height(), A.bloc_width());

			B = Matrix (Npiv, M.coldim() - Npiv, Matrix::ArrangementDownTop_LeftRight,
					true, B.bloc_height(), B.bloc_width()); //true: add fill blocs since they might be used

			C = Matrix (M.rowdim() - Npiv, Npiv, Matrix::ArrangementDownTop_RightLeft,
					false, C.bloc_height(), C.bloc_width());

			D = Matrix (M.rowdim() - Npiv, M.coldim() - Npiv, Matrix::ArrangementDownTop_LeftRight,
					true, D.bloc_height(), D.bloc_width()); //true: add fill blocs since they might be used
		}

		check_equal_or_raise_exception(A.rowdim(), A.coldim());
		check_equal_or_raise_exception(A.rowdim(), B.rowdim());
		check_equal_or_raise_exception(C.rowdim(), D.rowdim());
		check_equal_or_raise_exception(M.rowdim(), A.rowdim() + C.rowdim());

		check_equal_or_raise_exception(A.coldim(), C.coldim());
		check_equal_or_raise_exception(B.coldim(), D.coldim());
		check_equal_or_raise_exception(M.coldim(), A.coldim() + B.coldim());

		check_equal_or_raise_exception(M.rowdim(), this->rowdim);
		check_equal_or_raise_exception(M.coldim(), this->coldim);

		SparseMatrix<Element> s_A (A.rowdim(), A.coldim()),
							s_B (B.rowdim(), B.coldim()),
							s_C (C.rowdim(), C.coldim()),
							s_D (D.rowdim(), D.coldim());

		//TODO: change that immediately
		Indexer<uint32> one_line_indexer;
		one_line_indexer.processMatrix(M);
		one_line_indexer.constructSubMatrices(M, s_A, s_B, s_C, s_D, destruct_original_matrix);

		MatrixUtils::copy(s_A, A);
		MatrixUtils::copy(s_B, B);

		MatrixUtils::copy(s_C, C);
		MatrixUtils::copy(s_D, D);

	}

};


#endif /* INDEXER_H_ */
