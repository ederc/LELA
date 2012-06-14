/*
 * indexer.h
 *
 *  Created on: 12 juin 2012
 *      Author: martani
 */

#ifndef INDEXER_H_
#define INDEXER_H_

#include <assert.h>

using namespace LELA;
using namespace FG_Types;


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

public:


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
	void processMatrix(Matrix& M)
	{
		this->coldim = M.coldim ();
		this->rowdim = M.rowdim ();

		typename Matrix::ConstRowIterator i_M;
		ArrayType curr_row_idx = 0, entry;

		initArrays(this->rowdim, this->coldim);

		Npiv = 0;

		//Sweeps the rows to identify the row pivots and column pivots
		for (i_M = M.rowBegin (); i_M < M.rowEnd (); ++i_M)
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


	/*inline uint32 bloc_at_index(uint32 row_idx, uint32 col_idx)
	{
		return
	}*/
	//template <typename SparseMatrix>
	void constructSubMatrices(SparseMatrix<uint16>& M, SparseBlocMatrix<uint16, uint16>& A,
							  HybridBlocMatrix<uint16, uint16>& B, SparseBlocMatrix<uint16, uint16>& C,
							  HybridBlocMatrix<uint16, uint16>& D, bool destruct_original)
	{
		//Rpiv <- the values in pivots_positions
		//Cpiv <- the indexes 0..pivots_size

		std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

		assert(A.rowdim () == B.rowdim ());
		assert(A.bloc_rowdim () == B.bloc_rowdim ());

		typename SparseMatrix<uint16>::ConstRow *row;
		typename SparseMatrix<uint16>::Row::const_iterator it;
		typename SparseMatrix<uint16>::Row::const_reverse_iterator rit;

		uint32 curr_piv_AB = 0;
		uint32 nb_elts_A=0, nb_elts_B=0;

		uint32 bloc_row_idx = 0;
		uint32 bloc_idx_in_row=0;

		uint32 ii = 0;
		uint32 i = M.coldim ();
		do {
			--i;
			bloc_row_idx =  curr_piv_AB / A.bloc_rowdim ();

		//for (uint32 i = 0; i < M.coldim (); ++i) {
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			report << "PIV " << curr_piv_AB << std::endl;

			row = &(M[pivot_rows_idxs_by_entry[i]]);

			nb_elts_A = 0;
			nb_elts_B = 0;

			//1. count the number of the elements in row Rpiv[i] with column index in Cpiv
			for(it = row->begin (); it != row->end (); ++it) {
				if(pivot_columns_map[it->first] != MINUS_ONE)
					nb_elts_A++;
			}

			nb_elts_B = row->size () - nb_elts_A;

			//2. preallocate memory and add the elements
			//A[curr_piv_AB].reserve (nb_elts_A);
			//B[curr_piv_AB].reserve (nb_elts_B);

			if(curr_piv_AB % A.bloc_rowdim () == 0)
			{
				A[bloc_row_idx].resize(bloc_row_idx + 1);			//initialize row of blocs every A.bloc_rowdim () pivots

				for(uint32 j=0; j<bloc_row_idx + 1; ++j)
					A[bloc_row_idx][j].init (A.bloc_rowdim (), A.bloc_coldim ());
			}


			for(it = row->begin (); it != row->end (); ++it) {		//reverse sweep the row
				if(pivot_columns_map[it->first] != MINUS_ONE)
				{
					bloc_idx_in_row = (this->Npiv -1 - pivot_columns_map[it->first]) / A.bloc_coldim ();

					report << "bloc_idx_in_row " << bloc_idx_in_row << std::endl;
					/*report << "Bloc dim " << A[bloc_row_idx][bloc_idx_in_row].rowdim () << " x "
							<< A[bloc_row_idx][bloc_idx_in_row].coldim () << std::endl;*/

					A[bloc_row_idx][bloc_idx_in_row][curr_piv_AB % A.bloc_rowdim ()].push_back
							(typename SparseMatrix<uint16>::Row::value_type(pivot_columns_map[it->first], it->second));
				}
				else
				{
					//B[curr_piv_AB].push_back (typename Matrix::Row::value_type(non_pivot_columns_map[it->first], it->second));
				}
			}



			curr_piv_AB++;

			if(destruct_original)
				M[pivot_rows_idxs_by_entry[i]].free ();

		} while(i>0);
	}
};


#endif /* INDEXER_H_ */
