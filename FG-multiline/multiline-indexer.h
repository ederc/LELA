/*
 * multiline-index.h
 *
 *  Created on: 14 juin 2012
 *      Author: martani
 */

#ifndef MULTILINE_INDEX_H_
#define MULTILINE_INDEX_H_

using namespace LELA;

class MultiLineIndexer {
private:
	void initArrays(uint32 rowSize, uint32 colSize)
	{
		assert(colSize > 0);
		assert(rowSize > 0);

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

	uint32 coldim, rowdim;

	//static const uint32 MINUS_ONE = 0xFFFFFFFF;	//stable for uint32!! lookout for other types
	static const uint32 MINUS_ONE = 0 - 1 ;
	static const unsigned char MINUS_ONE_8 = 0xFF;

	//Array of size of the original matrix M.coldim() that maps a pivot (resp non pivot) column index to
	//its index in the sub matrix A (resp B)
	//If the value at any index i is MINUS_ONE, it means that the i'th column is not a pivot
	uint32 *pivot_columns_map, *non_pivot_columns_map;

	//the reverse map of pivot_columns_map and non_pivot_columns_map.
	//Maps columns from submatrix A and B to the original matrix M
	uint32 *pivot_columns_rev_map, *non_pivot_columns_rev_map;

	//the indexes of the non pivot rows
	uint32 *non_pivot_rows_idxs;

	//the indexes of pivot rows sorted by their entry (pivot column)
	uint32 *pivot_rows_idxs_by_entry;

public:

	//the number of pivots
	uint32 Npiv;

	MultiLineIndexer ()
	{
		pivot_columns_map = NULL;
		non_pivot_columns_map = NULL;
		pivot_columns_rev_map = NULL;
		non_pivot_columns_rev_map = NULL;
		non_pivot_rows_idxs = NULL;
		pivot_rows_idxs_by_entry = NULL;

	}

	~MultiLineIndexer ()
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
		uint32 curr_row_idx = 0, entry;

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
	}

	void constructSubMatrices(SparseMatrix<uint16>& M,  SparseMultilineMatrix<uint16>& A,
														SparseMultilineMatrix<uint16>& B,
														SparseMultilineMatrix<uint16>& C,
														SparseMultilineMatrix<uint16>& D,
														bool destruct_original)
	{
		//Rpiv <- the values in pivots_positions
		//Cpiv <- the indexes 0..pivots_size
		std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);


		typename SparseMatrix<uint16>::ConstRow *row;
		typename SparseMatrix<uint16>::Row::const_iterator it1, it2;
		//typename SparseVector<uint16>::const_iterator it1, it2;

		uint32 curr_piv_AB = 0;
		uint32 nb_elts_A=0, nb_elts_B=0;

		uint32 nb_lines_per_row = 2;		//A.nb_lines ()
		uint32 row1=-1, row2=-1;

		uint32 ii = 0;
		for (uint32 i = 0; i < M.coldim (); ++i) {
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)			//non pivot row
				continue;

			if(row1==-1)
			{
				row1 = pivot_rows_idxs_by_entry[i];
				continue;
			}

			if(row2 != -1)
				throw "ERROR";
			else
				row2 = pivot_rows_idxs_by_entry[i];

			it1 = M[row1].begin ();
			it2 = M[row2].begin ();

			while(it1 != M[row1].end () && it2 != M[row2].end ())
			{
				if(it1->first < it2->first)
				{
					if(pivot_columns_map[it1->first] != MINUS_ONE)
					{
						A[curr_piv_AB].IndexData.push_back (pivot_columns_map[it1->first]);
						A[curr_piv_AB].ValuesData.push_back (it1->second);
						A[curr_piv_AB].ValuesData.push_back (0);
					}
					else
					{
						B[curr_piv_AB].IndexData.push_back (non_pivot_columns_map[it1->first]);
						B[curr_piv_AB].ValuesData.push_back (it1->second);
						B[curr_piv_AB].ValuesData.push_back (0);
					}
					++it1;
				}
				else if(it1->first > it2->first)
				{
					if(pivot_columns_map[it2->first] != MINUS_ONE)
					{
						A[curr_piv_AB].IndexData.push_back (pivot_columns_map[it2->first]);
						A[curr_piv_AB].ValuesData.push_back (0);
						A[curr_piv_AB].ValuesData.push_back (it2->second);
					}
					else
					{
						B[curr_piv_AB].IndexData.push_back (non_pivot_columns_map[it2->first]);
						B[curr_piv_AB].ValuesData.push_back (0);
						B[curr_piv_AB].ValuesData.push_back (it2->second);
					}
					++it2;
				}
				else	// it1->first == it2->first
				{
					if(pivot_columns_map[it1->first] != MINUS_ONE)
					{
						A[curr_piv_AB].IndexData.push_back (pivot_columns_map[it1->first]);
						A[curr_piv_AB].ValuesData.push_back (it1->second);
						A[curr_piv_AB].ValuesData.push_back (it2->second);
					}
					else
					{
						B[curr_piv_AB].IndexData.push_back (non_pivot_columns_map[it1->first]);
						B[curr_piv_AB].ValuesData.push_back (it1->second);
						B[curr_piv_AB].ValuesData.push_back (it2->second);
					}
					++it1;
					++it2;
				}
			}

			while(it1 != M[row1].end ())
			{
				if(pivot_columns_map[it1->first] != MINUS_ONE)
				{
					A[curr_piv_AB].IndexData.push_back (pivot_columns_map[it1->first]);
					A[curr_piv_AB].ValuesData.push_back (it1->second);
					A[curr_piv_AB].ValuesData.push_back (0);
				}
				else
				{
					B[curr_piv_AB].IndexData.push_back (non_pivot_columns_map[it1->first]);
					B[curr_piv_AB].ValuesData.push_back (it1->second);
					B[curr_piv_AB].ValuesData.push_back (0);
				}
				++it1;
			}

			while(it2 != M[row2].end ())
			{
				if(pivot_columns_map[it2->first] != MINUS_ONE)
				{
					A[curr_piv_AB].IndexData.push_back (pivot_columns_map[it2->first]);
					A[curr_piv_AB].ValuesData.push_back (0);
					A[curr_piv_AB].ValuesData.push_back (it2->second);
				}
				else
				{
					B[curr_piv_AB].IndexData.push_back (non_pivot_columns_map[it2->first]);
					B[curr_piv_AB].ValuesData.push_back (0);
					B[curr_piv_AB].ValuesData.push_back (it2->second);
				}
				++it2;
			}

			curr_piv_AB++;
			row1 = row2 = -1;

			if(destruct_original)
			{
				M[row1].free ();
				M[row2].free ();
			}
		}

		if(row1 != -1 && row2 == -1){
			it1 = M[row1].begin ();

			while(it1 != M[row1].end ())
			{
				if(pivot_columns_map[it1->first] != MINUS_ONE)
				{
					A[curr_piv_AB].IndexData.push_back (pivot_columns_map[it1->first]);
					A[curr_piv_AB].ValuesData.push_back (it1->second);
					A[curr_piv_AB].ValuesData.push_back (0);
				}
				else
				{
					B[curr_piv_AB].IndexData.push_back (non_pivot_columns_map[it1->first]);
					B[curr_piv_AB].ValuesData.push_back (it1->second);
					B[curr_piv_AB].ValuesData.push_back (0);
				}
				++it1;
			}

			curr_piv_AB++;

			if(destruct_original)
			{
				M[row1].free ();
			}
		}
	}


};

#endif /* MULTILINE_INDEX_H_ */
