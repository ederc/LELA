/*
 * multiline-index.h
 *
 *  Created on: 14 juin 2012
 *      Author: martani
 */

#ifndef MULTILINE_INDEX_H_
#define MULTILINE_INDEX_H_
#include "FG-types.h"
#include <iterator>

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

	bool _index_maps_constructed;
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
		_index_maps_constructed = false;

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

		_index_maps_constructed = true;
	}

	void constructSubMatrices(const SparseMatrix<uint16>& M,
								SparseMultilineMatrix<uint16>& A,
								SparseMultilineMatrix<uint16>& B,
								SparseMultilineMatrix<uint16>& C,
								SparseMultilineMatrix<uint16>& D,
								bool destruct_original)
	{
		//Rpiv <- the values in pivots_positions
		//Cpiv <- the indexes 0..pivots_size
		std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

		if(_index_maps_constructed == false)
		{
			report << "Indexes not constructed yet, call Indexer.processMatrix() first" << std::endl;
			throw "Index maps not constructed, call Indexer.processMatrix() first";
		}

		typedef SparseMatrix<uint16>::Row::const_iterator const_iterator;
		//typename SparseMatrix<uint16>::ConstRow *row;
		const_iterator *it = new const_iterator[A.nb_lines_per_bloc ()];
		const_iterator *it_end = new const_iterator[A.nb_lines_per_bloc ()];
		//typename SparseVector<uint16>::const_iterator it1, it2;

		uint32 curr_piv_AB = 0;
		//uint32 nb_elts_A=0, nb_elts_B=0;

		//uint32 nb_lines_per_row = 2;		//A.nb_lines ()
		uint32 row[A.nb_lines_per_bloc()];
		uint16 last_filled_row=0; // row1=-1, row2=-1;

		register bool index_added = false;
		uint16 nb_pending_zeros = 0;

		for (uint32 i = 0; i < M.coldim (); ++i) {
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)			//non pivot row
				continue;

			if(last_filled_row < A.nb_lines_per_bloc())
			{
				row[last_filled_row] = pivot_rows_idxs_by_entry[i];
				++last_filled_row;
				
				//report << "Row added " <<row[last_filled_row-1] << std::endl;
				
				if(last_filled_row < A.nb_lines_per_bloc())	//on nblines - 1 we handle the elts
					continue;
			}
			
			for(uint16 j=0; j<last_filled_row; ++j)
			{
				it[j] = M[row[j]].begin ();
				it_end[j] = M[row[j]].end ();
			}

			for(uint32 j=0; j < M.coldim (); ++j)	//for all the columns
			{
				index_added = false;
				nb_pending_zeros = 0;
				
				if(pivot_columns_map[j] != MINUS_ONE)
					for(uint16 k=0; k<A.nb_lines_per_bloc (); ++k)		//for all lines in a multiline
					{
						if(it[k] != it_end[k] && it[k]->first == j)
						{
							if(!index_added)
							{
								A[curr_piv_AB].IndexData.push_back (pivot_columns_map[it[k]->first]);
								
								for(uint16 l=0; l<nb_pending_zeros; ++l)
									A[curr_piv_AB].ValuesData.push_back (0);	//add 0s for all null elements in lines before this one
								
								index_added = true;
							}
							A[curr_piv_AB].ValuesData.push_back (it[k]->second);
							++(it[k]);
						}
						else
							if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
								nb_pending_zeros++;
							else
								A[curr_piv_AB].ValuesData.push_back (0);
					}
				else
					for(uint16 k=0; k<A.nb_lines_per_bloc (); ++k)		//for all lines in a multiline
					{
						if(it[k] != it_end[k] && it[k]->first == j)
						{
							if(!index_added)
							{
								B[curr_piv_AB].IndexData.push_back (non_pivot_columns_map[it[k]->first]);
								for(uint16 l=0; l<nb_pending_zeros; ++l)
									B[curr_piv_AB].ValuesData.push_back (0);	//add 0s for all null elements in lines before this one
								
								index_added = true;
							}
							B[curr_piv_AB].ValuesData.push_back (it[k]->second);
							++(it[k]);
						}
						else
							if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
								nb_pending_zeros++;
							else
								B[curr_piv_AB].ValuesData.push_back (0);
					}
			}



			curr_piv_AB++;
			last_filled_row = 0;	//make sure we are ready to start fill for the next iteration

			/*if(destruct_original)
				for(uint16 j=0; j<A.nb_lines_per_bloc (); ++j)
					M[row[j]].free ();*/

		}

		//Fill the left rows
		report << "LEFT ROWS in AB: " << A.nb_lines_per_bloc() - last_filled_row << std::endl;
		
		if(last_filled_row != 0)
		{
			for(uint16 j=0; j<last_filled_row; ++j)
			{
				it[j] = M[row[j]].begin ();
				it_end[j] = M[row[j]].end ();
			}

			if(last_filled_row != A.nb_lines_per_bloc ())		//ADD 0s
				for(uint16 j=last_filled_row; j<A.nb_lines_per_bloc(); ++j)
				{
					it[j] = M[row[0]].end ();					//position to end
					it_end[j] = M[row[0]].end ();
				}

			for(uint32 j=0; j < M.coldim (); ++j)	//for all the columns
			{
				index_added = false;
				nb_pending_zeros = 0;

				if(pivot_columns_map[j] != MINUS_ONE)
					for(uint16 k=0; k<A.nb_lines_per_bloc (); ++k)		//for all lines in a multiline
					{
						if(it[k] != it_end[k] && it[k]->first == j)
						{
							if(!index_added)
							{
								A[curr_piv_AB].IndexData.push_back (pivot_columns_map[it[k]->first]);

								for(uint16 l=0; l<nb_pending_zeros; ++l)
									A[curr_piv_AB].ValuesData.push_back (0);	//add 0s for all null elements in lines before this one

								index_added = true;
							}
							A[curr_piv_AB].ValuesData.push_back (it[k]->second);
							++(it[k]);
						}
						else
							if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
								nb_pending_zeros++;
							else
								A[curr_piv_AB].ValuesData.push_back (0);
					}
				else
					for(uint16 k=0; k<A.nb_lines_per_bloc (); ++k)		//for all lines in a multiline
					{
						if(it[k] != it_end[k] && it[k]->first == j)
						{
							if(!index_added)
							{
								B[curr_piv_AB].IndexData.push_back (non_pivot_columns_map[it[k]->first]);
								for(uint16 l=0; l<nb_pending_zeros; ++l)
									B[curr_piv_AB].ValuesData.push_back (0);	//add 0s for all null elements in lines before this one

								index_added = true;
							}
							B[curr_piv_AB].ValuesData.push_back (it[k]->second);
							++(it[k]);
						}
						else
							if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
								nb_pending_zeros++;
							else
								B[curr_piv_AB].ValuesData.push_back (0);
					}
			}



			curr_piv_AB++;
			last_filled_row = 0;	//make sure we are ready to start fill for the next iteration

			/*if(destruct_original)
				for(uint16 j=0; j<A.nb_lines_per_bloc (); ++j)
					M[row[j]].free ();*/
		}
		
		report << "DONE AB\n";

		//MATRIX C and D
		uint32 curr_piv_CD = 0;
		for (uint32 i = 0; i < M.rowdim (); ++i) {
			if(non_pivot_rows_idxs[i] == MINUS_ONE)
				continue;

			if(last_filled_row < A.nb_lines_per_bloc())
			{
				row[last_filled_row] = non_pivot_rows_idxs[i];
				++last_filled_row;

				//report << "Row added " <<row[last_filled_row-1] << std::endl;

				if(last_filled_row < A.nb_lines_per_bloc())	//on nblines - 1 we handle the elts
					continue;
			}

			for(uint16 j=0; j<last_filled_row; ++j)
			{
				it[j] = M[row[j]].begin ();
				it_end[j] = M[row[j]].end ();
			}

			for(uint32 j=0; j < M.coldim (); ++j)	//for all the columns
			{
				index_added = false;
				nb_pending_zeros = 0;

				if(pivot_columns_map[j] != MINUS_ONE)
					for(uint16 k=0; k<C.nb_lines_per_bloc (); ++k)		//for all lines in a multiline
					{
						if(it[k] != it_end[k] && it[k]->first == j)
						{
							if(!index_added)
							{
								C[curr_piv_CD].IndexData.push_back (pivot_columns_map[it[k]->first]);

								for(uint16 l=0; l<nb_pending_zeros; ++l)
									C[curr_piv_CD].ValuesData.push_back (0);	//add 0s for all null elements in lines before this one

								index_added = true;
							}
							C[curr_piv_CD].ValuesData.push_back (it[k]->second);
							++(it[k]);
						}
						else
							if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
								nb_pending_zeros++;
							else
								C[curr_piv_CD].ValuesData.push_back (0);
					}
				else
					for(uint16 k=0; k<A.nb_lines_per_bloc (); ++k)		//for all lines in a multiline
					{
						if(it[k] != it_end[k] && it[k]->first == j)
						{
							if(!index_added)
							{
								D[curr_piv_CD].IndexData.push_back (non_pivot_columns_map[it[k]->first]);
								for(uint16 l=0; l<nb_pending_zeros; ++l)
									D[curr_piv_CD].ValuesData.push_back (0);	//add 0s for all null elements in lines before this one

								index_added = true;
							}
							D[curr_piv_CD].ValuesData.push_back (it[k]->second);
							++(it[k]);
						}
						else
							if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
								nb_pending_zeros++;
							else
								D[curr_piv_CD].ValuesData.push_back (0);
					}
			}

			curr_piv_CD++;
			last_filled_row = 0;	//make sure we are ready to start fill for the next iteration
		}

		//Fill the left rows
		report << "LEFT ROWS in CD: " << A.nb_lines_per_bloc() - last_filled_row << std::endl;

		if(last_filled_row != 0)
		{
			for(uint16 j=0; j<last_filled_row; ++j)
			{
				it[j] = M[row[j]].begin ();
				it_end[j] = M[row[j]].end ();
			}

			if(last_filled_row != C.nb_lines_per_bloc ())		//ADD 0s
				for(uint16 j=last_filled_row; j<C.nb_lines_per_bloc(); ++j)
				{
					it[j] = M[row[0]].end ();					//position to end
					it_end[j] = M[row[0]].end ();
				}

			for(uint32 j=0; j < M.coldim (); ++j)	//for all the columns
			{
				index_added = false;
				nb_pending_zeros = 0;

				if(pivot_columns_map[j] != MINUS_ONE)
					for(uint16 k=0; k<C.nb_lines_per_bloc (); ++k)		//for all lines in a multiline
					{
						if(it[k] != it_end[k] && it[k]->first == j)
						{
							if(!index_added)
							{
								C[curr_piv_CD].IndexData.push_back (pivot_columns_map[it[k]->first]);

								for(uint16 l=0; l<nb_pending_zeros; ++l)
									C[curr_piv_CD].ValuesData.push_back (0);	//add 0s for all null elements in lines before this one

								index_added = true;
							}
							C[curr_piv_CD].ValuesData.push_back (it[k]->second);
							++(it[k]);
						}
						else
							if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
								nb_pending_zeros++;
							else
								C[curr_piv_CD].ValuesData.push_back (0);
					}
				else
					for(uint16 k=0; k<C.nb_lines_per_bloc (); ++k)		//for all lines in a multiline
					{
						if(it[k] != it_end[k] && it[k]->first == j)
						{
							if(!index_added)
							{
								D[curr_piv_CD].IndexData.push_back (non_pivot_columns_map[it[k]->first]);
								for(uint16 l=0; l<nb_pending_zeros; ++l)
									D[curr_piv_CD].ValuesData.push_back (0);	//add 0s for all null elements in lines before this one

								index_added = true;
							}
							D[curr_piv_CD].ValuesData.push_back (it[k]->second);
							++(it[k]);
						}
						else
							if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
								nb_pending_zeros++;
							else
								D[curr_piv_CD].ValuesData.push_back (0);
					}
			}
		}
	}


};

#endif /* MULTILINE_INDEX_H_ */
