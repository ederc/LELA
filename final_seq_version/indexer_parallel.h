/*
 * indexer.h
 *
 *  Created on: 6 juil. 2012
 *      Author: martani
 */

#ifndef INDEXER_PARALLEL_H_
#define INDEXER_PARALLEL_H_

//#include "matrix-utils.h"
//#include "../only-D/indexer.h"

#include "consts-macros.h"

#include "lela/util/debug.h"
#include "lela/util/commentator.h"

using namespace LELA;

template <typename Element, typename Index>
class ParallelIndexer
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
	int NB_THREADS;

public:

	uint32 coldim, rowdim;

	/**
	 * The number of pivots in the matrix
	 */
	uint32 Npiv;

	ParallelIndexer ()
	{
		pivot_columns_map = NULL;
		non_pivot_columns_map = NULL;
		pivot_columns_rev_map = NULL;
		non_pivot_columns_rev_map = NULL;
		non_pivot_rows_idxs = NULL;
		pivot_rows_idxs_by_entry = NULL;

		_index_maps_constructed = false;

		NB_THREADS = NUM_THREADS_OMP;
	}

	ParallelIndexer (int nb_threads) : NB_THREADS(nb_threads)
	{
		pivot_columns_map = NULL;
		non_pivot_columns_map = NULL;
		pivot_columns_rev_map = NULL;
		non_pivot_columns_rev_map = NULL;
		non_pivot_rows_idxs = NULL;
		pivot_rows_idxs_by_entry = NULL;

		_index_maps_constructed = false;
	}

	~ParallelIndexer ()
	{
		this->freeMemory();
	}

	void setNumThreads (int num_threads)
	{
		NB_THREADS = num_threads;
	}

	void freeMemory()
	{
		if(pivot_columns_map != NULL)
		{
			delete [] pivot_columns_map;
			pivot_columns_map = NULL;
		}

		if(non_pivot_columns_map != NULL)
		{
			delete [] non_pivot_columns_map;
			non_pivot_columns_map = NULL;
		}

		if(pivot_columns_rev_map != NULL)
		{
			delete [] pivot_columns_rev_map;
			pivot_columns_rev_map = NULL;
		}

		if(non_pivot_columns_rev_map != NULL)
		{
			delete [] non_pivot_columns_rev_map;
			non_pivot_columns_rev_map = NULL;
		}

		if(non_pivot_rows_idxs != NULL)
		{
			delete [] non_pivot_rows_idxs;
			non_pivot_rows_idxs = NULL;
		}

		if(pivot_rows_idxs_by_entry != NULL)
		{
			delete [] pivot_rows_idxs_by_entry;
			pivot_rows_idxs_by_entry = NULL;
		}
	}

	/**
	 * Constructs the maps of indexes of the original matrix.
	 * Identifies the list of pivot/non pivot rows, pivot/non pivot columns and their reverse maps
	 * (reverse map is the correspondance from indexes in sub matrices to the original matrix)
	 */
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

	void processMatrix(SparseMultilineMatrix<Element>& M)
	{
		this->coldim = M.coldim ();
		this->rowdim = M.rowdim ();

		typename SparseMultilineMatrix<Element>::ConstRowIterator i_M;
		uint32 curr_row_idx = 0;
		long h1, h2;
		uint32 h1_idx, h2_idx;
		Element h1_val, h2_val;

		//In case there is a null line at the end of a multiline
		initArrays(this->rowdim + this->rowdim % M.nb_lines_per_bloc (), this->coldim);

		Npiv = 0;

		//Sweeps the rows to identify the row pivots and column pivots
		for (i_M = M.rowBegin (); i_M != M.rowEnd (); ++i_M)
		{
			if(!i_M->empty ())
			{
				h1 = Level1Ops::headMultiLineVectorHybrid(*i_M, 0, h1_val, h1_idx, M.coldim());
				h2 = Level1Ops::headMultiLineVectorHybrid(*i_M, 1, h2_val, h2_idx, M.coldim());

				//LINE 1
				if(h1 == -1)	//empty line
				{
					non_pivot_rows_idxs[curr_row_idx] = curr_row_idx;
				}
				else
				{
					if(pivot_rows_idxs_by_entry[h1] == MINUS_ONE)		//non empty
					{
						pivot_rows_idxs_by_entry[h1] = curr_row_idx;
						Npiv++;
					}
					else		//TODO New choose the least sparse row //ELGAB Sylvain
					{
						//check if the pivot is less dense. replace it with this one
						/*if(M[pivot_rows_idxs_by_entry[h1]].size () > i_M->size ())
						{
							non_pivot_rows_idxs[pivot_rows_idxs_by_entry[h1]] = pivot_rows_idxs_by_entry[h1];
							pivot_rows_idxs_by_entry[h1] = curr_row_idx;
							//Npiv++;
						}
						else*/
						non_pivot_rows_idxs[curr_row_idx] = curr_row_idx;
					}
				}

				++curr_row_idx;
				//LINE 2
				if(h2 == -1)	//empty line
				{
					non_pivot_rows_idxs[curr_row_idx] = curr_row_idx;
				}
				else
				{
					if(pivot_rows_idxs_by_entry[h2] == MINUS_ONE)
					{
						pivot_rows_idxs_by_entry[h2] = curr_row_idx;
						Npiv++;
					}
					else		//TODO New choose the least sparse row //ELGAB Sylvain
					{
						//check if the pivot is less dense. replace it with this one
						/*if(M[pivot_rows_idxs_by_entry[h2]].size () > i_M->size ())
						{
							non_pivot_rows_idxs[pivot_rows_idxs_by_entry[h2]] = pivot_rows_idxs_by_entry[h2];
							pivot_rows_idxs_by_entry[h2] = curr_row_idx;
							//Npiv++;
						}
						else*/
						non_pivot_rows_idxs[curr_row_idx] = curr_row_idx;
					}
				}

				++curr_row_idx;
			}
			else		//multiline is empty
			{
				non_pivot_rows_idxs[curr_row_idx] = curr_row_idx;
				non_pivot_rows_idxs[curr_row_idx+1] = curr_row_idx+1;
				curr_row_idx += 2;
			}
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

//		assert(Npiv == piv_col_idx);
//		assert(this->coldim - Npiv == non_piv_col_idx);
//
//		for(uint32 i = 0; i < this->coldim; ++i){
//			assert((pivot_columns_map[i] < Npiv) || (pivot_columns_map[i] == MINUS_ONE));
//		}
//
//		for(uint32 i = 0; i < this->coldim; ++i){
//			assert((non_pivot_columns_map[i] < (this->coldim - Npiv))  || (non_pivot_columns_map[i] == MINUS_ONE));
//		}

		_index_maps_constructed = true;
	}



	void write_row_blocs_to_Left_Right_matrix(const SparseMatrix<Element>& M,
												SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
												SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
												uint32 *rows_idxs,
												uint32 nb_rows,
												uint32 row_bloc_idx)
	{
		typename SparseVector<Element>::const_iterator it1, it2;

		A[row_bloc_idx].clear ();
		B[row_bloc_idx].clear ();

		std::vector<Index> idx_stack_elt_idx_in_line;
		std::vector<uint32> idx_stack_bloc;
		std::vector<Element> v1_stack, v2_stack;

		uint32 bloc_idx_in_row, elt_idx_in_line;

		uint32 line_idx_in_bloc;
		uint32 nb_blocs_per_dim;

		B.FirstBlocsColumIndexes[row_bloc_idx] = 0;
		A.FirstBlocsColumIndexes[row_bloc_idx] = 0;

		nb_blocs_per_dim = (uint32) std::ceil((float)A.coldim () / A.bloc_width ());
		A[row_bloc_idx].resize (nb_blocs_per_dim);

		for(uint32 k=0; k<nb_blocs_per_dim; ++k)
			A[row_bloc_idx][k].init(A.bloc_height (), A.bloc_width ());

		nb_blocs_per_dim = (uint32) std::ceil((float)B.coldim () / B.bloc_width ());
		B[row_bloc_idx].resize (nb_blocs_per_dim);

		for(uint32 k=0; k<nb_blocs_per_dim; ++k)
			B[row_bloc_idx][k].init(B.bloc_height (), B.bloc_width ());

		//2. write elements to blocs
		uint32 i=0;
		for(; i< ROUND_DOWN(nb_rows, 2); i+=2)	//A.bloc_height () must be even
		{
			it1 = M[rows_idxs[i]].begin ();
			it2 = M[rows_idxs[i+1]].begin ();

			idx_stack_elt_idx_in_line.clear();
			idx_stack_bloc.clear();
			v1_stack.clear();
			v2_stack.clear();

			line_idx_in_bloc = i/2;

			while (it1 != M[rows_idxs[i]].end () && it2 != M[rows_idxs[i+1]].end ())
			{
				if(it1->first < it2->first)
				{
					if(pivot_columns_map[it1->first] != MINUS_ONE)		//pivot column
					{
						bloc_idx_in_row = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) / A.bloc_width();
						elt_idx_in_line = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) % A.bloc_width();

						idx_stack_elt_idx_in_line.push_back (elt_idx_in_line);
						idx_stack_bloc.push_back(bloc_idx_in_row);
						v1_stack.push_back (it1->second);
						v2_stack.push_back (0);
					}
					else
					{
						bloc_idx_in_row = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
						elt_idx_in_line = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it1->second);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
					}

					++it1;
				}
				else if(it1->first > it2->first)
				{
					if(pivot_columns_map[it2->first] != MINUS_ONE)		//pivot column
					{
						bloc_idx_in_row = (A.coldim () - 1 - pivot_columns_map[it2->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) / A.bloc_width();
						elt_idx_in_line = (A.coldim () - 1 - pivot_columns_map[it2->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) % A.bloc_width();

						idx_stack_elt_idx_in_line.push_back (elt_idx_in_line);
						idx_stack_bloc.push_back(bloc_idx_in_row);
						v1_stack.push_back (0);
						v2_stack.push_back (it2->second);
					}
					else
					{
						bloc_idx_in_row = (non_pivot_columns_map[it2->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
						elt_idx_in_line = (non_pivot_columns_map[it2->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it2->second);
					}

					++it2;
				}
				else	// it1->first == it2->first
				{
					if(pivot_columns_map[it1->first] != MINUS_ONE)		//pivot column
					{
						bloc_idx_in_row = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) / A.bloc_width();
						elt_idx_in_line = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) % A.bloc_width();

						idx_stack_elt_idx_in_line.push_back (elt_idx_in_line);
						idx_stack_bloc.push_back(bloc_idx_in_row);
						v1_stack.push_back (it1->second);
						v2_stack.push_back (it2->second);
					}
					else
					{
						bloc_idx_in_row = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
						elt_idx_in_line = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it1->second);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it2->second);
					}

					++it1;
					++it2;
				}
			}

			while (it1 != M[rows_idxs[i]].end())
			{
				if(pivot_columns_map[it1->first] != MINUS_ONE)		//pivot column
				{
					bloc_idx_in_row = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) / A.bloc_width();
					elt_idx_in_line = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) % A.bloc_width();

					idx_stack_elt_idx_in_line.push_back (elt_idx_in_line);
					idx_stack_bloc.push_back(bloc_idx_in_row);
					v1_stack.push_back (it1->second);
					v2_stack.push_back (0);
				}
				else
				{
					bloc_idx_in_row = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
					elt_idx_in_line = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it1->second);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
				}

				++it1;
			}

			while (it2 != M[rows_idxs[i+1]].end())
			{
				if(pivot_columns_map[it2->first] != MINUS_ONE)		//pivot column
				{
					bloc_idx_in_row = (A.coldim () - 1 - pivot_columns_map[it2->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) / A.bloc_width();
					elt_idx_in_line = (A.coldim () - 1 - pivot_columns_map[it2->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) % A.bloc_width();

					idx_stack_elt_idx_in_line.push_back (elt_idx_in_line);
					idx_stack_bloc.push_back(bloc_idx_in_row);
					v1_stack.push_back (0);
					v2_stack.push_back (it2->second);
				}
				else
				{
					bloc_idx_in_row = (non_pivot_columns_map[it2->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
					elt_idx_in_line = (non_pivot_columns_map[it2->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it2->second);
				}

				++it2;
			}

			//write into A (Right To Left)
			for (int k = idx_stack_elt_idx_in_line.size() - 1; k >= 0; --k)
			{
				A[row_bloc_idx][idx_stack_bloc[k]][line_idx_in_bloc].IndexData.push_back (idx_stack_elt_idx_in_line[k]);
				A[row_bloc_idx][idx_stack_bloc[k]][line_idx_in_bloc].ValuesData.push_back (v1_stack[k]);
				A[row_bloc_idx][idx_stack_bloc[k]][line_idx_in_bloc].ValuesData.push_back (v2_stack[k]);
			}
		}

		if(i<nb_rows)
		{
			it1 = M[rows_idxs[i]].begin ();

			idx_stack_elt_idx_in_line.clear();
			idx_stack_bloc.clear();
			v1_stack.clear();
			v2_stack.clear();

			line_idx_in_bloc = i/2;

			while (it1 != M[rows_idxs[i]].end())
			{
				if(pivot_columns_map[it1->first] != MINUS_ONE)		//pivot column
				{
					bloc_idx_in_row = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) / A.bloc_width();
					elt_idx_in_line = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) % A.bloc_width();

					idx_stack_elt_idx_in_line.push_back (elt_idx_in_line);
					idx_stack_bloc.push_back(bloc_idx_in_row);
					v1_stack.push_back (it1->second);
					v2_stack.push_back (0);
				}
				else
				{
					bloc_idx_in_row = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
					elt_idx_in_line = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it1->second);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
				}

				++it1;
			}

			//write into A (Right To Left)
			for (int k = idx_stack_elt_idx_in_line.size() - 1; k >= 0; --k)
			{
				A[row_bloc_idx][idx_stack_bloc[k]][line_idx_in_bloc].IndexData.push_back (idx_stack_elt_idx_in_line[k]);
				A[row_bloc_idx][idx_stack_bloc[k]][line_idx_in_bloc].ValuesData.push_back (v1_stack[k]);
				A[row_bloc_idx][idx_stack_bloc[k]][line_idx_in_bloc].ValuesData.push_back (v2_stack[k]);
			}
		}

		//Make hybrid multirows
		if(!B.acceptRowsHybrid)
			return;

		MultiLineVector<Element, Index> tmp;

		for(uint32 j=0; j<B[row_bloc_idx].size (); ++j)
		{
			for(uint16 k=0; k<B[row_bloc_idx][j].bloc_height (); ++k)
			{
				if(B[row_bloc_idx][j][k].empty ())
					continue;

				if((float)B[row_bloc_idx][j][k].size () / (float)B.bloc_width() < B.get_HYBRID_REPRESENTATION_THRESHOLD ())
					continue;

				tmp.clear();
				uint32 idx=0;

				for(uint32 p=0; p<B.bloc_width(); ++p) // p<A[i][j][k].size (); ++p)
				{
					if(idx < B[row_bloc_idx][j][k].size () && B[row_bloc_idx][j][k].IndexData[idx] == p)
					{
						tmp.ValuesData.push_back(B[row_bloc_idx][j][k].at_unchecked(0, idx));
						tmp.ValuesData.push_back(B[row_bloc_idx][j][k].at_unchecked(1, idx));
						++idx;
					}
					else
					{
						tmp.ValuesData.push_back(0);
						tmp.ValuesData.push_back(0);
					}
				}

				B[row_bloc_idx][j][k].swap(tmp);
			}
		}
	}

	/**
	 * Constructs the submatrices A, B, C, D from matrix M
	 * In the case where the index maps are not yet constructed, the functions constructs the maps
	 * implicitly; in that case, the matrices are re_initiliazed with the new corresponding dimentions
	 */

	void constructSubMatrices(SparseMatrix<Element>& M,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& C,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& D,
			bool destruct_original_matrix)
	{
		//std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

		typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

		if (!_index_maps_constructed)
		{
			commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "Warning: Maps are not yet constructed. Constructing maps and reinitializing submatrices"
					<< std::endl;

			//throw std::runtime_error ("Indexes not constructed yet! call Indexer.processMatrix() first");
			processMatrix(M);

			A = Matrix(Npiv, Npiv, Matrix::ArrangementDownTop_RightLeft, false,
					A.bloc_height(), A.bloc_width());
			A.acceptRowsHybrid = false;

			B = Matrix(Npiv, M.coldim() - Npiv,
					Matrix::ArrangementDownTop_LeftRight, true, B.bloc_height(),
					B.bloc_width()); //true: add fill blocs since they might be used
			B.acceptRowsHybrid = true;

			C = Matrix(M.rowdim() - Npiv, Npiv,
					Matrix::ArrangementDownTop_RightLeft, false,
					C.bloc_height(), C.bloc_width());
			C.acceptRowsHybrid = false;

			D = Matrix(M.rowdim() - Npiv, M.coldim() - Npiv,
					Matrix::ArrangementDownTop_LeftRight, true, D.bloc_height(),
					D.bloc_width()); //true: add fill blocs since they might be used
			D.acceptRowsHybrid = true;
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

		uint32 nb_pivs = 0;

		uint32 piv_starting_idxs[this->Npiv / B.bloc_height() + 2];
		piv_starting_idxs[0] = M.coldim();
		uint32 bloc;

		for (int i = (int) M.coldim() - 1; i >= 0; --i)
		{
			if (pivot_rows_idxs_by_entry[i] != MINUS_ONE) //non pivot row
				nb_pivs++;

			if (nb_pivs % B.bloc_height() == 0)
				piv_starting_idxs[nb_pivs / B.bloc_height()] = i;
		}

		piv_starting_idxs[this->Npiv / B.bloc_height() + 1] = 0;

omp_set_dynamic(0);
#pragma omp parallel num_threads(NB_THREADS)
		{
			uint32 rows_idxs_horizontal_bloc[B.bloc_height()];
			Index curr_vect_in_bloc = 0;

#pragma omp for schedule(dynamic) nowait
			for (bloc = 0; bloc <= this->Npiv / B.bloc_height(); ++bloc)
			{
				//MATRIX A and B
				for (int i = (int) piv_starting_idxs[bloc] - 1;
						i >= (int) piv_starting_idxs[bloc + 1]; --i)
				{
					if (pivot_rows_idxs_by_entry[i] != MINUS_ONE) //non pivot row
					{
						rows_idxs_horizontal_bloc[curr_vect_in_bloc] =
								pivot_rows_idxs_by_entry[i];
						curr_vect_in_bloc++;
					}

					if (curr_vect_in_bloc == A.bloc_height() || i == 0)
					{
						write_row_blocs_to_Left_Right_matrix(M, A, B,
								rows_idxs_horizontal_bloc, curr_vect_in_bloc,
								bloc);

						if (destruct_original_matrix)
						{
							for (uint32 j = 0; j < curr_vect_in_bloc; ++j)
							{
								M[rows_idxs_horizontal_bloc[j]].free();
							}
						}

						//row_bloc_idx++;
						curr_vect_in_bloc = 0;
					}
				}
			}
		}

		uint32 piv_starting_idxs2[C.rowdim() / B.bloc_height() + 2];
		piv_starting_idxs2[0] = M.rowdim();
		nb_pivs = 0;
		for (int i = (int) M.rowdim() - 1; i >= 0; --i)
		{
			if (non_pivot_rows_idxs[i] != MINUS_ONE) //non pivot row
				nb_pivs++;

			if (nb_pivs % B.bloc_height() == 0)
				piv_starting_idxs2[nb_pivs / B.bloc_height()] = i;
		}

		piv_starting_idxs2[C.rowdim() / B.bloc_height() + 1] = 0;


omp_set_dynamic(0);
#pragma omp parallel num_threads(NB_THREADS)
		{
			uint32 rows_idxs_horizontal_bloc[B.bloc_height()];
			Index curr_vect_in_bloc = 0;

#pragma omp for schedule(dynamic) nowait
			for (bloc = 0; bloc <= C.rowdim() / B.bloc_height(); ++bloc)
			{
				//MATRIX C and D
				for (int i = (int) piv_starting_idxs2[bloc] - 1;
						i >= (int) piv_starting_idxs2[bloc + 1]; --i) //non pivot rows
				{
					if (non_pivot_rows_idxs[i] != MINUS_ONE) //non pivot row
					{
						rows_idxs_horizontal_bloc[curr_vect_in_bloc] =
								non_pivot_rows_idxs[i];
						curr_vect_in_bloc++;
					}

					if (curr_vect_in_bloc == C.bloc_height() || i == 0)
					{
						write_row_blocs_to_Left_Right_matrix(M, C, D,
								rows_idxs_horizontal_bloc, curr_vect_in_bloc,
								bloc);

						if (destruct_original_matrix)
						{
							for (uint32 j = 0; j < curr_vect_in_bloc; ++j)
							{
								M[rows_idxs_horizontal_bloc[j]].free();
							}
						}

						curr_vect_in_bloc = 0;
					}

				}
			}
		}

	}

	

	void write_row_blocs_to_LeftMultiline_RightBloc_matrix(const SparseMatrix<Element>& M,
							        SparseMultilineMatrix<Element>& A,
								SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
								uint32 *rows_idxs,
								uint32 nb_rows,
								uint32 row_bloc_idx)
	{
		typename SparseVector<Element>::const_iterator it1, it2;

		B[row_bloc_idx].clear ();

		uint32 bloc_idx_in_row, elt_idx_in_line;

		uint32 line_idx_in_bloc;
		uint32 nb_blocs_per_dim;

		B.FirstBlocsColumIndexes[row_bloc_idx] = 0;

		nb_blocs_per_dim = (uint32) std::ceil((float)B.coldim () / B.bloc_width ());
		B[row_bloc_idx].resize (nb_blocs_per_dim);

		for(uint32 k=0; k<nb_blocs_per_dim; ++k)
			B[row_bloc_idx][k].init(B.bloc_height (), B.bloc_width ());

		uint32 row_in_multiline_A;
		
		//2. write elements to blocs
		uint32 i=0;
		for(; i< ROUND_DOWN(nb_rows, 2); i+=2)	//A.bloc_height () must be even
		{
			it1 = M[rows_idxs[i]].begin ();
			it2 = M[rows_idxs[i+1]].begin ();

			line_idx_in_bloc = i/2;

			row_in_multiline_A = //(A.rowdim () / NB_ROWS_PER_MULTILINE) + (A.rowdim () % NB_ROWS_PER_MULTILINE) - 1 - 
						row_bloc_idx * B.bloc_height () / NB_ROWS_PER_MULTILINE + line_idx_in_bloc;
			
			while (it1 != M[rows_idxs[i]].end () && it2 != M[rows_idxs[i+1]].end ())
			{
				if(it1->first < it2->first)
				{
					if(pivot_columns_map[it1->first] != MINUS_ONE)		//pivot column
					{
						//lines are inverted in this case
						A[row_in_multiline_A].IndexData.push_back (pivot_columns_map[it1->first]);
						A[row_in_multiline_A].ValuesData.push_back (it1->second);
						A[row_in_multiline_A].ValuesData.push_back (0);
						
// 						bloc_idx_in_row = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) / A.bloc_width();
// 						elt_idx_in_line = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) % A.bloc_width();
// 
// 						idx_stack_elt_idx_in_line.push_back (elt_idx_in_line);
// 						idx_stack_bloc.push_back(bloc_idx_in_row);
// 						v1_stack.push_back (it1->second);
// 						v2_stack.push_back (0);
					}
					else
					{
						bloc_idx_in_row = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
						elt_idx_in_line = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it1->second);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
					}

					++it1;
				}
				else if(it1->first > it2->first)
				{
					if(pivot_columns_map[it2->first] != MINUS_ONE)		//pivot column
					{
						A[row_in_multiline_A].IndexData.push_back (pivot_columns_map[it2->first]);
						A[row_in_multiline_A].ValuesData.push_back (0);
						A[row_in_multiline_A].ValuesData.push_back (it2->second);
						
						
// 						bloc_idx_in_row = (A.coldim () - 1 - pivot_columns_map[it2->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) / A.bloc_width();
// 						elt_idx_in_line = (A.coldim () - 1 - pivot_columns_map[it2->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) % A.bloc_width();
// 
// 						idx_stack_elt_idx_in_line.push_back (elt_idx_in_line);
// 						idx_stack_bloc.push_back(bloc_idx_in_row);
// 						v1_stack.push_back (0);
// 						v2_stack.push_back (it2->second);
					}
					else
					{
						bloc_idx_in_row = (non_pivot_columns_map[it2->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
						elt_idx_in_line = (non_pivot_columns_map[it2->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it2->second);
					}

					++it2;
				}
				else	// it1->first == it2->first
				{
					if(pivot_columns_map[it1->first] != MINUS_ONE)		//pivot column
					{
						A[row_in_multiline_A].IndexData.push_back (pivot_columns_map[it1->first]);
						A[row_in_multiline_A].ValuesData.push_back (it1->second);
						A[row_in_multiline_A].ValuesData.push_back (it2->second);
						
// 						bloc_idx_in_row = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) / A.bloc_width();
// 						elt_idx_in_line = (A.coldim () - 1 - pivot_columns_map[it1->first] - A.FirstBlocsColumIndexes[row_bloc_idx]) % A.bloc_width();
// 
// 						idx_stack_elt_idx_in_line.push_back (elt_idx_in_line);
// 						idx_stack_bloc.push_back(bloc_idx_in_row);
// 						v1_stack.push_back (it1->second);
// 						v2_stack.push_back (it2->second);
					}
					else
					{
						bloc_idx_in_row = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
						elt_idx_in_line = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it1->second);
						B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it2->second);
					}

					++it1;
					++it2;
				}
			}

			while (it1 != M[rows_idxs[i]].end())
			{
				if(pivot_columns_map[it1->first] != MINUS_ONE)		//pivot column
				{
					A[row_in_multiline_A].IndexData.push_back (pivot_columns_map[it1->first]);
					A[row_in_multiline_A].ValuesData.push_back (it1->second);
					A[row_in_multiline_A].ValuesData.push_back (0);
				}
				else
				{
					bloc_idx_in_row = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
					elt_idx_in_line = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it1->second);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
				}

				++it1;
			}

			while (it2 != M[rows_idxs[i+1]].end())
			{
				if(pivot_columns_map[it2->first] != MINUS_ONE)		//pivot column
				{
					A[row_in_multiline_A].IndexData.push_back (pivot_columns_map[it2->first]);
					A[row_in_multiline_A].ValuesData.push_back (0);
					A[row_in_multiline_A].ValuesData.push_back (it2->second);
				}
				else
				{
					bloc_idx_in_row = (non_pivot_columns_map[it2->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
					elt_idx_in_line = (non_pivot_columns_map[it2->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it2->second);
				}

				++it2;
			}
		}

		if(i<nb_rows)
		{
			it1 = M[rows_idxs[i]].begin ();
			line_idx_in_bloc = i/2;

			row_in_multiline_A = //A.rowdim () / NB_ROWS_PER_MULTILINE + A.rowdim () % NB_ROWS_PER_MULTILINE - 1 - 
					row_bloc_idx * B.bloc_height () / NB_ROWS_PER_MULTILINE + line_idx_in_bloc;
			
			while (it1 != M[rows_idxs[i]].end())
			{
				if(pivot_columns_map[it1->first] != MINUS_ONE)		//pivot column
				{
					A[row_in_multiline_A].IndexData.push_back (pivot_columns_map[it1->first]);
					A[row_in_multiline_A].ValuesData.push_back (it1->second);
					A[row_in_multiline_A].ValuesData.push_back (0);
				}
				else
				{
					bloc_idx_in_row = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
					elt_idx_in_line = (non_pivot_columns_map[it1->first] - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (it1->second);
					B[row_bloc_idx][bloc_idx_in_row][line_idx_in_bloc].ValuesData.push_back (0);
				}

				++it1;
			}
		}
	}
	
	//C multiline

	void constructSubMatrices(SparseMatrix<Element>& M,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
			SparseMultilineMatrix<Element>& C,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& D,
			bool destruct_original_matrix)
	{
		typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

		if (!_index_maps_constructed)
		{
			commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "Warning: Maps are not yet constructed. Constructing maps and reinitializing submatrices"
					<< std::endl;

			//throw std::runtime_error ("Indexes not constructed yet! call Indexer.processMatrix() first");
			processMatrix(M);

			A = Matrix(Npiv, Npiv, Matrix::ArrangementDownTop_RightLeft, false,
					A.bloc_height(), A.bloc_width());
			A.acceptRowsHybrid = false;

			B = Matrix(Npiv, M.coldim() - Npiv,
					Matrix::ArrangementDownTop_LeftRight, true, B.bloc_height(),
					B.bloc_width()); //true: add fill blocs since they might be used
			B.acceptRowsHybrid = true;

			C = SparseMultilineMatrix<Element>(M.rowdim() - Npiv, Npiv);

			D = Matrix(M.rowdim() - Npiv, M.coldim() - Npiv,
					Matrix::ArrangementDownTop_LeftRight, true, D.bloc_height(),
					D.bloc_width()); //true: add fill blocs since they might be used
			D.acceptRowsHybrid = true;
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

		uint32 nb_pivs = 0;

		uint32 piv_starting_idxs[this->Npiv / B.bloc_height() + 2];
		piv_starting_idxs[0] = M.coldim();
		uint32 bloc;

		for (int i = (int) M.coldim() - 1; i >= 0; --i)
		{
			if (pivot_rows_idxs_by_entry[i] != MINUS_ONE) //non pivot row
				nb_pivs++;

			if (nb_pivs % B.bloc_height() == 0)
				piv_starting_idxs[nb_pivs / B.bloc_height()] = i;
		}

		piv_starting_idxs[this->Npiv / B.bloc_height() + 1] = 0;

		omp_set_dynamic(0);
#pragma omp parallel num_threads(NB_THREADS)
		{
			uint32 rows_idxs_horizontal_bloc[B.bloc_height()];
			Index curr_vect_in_bloc = 0;

#pragma omp for schedule(dynamic) nowait
			for (bloc = 0; bloc <= this->Npiv / B.bloc_height(); ++bloc)
			{
				//MATRIX A and B
				for (int i = (int) piv_starting_idxs[bloc] - 1;
						i >= (int) piv_starting_idxs[bloc + 1]; --i)
				{
					if (pivot_rows_idxs_by_entry[i] != MINUS_ONE) //non pivot row
					{
						rows_idxs_horizontal_bloc[curr_vect_in_bloc] =
								pivot_rows_idxs_by_entry[i];
						curr_vect_in_bloc++;
					}

					if (curr_vect_in_bloc == A.bloc_height() || i == 0)
					{
						write_row_blocs_to_Left_Right_matrix(M, A, B,
								rows_idxs_horizontal_bloc, curr_vect_in_bloc,
								bloc);

						if (destruct_original_matrix)
						{
							for (uint32 j = 0; j < curr_vect_in_bloc; ++j)
							{
								M[rows_idxs_horizontal_bloc[j]].free();
							}
						}

						curr_vect_in_bloc = 0;
					}
				}
			}
		}

		//MATRIX C and D
		uint32 piv_starting_idxs2[C.rowdim() / B.bloc_height() + 2];
		piv_starting_idxs2[0] = M.rowdim();
		nb_pivs = 0;
		for (int i = (int) M.rowdim() - 1; i >= 0; --i)
		{
			if (non_pivot_rows_idxs[i] != MINUS_ONE) //non pivot row
				nb_pivs++;

			if (nb_pivs % B.bloc_height() == 0)
				piv_starting_idxs2[nb_pivs / B.bloc_height()] = i;
		}

		piv_starting_idxs2[C.rowdim() / B.bloc_height() + 1] = 0;

		omp_set_dynamic(0);
#pragma omp parallel num_threads(NB_THREADS)
		{

			uint32 rows_idxs_horizontal_bloc[B.bloc_height()];
			Index curr_vect_in_bloc = 0;

#pragma omp for schedule(dynamic) nowait
			for (bloc = 0; bloc <= C.rowdim() / B.bloc_height(); ++bloc)
			{
				//MATRIX C and D
				for (int i = (int) piv_starting_idxs2[bloc] - 1;
						i >= (int) piv_starting_idxs2[bloc + 1]; --i) //non pivot rows
				{
					if (non_pivot_rows_idxs[i] != MINUS_ONE) //non pivot row
					{
						rows_idxs_horizontal_bloc[curr_vect_in_bloc] =
								non_pivot_rows_idxs[i];
						curr_vect_in_bloc++;
					}

					if (curr_vect_in_bloc == D.bloc_height() || i == 0)
					{
						write_row_blocs_to_LeftMultiline_RightBloc_matrix(M, C,
								D, rows_idxs_horizontal_bloc, curr_vect_in_bloc,
								bloc);

						if (destruct_original_matrix)
						{
							for (uint32 j = 0; j < curr_vect_in_bloc; ++j)
							{
								M[rows_idxs_horizontal_bloc[j]].free();
							}
						}

						curr_vect_in_bloc = 0;
					}

				}
			}
		}

	}
	

	//A and C multiline

	void constructSubMatrices(SparseMatrix<Element>& M,
			SparseMultilineMatrix<Element>& A,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
			SparseMultilineMatrix<Element>& C,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& D,
			bool destruct_original_matrix)
	{
		typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

		if (!_index_maps_constructed)
		{
			commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "Warning: Maps are not yet constructed. Constructing maps and reinitializing submatrices"
					<< std::endl;

			//throw std::runtime_error ("Indexes not constructed yet! call Indexer.processMatrix() first");
			processMatrix(M);

			A = SparseMultilineMatrix<Element>(Npiv, Npiv);

			B = Matrix(Npiv, M.coldim() - Npiv,
					Matrix::ArrangementDownTop_LeftRight, true, B.bloc_height(),
					B.bloc_width()); //true: add fill blocs since they might be used
			B.acceptRowsHybrid = true;

			C = SparseMultilineMatrix<Element>(M.rowdim() - Npiv, Npiv);

			D = Matrix(M.rowdim() - Npiv, M.coldim() - Npiv,
					Matrix::ArrangementDownTop_LeftRight, true, D.bloc_height(),
					D.bloc_width()); //true: add fill blocs since they might be used
			D.acceptRowsHybrid = true;
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

		//A & B
		uint32 nb_pivs = 0;
		uint32 piv_starting_idxs[this->Npiv / B.bloc_height() + 2];
		piv_starting_idxs[0] = M.coldim();
		uint32 bloc;

		for (int i = (int) M.coldim() - 1; i >= 0; --i)
		{
			if (pivot_rows_idxs_by_entry[i] != MINUS_ONE) //non pivot row
				nb_pivs++;

			if (nb_pivs % B.bloc_height() == 0)
				piv_starting_idxs[nb_pivs / B.bloc_height()] = i;
		}

		piv_starting_idxs[this->Npiv / B.bloc_height() + 1] = 0;

		omp_set_dynamic(0);
#pragma omp parallel num_threads(NB_THREADS)
		{

			uint32 rows_idxs_horizontal_bloc[B.bloc_height()];
			Index curr_vect_in_bloc = 0;

#pragma omp for schedule(dynamic) nowait
			for (bloc = 0; bloc <= this->Npiv / B.bloc_height(); ++bloc)
			{
				//MATRIX A and B
				for (int i = (int) piv_starting_idxs[bloc] - 1;
						i >= (int) piv_starting_idxs[bloc + 1]; --i)
				{
					if (pivot_rows_idxs_by_entry[i] != MINUS_ONE) //non pivot row
					{
						rows_idxs_horizontal_bloc[curr_vect_in_bloc] =
								pivot_rows_idxs_by_entry[i];
						curr_vect_in_bloc++;
					}

					if (curr_vect_in_bloc == B.bloc_height() || i == 0)
					{
						write_row_blocs_to_LeftMultiline_RightBloc_matrix(M, A,
								B, rows_idxs_horizontal_bloc, curr_vect_in_bloc,
								bloc);

						if (destruct_original_matrix)
						{
							for (uint32 j = 0; j < curr_vect_in_bloc; ++j)
							{
								M[rows_idxs_horizontal_bloc[j]].free();
							}
						}

						curr_vect_in_bloc = 0;
					}
				}
			}
		}

		//MATRIX C and D
		//MATRIX C and D
		uint32 piv_starting_idxs2[C.rowdim() / B.bloc_height() + 2];
		piv_starting_idxs2[0] = M.rowdim();
		nb_pivs = 0;
		for (int i = (int) M.rowdim() - 1; i >= 0; --i)
		{
			if (non_pivot_rows_idxs[i] != MINUS_ONE) //non pivot row
				nb_pivs++;

			if (nb_pivs % B.bloc_height() == 0)
				piv_starting_idxs2[nb_pivs / B.bloc_height()] = i;
		}

		piv_starting_idxs2[C.rowdim() / B.bloc_height() + 1] = 0;

		omp_set_dynamic(0);
#pragma omp parallel num_threads(NB_THREADS)
		{

			uint32 rows_idxs_horizontal_bloc[B.bloc_height()];
			Index curr_vect_in_bloc = 0;

#pragma omp for schedule(dynamic) nowait
			for (bloc = 0; bloc <= C.rowdim() / B.bloc_height(); ++bloc)
			{
				//MATRIX C and D
				for (int i = (int) piv_starting_idxs2[bloc] - 1;
						i >= (int) piv_starting_idxs2[bloc + 1]; --i) //non pivot rows
				{
					if (non_pivot_rows_idxs[i] != MINUS_ONE) //non pivot row
					{
						rows_idxs_horizontal_bloc[curr_vect_in_bloc] =
								non_pivot_rows_idxs[i];
						curr_vect_in_bloc++;
					}

					if (curr_vect_in_bloc == D.bloc_height() || i == 0)
					{
						write_row_blocs_to_LeftMultiline_RightBloc_matrix(M, C,
								D, rows_idxs_horizontal_bloc, curr_vect_in_bloc,
								bloc);

						if (destruct_original_matrix)
						{
							for (uint32 j = 0; j < curr_vect_in_bloc; ++j)
							{
								M[rows_idxs_horizontal_bloc[j]].free();
							}
						}

						curr_vect_in_bloc = 0;
					}

				}
			}
		}

	}



	//inneficient computations

	void constructSubMatrices(SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
			SparseMultilineMatrix<Element>& D,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B1,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B2,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& D1,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& D2,
			bool destruct_original_matrix)
	{
		typedef SparseBlocMatrix<SparseMultilineBloc<uint16, Index> > Matrix;

		B1 = Matrix(B.rowdim(), Npiv,
				Matrix::ArrangementDownTop_RightLeft,
				false, B.bloc_height(), B.bloc_width());
		B1.acceptRowsHybrid = false;

		B2 = Matrix(B.rowdim(), B.coldim() - Npiv,
				Matrix::ArrangementDownTop_LeftRight,
				true, B.bloc_height(), B.bloc_width()); //true: add fill blocs since they might be used
		B2.acceptRowsHybrid = true;

		D1 = Matrix(Npiv, Npiv,
				Matrix::ArrangementDownTop_RightLeft,
				false, B.bloc_height(), B.bloc_width());
		D1.acceptRowsHybrid = false;

		D2 = Matrix(Npiv, D.coldim() - Npiv,
				Matrix::ArrangementDownTop_LeftRight,
				true, B.bloc_height(), B.bloc_width()); //true: add fill blocs since they might be used
		D2.acceptRowsHybrid = true;


//		MatrixUtils::show_mem_usage("before copy");
//		commentator.start("copy matrix");
//			SparseMatrix<Element> _tmp_B (B.rowdim(), B.coldim());
//			MatrixUtils::copy(B, _tmp_B, destruct_original_matrix);
//		commentator.stop(MSG_DONE);
//		MatrixUtils::show_mem_usage("after copy");

		SparseMatrix<Element> _tmp_B (B.bloc_height(), B.coldim());

		uint8 *vectors_to_free_B;
		if(destruct_original_matrix)	//free data on the go
		{
			posix_memalign((void**)&vectors_to_free_B, 16, (B.rowdim() / NB_ROWS_PER_MULTILINE + 1)  * sizeof(uint8));
			Level1Ops::memsetToZero(&vectors_to_free_B, 1, B.rowdim() / NB_ROWS_PER_MULTILINE + 1);
		}

		//B1 B2
		uint32 row_bloc_idx = 0;
		Index curr_vect_in_bloc = 0;

		uint32 rows_idxs_horizontal_bloc[B1.bloc_height ()];

		//MATRIX B
		for (int i = B.rowdim () - 1; i >= 0; --i)
		{
			//-------------------copy B.bloc_height() at a time ------------------------------------
			//for each bloc of A going horizontally
			uint32 blc_row_idx = (B.rowdim() - 1 - i) / B.bloc_height();

			uint32 row_idx_in_blc = ((B.rowdim() - 1 - i) % B.bloc_height()) / NB_ROWS_PER_MULTILINE;
			uint32 bloc_start_idx;
			Element val;
			Index idx;
			uint32 line = (B.rowdim() - 1 - i) % NB_ROWS_PER_MULTILINE;

			_tmp_B[curr_vect_in_bloc].clear ();

			for(int j=0; j<(int)B[blc_row_idx].size (); ++j)
			{
				if(B[blc_row_idx][j].empty ())
					continue;

				if(B[blc_row_idx][j][row_idx_in_blc].empty ())
					continue;

				bloc_start_idx = B.FirstBlocsColumIndexes[blc_row_idx] + (B.bloc_width() * j);

				if(B[blc_row_idx][j][row_idx_in_blc].is_sparse ())
				{
					for (int p = 0; p<(int)B[blc_row_idx][j][row_idx_in_blc].size(); ++p)
					{
						idx = B[blc_row_idx][j][row_idx_in_blc].IndexData[p];
						val = B[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);

						if(val != 0 )
							_tmp_B[curr_vect_in_bloc].push_back(typename SparseVector<Element>::value_type(bloc_start_idx + idx, val));
					}
				}
				else
				{
					for (int p = 0; p < B.bloc_width (); ++p)
					{
						val = B[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);

						if(val != 0)
						{
							_tmp_B[curr_vect_in_bloc].push_back(typename SparseVector<Element>::value_type(bloc_start_idx + p, val));
						}
					}
				}

				if(destruct_original_matrix)
				{
					if(vectors_to_free_B[(B.rowdim() - 1 - i)/ NB_ROWS_PER_MULTILINE] == 1)
						B[blc_row_idx][j][row_idx_in_blc].free ();
				}
			}

			if(destruct_original_matrix)
			{
				if(vectors_to_free_B[(B.rowdim() - 1 - i)/ NB_ROWS_PER_MULTILINE] == 0)
					vectors_to_free_B[(B.rowdim() - 1 - i)/ NB_ROWS_PER_MULTILINE] = 1;
			}


			//--------------------------------------------------------------------------------------
			rows_idxs_horizontal_bloc[curr_vect_in_bloc] = curr_vect_in_bloc;
			curr_vect_in_bloc++;

			if(curr_vect_in_bloc == B1.bloc_height ())
			{
				write_row_blocs_to_Left_Right_matrix(_tmp_B, B1, B2, rows_idxs_horizontal_bloc, curr_vect_in_bloc, row_bloc_idx);

				row_bloc_idx++;
				curr_vect_in_bloc = 0;
			}
		}

		//write the left rows of A and B
		if(curr_vect_in_bloc != 0)
			write_row_blocs_to_Left_Right_matrix(_tmp_B, B1, B2, rows_idxs_horizontal_bloc, curr_vect_in_bloc, row_bloc_idx);

		if(destruct_original_matrix)
		{
			for(uint32 j=0; j<curr_vect_in_bloc; ++j)
			{
				_tmp_B[rows_idxs_horizontal_bloc[j]].free ();
			}

			free(vectors_to_free_B);
		}

		//D1 D2
		//TODO: unefficient code
		//TODO: free D on the go

		uint8 *vectors_to_free_D;
		if(destruct_original_matrix)	//free data on the go
		{
			posix_memalign((void**)&vectors_to_free_D, 16, (D.rowdim() / NB_ROWS_PER_MULTILINE + 1)  * sizeof(uint8));
			Level1Ops::memsetToZero(&vectors_to_free_D, 1, D.rowdim() / NB_ROWS_PER_MULTILINE + 1);
		}

		SparseMatrix<Element> _tmp_D (B.bloc_height(), D.coldim());
		MultiLineVector<Element, uint32> *rowD;
		uint32 idx, val;

		curr_vect_in_bloc = 0;
		row_bloc_idx = 0;
		for (int i = D.coldim () - 1; i >= 0 ; --i)
		{
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)			//non pivot row
				continue;

			rowD = &(D[pivot_rows_idxs_by_entry[i]/NB_ROWS_PER_MULTILINE]);

			//copy from multiline to sparse vector
			if(rowD->is_sparse (D.coldim ()))
				for(uint32 j=0; j<rowD->size (); ++j)
				{
					idx = rowD->IndexData[j];
					val = rowD->at_unchecked(pivot_rows_idxs_by_entry[i] % NB_ROWS_PER_MULTILINE, j);

					if(val != 0)
						_tmp_D[curr_vect_in_bloc].push_back(typename SparseMatrix<Element>::Row::value_type(idx, val));
				}
			else
				for(uint32 j=0; j<D[pivot_rows_idxs_by_entry[i]/NB_ROWS_PER_MULTILINE].size (); ++j)
				{
					val = rowD->at_unchecked(pivot_rows_idxs_by_entry[i] % NB_ROWS_PER_MULTILINE, j);

					if(val != 0)
						_tmp_D[curr_vect_in_bloc].push_back(typename SparseMatrix<Element>::Row::value_type(j, val));
				}

			if(destruct_original_matrix)
			{
				if(vectors_to_free_D[pivot_rows_idxs_by_entry[i]/NB_ROWS_PER_MULTILINE] == 0)
					vectors_to_free_D[pivot_rows_idxs_by_entry[i]/NB_ROWS_PER_MULTILINE] = 1;
				else
					rowD->free ();
			}

			rows_idxs_horizontal_bloc[curr_vect_in_bloc] = curr_vect_in_bloc;
			curr_vect_in_bloc++;

			if(curr_vect_in_bloc == B1.bloc_height ())
			{
				write_row_blocs_to_Left_Right_matrix(_tmp_D, D1, D2, rows_idxs_horizontal_bloc, curr_vect_in_bloc, row_bloc_idx);
				for(uint32 j=0; j<curr_vect_in_bloc; ++j)
				{
					_tmp_D[j].clear ();
				}

				row_bloc_idx++;
				curr_vect_in_bloc = 0;
			}
		}


		if(curr_vect_in_bloc != 0)
			write_row_blocs_to_Left_Right_matrix(_tmp_D, D1, D2, rows_idxs_horizontal_bloc, curr_vect_in_bloc, row_bloc_idx);

		for(uint32 j=0; j<curr_vect_in_bloc; ++j)
		{
			_tmp_D[j].clear ();
		}

		if(destruct_original_matrix)
		{
			free(vectors_to_free_D);
			D.free();
			B.free();
		}
	}


	template <typename Multiline, typename SparseRow>
	void push_rowsAB_to_rowM(const Multiline& rowA, const Multiline& rowB, const uint16 line, SparseRow& rowM)
	{
		uint32 currA, currB, idxA, idxB;
		currA = 0;
		currB = 0;

		while(currA < rowA.size () && currB < rowB.size ())
		{
			idxA = rowA.IndexData[currA];
			idxB = rowB.IndexData[currB];

			//cout << "A " << idxA << " B " << idxB << endl;

			if(rowA.at_unchecked(line, currA) == 0)
			{
				++currA;
				continue;
			}

			if(rowB.at_unchecked(line, currB) == 0)
			{
				++currB;
				continue;
			}

			if(pivot_columns_rev_map[idxA] < non_pivot_columns_rev_map[idxB])
			{
				rowM.push_back (typename SparseRow::value_type(pivot_columns_rev_map[idxA], rowA.at_unchecked(line, currA)));
				++currA;
			}
			else if(pivot_columns_rev_map[idxA] > non_pivot_columns_rev_map[idxB])
			{
				rowM.push_back (typename SparseRow::value_type(non_pivot_columns_rev_map[idxB], rowB.at_unchecked(line, currB)));
				++currB;
			}
			else
				throw "ERROR IN INDEX";
		}

		while(currA < rowA.size ())
		{
			idxA = rowA.IndexData[currA];

			if(rowA.at_unchecked(line, currA) == 0)
			{
				++currA;
				continue;
			}

			rowM.push_back (typename SparseRow::value_type(pivot_columns_rev_map[idxA], rowA.at_unchecked(line, currA)));
			++currA;
		}

		while(currB < rowB.size ())
		{
			idxB = rowB.IndexData[currB];

			if(rowB.at_unchecked(line, currB) == 0)
			{
				++currB;
				continue;
			}

			rowM.push_back (typename SparseRow::value_type(non_pivot_columns_rev_map[idxB], rowB.at_unchecked(line, currB)));
			++currB;
		}
	}


	template <typename SparseRow>
	void push_rowsAB_to_rowM(const SparseRow& rowA, const SparseRow& rowB, SparseRow& rowM)
	{
		typename SparseRow::const_iterator it1, it2;
		it1 = rowA.begin ();
		it2 = rowB.begin ();

		while(it1 != rowA.end () && it2 != rowB.end ())
		{
			if(pivot_columns_rev_map[it1->first] < non_pivot_columns_rev_map[it2->first])
			{
				rowM.push_back (typename SparseRow::value_type(pivot_columns_rev_map[it1->first], it1->second));
				++it1;
			}
			else if(pivot_columns_rev_map[it1->first] > non_pivot_columns_rev_map[it2->first])
			{
				rowM.push_back (typename SparseRow::value_type(non_pivot_columns_rev_map[it2->first], it2->second));
				++it2;
			}
			else
				throw "ERROR IN INDEX";
		}

		while(it1 != rowA.end ())
		{
			rowM.push_back (typename SparseRow::value_type(pivot_columns_rev_map[it1->first], it1->second));
			++it1;
		}

		while(it2 != rowB.end ())
		{
			rowM.push_back (typename SparseRow::value_type(non_pivot_columns_rev_map[it2->first], it2->second));
			++it2;
		}
	}



	void reconstructMatrix(SparseMatrix<Element>& M, SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
									  SparseMultilineMatrix<Element>& D, bool free_matrices = false)
	{
		SparseBlocMatrix<SparseMultilineBloc<Element, Index> > _A_dummy;
		return reconstructMatrix(M, _A_dummy, B, D, free_matrices, true, false);
	}




	void reconstructMatrix(SparseMatrix<Element>& M,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
			SparseMultilineMatrix<Element>& D,
			bool free_matrices = false,
			bool A_is_null = false,
			bool D_is_null = false)
	{
		if(!A_is_null)
		{
			check_equal_or_raise_exception(A.rowdim(), B.rowdim());
			check_equal_or_raise_exception(A.rowdim(), A.coldim());
		}

		if(!D_is_null)
			check_equal_or_raise_exception(D.coldim(), B.coldim());


		uint8 *vectors_to_free_AB, *vectors_to_free_D;
		if(free_matrices)	//free data on the go
		{
			posix_memalign((void**)&vectors_to_free_AB, 16, (B.rowdim() / NB_ROWS_PER_MULTILINE + 1)  * sizeof(uint8));
			Level1Ops::memsetToZero(&vectors_to_free_AB, 1, B.rowdim() / NB_ROWS_PER_MULTILINE + 1);

			posix_memalign((void**)&vectors_to_free_D, 16, (D.rowdim() / NB_ROWS_PER_MULTILINE + 1)  * sizeof(uint8));
			Level1Ops::memsetToZero(&vectors_to_free_D, 1, D.rowdim() / NB_ROWS_PER_MULTILINE + 1);
		}

		uint32 new_piv_to_row[M.coldim()];	//To which piv in the final matrix this row it assigned
		uint32 new_piv = 0;
		for(uint32 i=0; i < this->coldim; i++)
		{
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			new_piv_to_row[i] = new_piv;
			new_piv++;
		}


omp_set_dynamic(0);
#pragma omp parallel num_threads(NB_THREADS)
		{
		typename SparseMatrix<Element>::Row *row_M;
		typename SparseMultilineMatrix<Element>::Row *rowD;
		SparseVector<Element> _tmpA, _tmpB;
		uint32 local_piv;
		uint16 line;

#pragma omp for schedule(dynamic)
		for(uint32 i=0; i < this->coldim; i++)
		{
			//TODO: handle successive rows in same multiline
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			local_piv = new_piv_to_row[i];

			if(this->pivot_rows_idxs_by_entry[i] >= this->Npiv)		//D
			{
				rowD = &(D[(this->pivot_rows_idxs_by_entry[i] - this->Npiv)/NB_ROWS_PER_MULTILINE]);
				line = (this->pivot_rows_idxs_by_entry[i] - this->Npiv)%NB_ROWS_PER_MULTILINE;

				row_M = &(M[local_piv]);
				row_M->free ();

				uint32 idxD;
				if (rowD->is_sparse (D.coldim()))
				{
					for(uint32 j=0; j<rowD->size (); ++j)
					{
						idxD = rowD->IndexData[j];
						if(rowD->at_unchecked(line, j) == 0)
							continue;

						row_M->push_back (typename SparseMatrix<Element>::Row::value_type(
											non_pivot_columns_rev_map[idxD],
											rowD->at_unchecked(line, j)));
					}
				}
				else
				{
					for(uint32 j=0; j<D.coldim (); ++j)
					{
						if(rowD->at_unchecked(line, j) == 0)
							continue;

						row_M->push_back (typename SparseMatrix<Element>::Row::value_type(
											non_pivot_columns_rev_map[j],
											rowD->at_unchecked(line, j)));
					}
				}

				if(free_matrices)
				{
					if(vectors_to_free_D[(this->pivot_rows_idxs_by_entry[i] - this->Npiv)/NB_ROWS_PER_MULTILINE] == 0)
					{
						vectors_to_free_D[(this->pivot_rows_idxs_by_entry[i] - this->Npiv)/NB_ROWS_PER_MULTILINE] = 1;	//mark it for deletion
					}
					else	//handled all the rows in the multiline, FREE now
					{
						rowD->free ();
					}
				}
			}
			else															//A & B
			{
				_tmpA.clear ();
				_tmpB.clear ();

				//for each bloc of A going horizontally
				uint32 blc_row_idx = (B.rowdim() - 1
						- this->pivot_rows_idxs_by_entry[i]) / B.bloc_height();

				uint32 row_idx_in_blc = ((B.rowdim() - 1
						- this->pivot_rows_idxs_by_entry[i]) % B.bloc_height())
						/ NB_ROWS_PER_MULTILINE;

				uint32 bloc_start_idx;
				uint16 val;
				uint32 idx;

				line = (B.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i]) % NB_ROWS_PER_MULTILINE;

				if(!A_is_null)
				{
					for(int j=A[blc_row_idx].size () - 1; j>=0; --j)
					{
						if(A[blc_row_idx][j].empty ())
							continue;

						if(A[blc_row_idx][j][row_idx_in_blc].empty ())
							continue;

						bloc_start_idx = A.coldim() - 1 - A.FirstBlocsColumIndexes[blc_row_idx] - A.bloc_width() * j;

						if(A[blc_row_idx][j][row_idx_in_blc].is_sparse ())
						{
							for (int p = A[blc_row_idx][j][row_idx_in_blc].size() - 1; p >= 0; --p)
							{
								idx = A[blc_row_idx][j][row_idx_in_blc].IndexData[p];
								val = A[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);

								if(val != 0 )
									_tmpA.push_back(typename SparseVector<Element>::value_type(bloc_start_idx - idx, val));
							}
						}
						else
						{
							for (int p = A[blc_row_idx][j][row_idx_in_blc].size() - 1; p >= 0; --p)
							{
								val = A[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);

								if(val != 0 )
									_tmpA.push_back(typename SparseVector<Element>::value_type(bloc_start_idx - p, val));
							}
						}

						if(free_matrices)
						{
							if(vectors_to_free_AB[(B.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] == 1)
							{
								A[blc_row_idx][j][row_idx_in_blc].free ();
							}
						}
					}
				}
				else	//add the identity row
				{
					_tmpA.push_back(typename SparseVector<Element>::value_type(pivot_columns_map[i], 1));
				}

				for(int j=0; j<(int)B[blc_row_idx].size (); ++j)
				{
					if(B[blc_row_idx][j].empty ())
						continue;

					if(B[blc_row_idx][j][row_idx_in_blc].empty ())
						continue;

					bloc_start_idx = B.FirstBlocsColumIndexes[blc_row_idx] + (B.bloc_width() * j);

					if(B[blc_row_idx][j][row_idx_in_blc].is_sparse ())
					{
						for (int p = 0; p<(int)B[blc_row_idx][j][row_idx_in_blc].size(); ++p)
						{
							idx = B[blc_row_idx][j][row_idx_in_blc].IndexData[p];
							val = B[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);

							if(val != 0 )
								_tmpB.push_back(typename SparseVector<Element>::value_type(bloc_start_idx + idx, val));
						}
					}
					else
					{
						for (int p = 0; p < B.bloc_width (); ++p)
						{
							val = B[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);

							if(val != 0)
							{
								_tmpB.push_back(typename SparseVector<Element>::value_type(bloc_start_idx + p, val));
							}
						}
					}

					if(free_matrices)
					{
						if(vectors_to_free_AB[(B.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] == 1)
						{
							B[blc_row_idx][j][row_idx_in_blc].free ();
						}
					}
				}

				if(free_matrices)
				{
					if(vectors_to_free_AB[(B.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] == 0)
					{
						vectors_to_free_AB[(B.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] = 1;	//mark it for deletion
					}
				}

				row_M = &(M[local_piv]);
				row_M->free ();
				
				push_rowsAB_to_rowM(_tmpA, _tmpB, *row_M);
			}
		}
	}


		if(free_matrices)
		{
			A.free ();
			B.free ();
			D.free ();

			//free matrix blocs here
			free(vectors_to_free_AB);
			free(vectors_to_free_D);
		}
	}



	void reconstructMatrix(SparseMatrix<Element>& M,
			SparseMultilineMatrix<Element>& A,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
			SparseMultilineMatrix<Element>& D,
			bool free_matrices = false)
	{
		check_equal_or_raise_exception(A.rowdim(), B.rowdim());
		check_equal_or_raise_exception(A.rowdim(), A.coldim());
		check_equal_or_raise_exception(D.coldim(), B.coldim());

		uint8 *vectors_to_free_AB, *vectors_to_free_D;
		if(free_matrices)	//free data on the go
		{
			posix_memalign((void**)&vectors_to_free_AB, 16, (B.rowdim() / NB_ROWS_PER_MULTILINE + 1)  * sizeof(uint8));
			Level1Ops::memsetToZero(&vectors_to_free_AB, 1, B.rowdim() / NB_ROWS_PER_MULTILINE + 1);

			posix_memalign((void**)&vectors_to_free_D, 16, (D.rowdim() / NB_ROWS_PER_MULTILINE + 1)  * sizeof(uint8));
			Level1Ops::memsetToZero(&vectors_to_free_D, 1, D.rowdim() / NB_ROWS_PER_MULTILINE + 1);
		}

uint32 new_piv_to_row[M.coldim()];	//To which piv in the final matrix this row it assigned
		uint32 new_piv = 0;
		for(uint32 i=0; i < this->coldim; i++)
		{
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			new_piv_to_row[i] = new_piv;
			new_piv++;
		}


omp_set_dynamic(0);
#pragma omp parallel num_threads(NB_THREADS)
		{
		typename SparseMatrix<Element>::Row *row_M;
		typename SparseMultilineMatrix<Element>::Row *rowD, *rowA;
		SparseVector<Element> _tmpA, _tmpB;
		uint32 local_piv;
		uint16 line;


#pragma omp for schedule(dynamic)
		for(uint32 i=0; i < this->coldim; i++)
		{
			//TODO: handle successive rows in same multiline
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			local_piv = new_piv_to_row[i];

			if(this->pivot_rows_idxs_by_entry[i] >= this->Npiv)		//D
			{
				rowD = &(D[(this->pivot_rows_idxs_by_entry[i] - this->Npiv)/NB_ROWS_PER_MULTILINE]);
				line = (this->pivot_rows_idxs_by_entry[i] - this->Npiv)%NB_ROWS_PER_MULTILINE;

				row_M = &(M[local_piv]);
				row_M->clear ();

				uint32 idxD;
				if (rowD->is_sparse (D.coldim()))
				{
					for(uint32 j=0; j<rowD->size (); ++j)
					{
						idxD = rowD->IndexData[j];
						if(rowD->at_unchecked(line, j) == 0)
							continue;

						row_M->push_back (typename SparseMatrix<Element>::Row::value_type(
											non_pivot_columns_rev_map[idxD],
											rowD->at_unchecked(line, j)));
					}
				}
				else
				{
					for(uint32 j=0; j<D.coldim (); ++j)
					{
						if(rowD->at_unchecked(line, j) == 0)
							continue;

						row_M->push_back (typename SparseMatrix<Element>::Row::value_type(
											non_pivot_columns_rev_map[j],
											rowD->at_unchecked(line, j)));
					}
				}

				if(free_matrices)
				{
					if(vectors_to_free_D[(this->pivot_rows_idxs_by_entry[i] - this->Npiv)/NB_ROWS_PER_MULTILINE] == 0)
					{
						vectors_to_free_D[(this->pivot_rows_idxs_by_entry[i] - this->Npiv)/NB_ROWS_PER_MULTILINE] = 1;	//mark it for deletion
					}
					else	//handled all the rows in the multiline, FREE now
					{
						rowD->free ();
					}
				}
			}
			else															//A & B
			{
				_tmpA.clear ();
				_tmpB.clear ();

				//for each bloc of A going horizontally
				uint32 blc_row_idx = (B.rowdim() - 1
						- this->pivot_rows_idxs_by_entry[i]) / B.bloc_height();

				uint32 row_idx_in_blc = ((B.rowdim() - 1
						- this->pivot_rows_idxs_by_entry[i]) % B.bloc_height())
						/ NB_ROWS_PER_MULTILINE;

				uint32 bloc_start_idx;
				uint16 val;
				uint32 idx;



				line = (B.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i]) % NB_ROWS_PER_MULTILINE;
				rowA = &(A[(A.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i]) / NB_ROWS_PER_MULTILINE]);
				//line = (this->pivot_rows_idxs_by_entry[i] - this->Npiv)%NB_ROWS_PER_MULTILINE;

				if(rowA->is_sparse (A.coldim ()))
				{
					for(uint32 j=0; j<rowA->size (); ++j)
					{
						idx = rowA->IndexData[j];
						val = rowA->at_unchecked (line, j);

						if(val != 0)
							_tmpA.push_back(typename SparseVector<Element>::value_type(idx, val));
					}
				}
				else
				{
					for(uint32 j=0; j<rowA->size (); ++j)
					{
						val = rowA->at_unchecked (line, j);

						if(val != 0)
							_tmpA.push_back(typename SparseVector<Element>::value_type(j, val));
					}
				}


				for(int j=0; j<(int)B[blc_row_idx].size (); ++j)
				{
					if(B[blc_row_idx][j].empty ())
						continue;

					if(B[blc_row_idx][j][row_idx_in_blc].empty ())
						continue;

					bloc_start_idx = B.FirstBlocsColumIndexes[blc_row_idx] + (B.bloc_width() * j);

					if(B[blc_row_idx][j][row_idx_in_blc].is_sparse ())
					{
						for (int p = 0; p<(int)B[blc_row_idx][j][row_idx_in_blc].size(); ++p)
						{
							idx = B[blc_row_idx][j][row_idx_in_blc].IndexData[p];
							val = B[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);

							if(val != 0 )
								_tmpB.push_back(typename SparseVector<Element>::value_type(bloc_start_idx + idx, val));
						}
					}
					else
					{
						for (int p = 0; p < B.bloc_width (); ++p)
						{
							val = B[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);

							if(val != 0)
							{
								_tmpB.push_back(typename SparseVector<Element>::value_type(bloc_start_idx + p, val));
							}
						}
					}

					if(free_matrices)
					{
						if(vectors_to_free_AB[(B.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] == 1)
						{
							B[blc_row_idx][j][row_idx_in_blc].free ();
						}
					}
				}

				if(free_matrices)
				{
					if(vectors_to_free_AB[(B.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] == 0)
					{
						vectors_to_free_AB[(B.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] = 1;	//mark it for deletion
					}
				}

				row_M = &(M[local_piv]);
				row_M->clear ();

				push_rowsAB_to_rowM(_tmpA, _tmpB, *row_M);
			}

		}

		}


		if(free_matrices)
		{
			A.free ();
			B.free ();
			D.free ();

			//free matrix blocs here
			free(vectors_to_free_AB);
			free(vectors_to_free_D);
		}
	}




	void reconstructMatrix(SparseMatrix<Element>& M,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B,
			bool free_matrices = false)
	{
		SparseBlocMatrix<SparseMultilineBloc<Element, Index> > _A_dummy;
		SparseMultilineMatrix<Element> _D_dummy;
		return reconstructMatrix(M, _A_dummy, B, _D_dummy, free_matrices, true, true);
	}


	void reconstructMatrix(SparseMatrix<Element>& M,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B2,
			SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& D2,
			bool free_matrices = false)
	{

		//check_equal_or_raise_exception(M.rowdim(), B2.rowdim() + D2.rowdim()); True only for Full Rank matrices

		typename SparseMatrix<Element>::Row *row_M;
		SparseVector<Element> _tmp_Left, _tmp_Right;

		uint8 *vectors_to_free_B2, *vectors_to_free_D2;
		if(free_matrices)	//free data on the go
		{
			posix_memalign((void**)&vectors_to_free_B2, 16, (B2.rowdim() / NB_ROWS_PER_MULTILINE + 1)  * sizeof(uint8));
			Level1Ops::memsetToZero(&vectors_to_free_B2, 1, B2.rowdim() / NB_ROWS_PER_MULTILINE + 1);

			posix_memalign((void**)&vectors_to_free_D2, 16, (D2.rowdim() / NB_ROWS_PER_MULTILINE + 1)  * sizeof(uint8));
			Level1Ops::memsetToZero(&vectors_to_free_D2, 1, D2.rowdim() / NB_ROWS_PER_MULTILINE + 1);
		}

		uint32 piv=0;
		uint16 line;

		//for each bloc of A going horizontally
		for(uint32 i=0; i < this->coldim; i++)
		{
			//TODO: handle successive rows in same multiline
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			row_M = &(M[piv]);
			row_M->free ();

			if(this->pivot_rows_idxs_by_entry[i] >= this->Npiv)		//D
			{
				_tmp_Left.clear ();
				_tmp_Right.clear ();

				//for each bloc of A going horizontally
				uint32 blc_row_idx = (D2.rowdim() - 1
						- (this->pivot_rows_idxs_by_entry[i] - this->Npiv)) / D2.bloc_height();

				uint32 row_idx_in_blc = ((D2.rowdim() - 1
						- (this->pivot_rows_idxs_by_entry[i] - this->Npiv)) % D2.bloc_height())
						/ NB_ROWS_PER_MULTILINE;

				uint32 bloc_start_idx;
				uint16 val;
				uint32 idx;

				line = (D2.rowdim() - 1 - (this->pivot_rows_idxs_by_entry[i] - this->Npiv)) % NB_ROWS_PER_MULTILINE;

				row_M->push_back (typename SparseMatrix<Element>::Row::value_type(i, 1));

				//D
				for(int j=0; j<(int)D2[blc_row_idx].size (); ++j)
				{
					if(D2[blc_row_idx][j].empty ())
						continue;

					if(D2[blc_row_idx][j][row_idx_in_blc].empty ())
						continue;

					bloc_start_idx = D2.FirstBlocsColumIndexes[blc_row_idx] + (D2.bloc_width() * j);

					if(D2[blc_row_idx][j][row_idx_in_blc].is_sparse ())
					{
						for (int p = 0; p<(int)D2[blc_row_idx][j][row_idx_in_blc].size(); ++p)
						{
							idx = D2[blc_row_idx][j][row_idx_in_blc].IndexData[p];
							val = D2[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);
							if(val != 0 )
								_tmp_Right.push_back(typename SparseVector<Element>::value_type(bloc_start_idx + idx, val));
						}
					}
					else
					{
						for (int p = 0; p < D2.bloc_width (); ++p)
						{
							val = D2[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);
							if(val != 0)
								_tmp_Right.push_back(typename SparseVector<Element>::value_type(bloc_start_idx + p, val));
						}
					}

					if(free_matrices)
					{
						if(vectors_to_free_D2[(D2.rowdim() - 1 - (this->pivot_rows_idxs_by_entry[i] - this->Npiv))/ NB_ROWS_PER_MULTILINE] == 1)
							D2[blc_row_idx][j][row_idx_in_blc].free ();
					}
				}

				if(free_matrices)
				{
					if(vectors_to_free_D2[(D2.rowdim() - 1 - (this->pivot_rows_idxs_by_entry[i] - this->Npiv))/ NB_ROWS_PER_MULTILINE] == 0)
						vectors_to_free_D2[(D2.rowdim() - 1 - (this->pivot_rows_idxs_by_entry[i] - this->Npiv))/ NB_ROWS_PER_MULTILINE] = 1;	//mark it for deletion
				}
			}
			else															//A & B
			{
				_tmp_Left.clear ();
				_tmp_Right.clear ();

				//for each bloc of A going horizontally
				uint32 blc_row_idx = (B2.rowdim() - 1
						- this->pivot_rows_idxs_by_entry[i]) / B2.bloc_height();

				uint32 row_idx_in_blc = ((B2.rowdim() - 1
						- this->pivot_rows_idxs_by_entry[i]) % B2.bloc_height())
						/ NB_ROWS_PER_MULTILINE;

				uint32 bloc_start_idx;
				uint16 val;
				uint32 idx;

				line = (B2.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i]) % NB_ROWS_PER_MULTILINE;

				//A
				_tmp_Left.push_back(typename SparseVector<Element>::value_type(pivot_columns_map[i], 1));
				//equivalent row_M->push_back (typename SparseMatrix<Element>::Row::value_type(i, 1));

				//B
				for(int j=0; j<(int)B2[blc_row_idx].size (); ++j)
				{
					if(B2[blc_row_idx][j].empty ())
						continue;

					if(B2[blc_row_idx][j][row_idx_in_blc].empty ())
						continue;

					bloc_start_idx = B2.FirstBlocsColumIndexes[blc_row_idx] + (B2.bloc_width() * j);

					if(B2[blc_row_idx][j][row_idx_in_blc].is_sparse ())
					{
						for (int p = 0; p<(int)B2[blc_row_idx][j][row_idx_in_blc].size(); ++p)
						{
							idx = B2[blc_row_idx][j][row_idx_in_blc].IndexData[p];
							val = B2[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);
							if(val != 0 )
								_tmp_Right.push_back(typename SparseVector<Element>::value_type(bloc_start_idx + idx, val));
						}
					}
					else
					{
						for (int p = 0; p < B2.bloc_width (); ++p)
						{
							val = B2[blc_row_idx][j][row_idx_in_blc].at_unchecked(line, p);
							if(val != 0)
								_tmp_Right.push_back(typename SparseVector<Element>::value_type(bloc_start_idx + p, val));
						}
					}

					if(free_matrices)
					{
						if(free_matrices)
						{
							if(vectors_to_free_B2[(B2.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] == 1)
								B2[blc_row_idx][j][row_idx_in_blc].free ();
						}
					}
				}

				if(free_matrices)
				{
					if(vectors_to_free_B2[(B2.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] == 0)
						vectors_to_free_B2[(B2.rowdim() - 1 - this->pivot_rows_idxs_by_entry[i])/ NB_ROWS_PER_MULTILINE] = 1;	//mark it for deletion
				}
			}

			push_rowsAB_to_rowM(_tmp_Left, _tmp_Right, *row_M);

			++piv;
		}

		if(free_matrices)
		{
			B2.free (true);
			D2.free (true);

			//free matrix blocs here
			free(vectors_to_free_B2);
			free(vectors_to_free_D2);
		}
	}


	///Param only_rows is used to remap the new pivots only, no columns are considered;
	///This is basically used when the echelon form is wanted, and not the RREF
	void combineInnerIndexer(ParallelIndexer<Element, Index>& inner_idxr, bool only_rows = false)
	{
		uint32 i, idx_in_B, idx_in_B2, rev_idx_in_A;

		//remap pivot rows indexes [lines of matrix A]
		uint32 next_piv = 0;
		for (i = 0; i < this->coldim; ++i) {
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			pivot_rows_idxs_by_entry[i] = next_piv++;
		}

		if(only_rows)
		{
			//pivots from matrix D
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
		}
	}















};


#endif /* INDEXER_H_ */
