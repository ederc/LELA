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

//	template <typename Element, typename Index>
//	void push_back_element(MultiLineVector<Element, Index>& v_piv, MultiLineVector<Element, Index>& v_non_piv,
//								const Index idx,
//								const Element val1,
//								const Element val2)
//	{
//		if(pivot_columns_map[idx] != MINUS_ONE)
//		{
//			v_piv.IndexData.push_back (pivot_columns_map[idx]);
//			v_piv.ValuesData.push_back (val1);
//			v_piv.ValuesData.push_back (val2);
//		}
//		else
//		{
//			v_non_piv.IndexData.push_back (non_pivot_columns_map[idx]);
//			v_non_piv.ValuesData.push_back (val1);
//			v_non_piv.ValuesData.push_back (val2);
//		}
//	}

#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

	template <typename Element, typename Index>
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

		std::vector<Index> idx_stack_elt_idx_in_line, idx_stack_bloc;
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

//	template <typename Element, typename Index>
//	void write_row_blocs_to_LeftToRight_matrix(SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
//												std::vector<SparseVector<Element> >& tmp_vect_A,
//												uint32 row_bloc_idx)
//	{
//		typename SparseVector<Element>::const_iterator it1, it2;
//
//		{
//			A[row_bloc_idx].clear ();
//			uint32 nb_blocs_per_A_dim;
//
//			if(A.isFilledWithEmptyBlocs())		//In case empty blocs should be appended
//			{
//				nb_blocs_per_A_dim = (uint32) std::ceil((float)A.coldim () / A.bloc_width ());
//				A.FirstBlocsColumIndexes[row_bloc_idx] = 0;
//			}
//			else
//			{
//				//search for the smallest/biggest column entry in the following B.bloc_height () rows
//				uint32 entry = (uint32)-1;
//				uint32 end = (uint32)0;
//
//				for(uint32 j=0; j < A.bloc_height (); ++j)
//				{
//					if(!(tmp_vect_A[j].empty ()))
//					{
//						entry = min(entry, tmp_vect_A[j].front ().first);
//						end = max(end, tmp_vect_A[j].back ().first);
//					}
//				}
//
//				nb_blocs_per_A_dim = (uint32) std::ceil((float)A.coldim () / A.bloc_width ());
//
//				if(entry != (uint32)-1)			//number of blocs in this row
//				{
//					nb_blocs_per_A_dim -= entry / A.bloc_width();
//					nb_blocs_per_A_dim -= (A.coldim() - end) / A.bloc_width();
//					A.FirstBlocsColumIndexes[row_bloc_idx] = entry / A.bloc_width() * A.bloc_width();
//				}
//				else
//				{
//					nb_blocs_per_A_dim = 0;
//					A.FirstBlocsColumIndexes[row_bloc_idx] = 0;
//				}
//			}
//
//			A[row_bloc_idx].resize (nb_blocs_per_A_dim);		//init memory of nb_blocs_per_A_dim blocs
//
//			for(uint32 k=0; k<nb_blocs_per_A_dim; ++k)			//initiliaze blocs in a row
//				A[row_bloc_idx][k].init(A.bloc_height (), A.bloc_width ());
//		}
//	}

	/**
	 * Constructs the submatrices A, B, C, D from matrix M
	 * In the case where the index maps are not yet constructed, the functions constructs the maps
	 * implicitly; in that case, the matrices are re_initiliazed with the new corresponding dimentions
	 */
	template<typename Element, typename Index>
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
					<< endl;

			//throw std::runtime_error ("Indexes not constructed yet! call Indexer.processMatrix() first");
			processMatrix(M);

			A = Matrix(Npiv, Npiv,
					Matrix::ArrangementDownTop_RightLeft,
					false, A.bloc_height(), A.bloc_width());
			A.acceptRowsHybrid = false;

			B = Matrix(Npiv, M.coldim() - Npiv,
					Matrix::ArrangementDownTop_LeftRight,
					true, B.bloc_height(), B.bloc_width()); //true: add fill blocs since they might be used
			B.acceptRowsHybrid = true;

			C = Matrix(M.rowdim() - Npiv, Npiv,
					Matrix::ArrangementDownTop_RightLeft,
					false, C.bloc_height(), C.bloc_width());
			C.acceptRowsHybrid = false;

			D = Matrix(M.rowdim() - Npiv, M.coldim() - Npiv,
					Matrix::ArrangementDownTop_LeftRight,
					true, D.bloc_height(), D.bloc_width()); //true: add fill blocs since they might be used
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

		SparseMatrix<Element> s_A (A.rowdim(), A.coldim()),
							s_B (B.rowdim(), B.coldim()),
							s_C (C.rowdim(), C.coldim()),
							s_D (D.rowdim(), D.coldim());

		//TODO: change that immediately
//		Indexer<uint32> one_line_indexer;
//		one_line_indexer.processMatrix(M);
//		one_line_indexer.constructSubMatrices(M, s_A, s_B, s_C, s_D, destruct_original_matrix);
//
//		MatrixUtils::copy(s_A, A);
//		MatrixUtils::copy(s_B, B);
//
//		MatrixUtils::copy(s_C, C);
//		MatrixUtils::copy(s_D, D);
//
//		return;

		SparseMatrix<uint16>::Row::const_iterator it;
		uint32 row_bloc_idx = 0;

		uint32 rows_idxs_horizontal_bloc[A.bloc_height ()];

		Index curr_vect_in_bloc = 0;

		//MATRIX A and B
		for (int i = M.coldim () - 1; i >= 0 ; --i)
		{
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)			//non pivot row
				continue;

			rows_idxs_horizontal_bloc[curr_vect_in_bloc] = pivot_rows_idxs_by_entry[i];
			curr_vect_in_bloc++;

			if(curr_vect_in_bloc == A.bloc_height ())
			{
				write_row_blocs_to_Left_Right_matrix(M, A, B, rows_idxs_horizontal_bloc, curr_vect_in_bloc, row_bloc_idx);
				//write_row_blocs_to_LeftToRight_matrix(B, tmp_vect_right, row_bloc_idx);

				row_bloc_idx++;
				curr_vect_in_bloc = 0;
			}

			if(destruct_original_matrix)
				M[pivot_rows_idxs_by_entry[i]].free ();

		}

		//write the left rows of A and B
		if(curr_vect_in_bloc != 0)
			write_row_blocs_to_Left_Right_matrix(M, A, B, rows_idxs_horizontal_bloc, curr_vect_in_bloc, row_bloc_idx);


		//MATRIX C and D
		curr_vect_in_bloc = 0;
		row_bloc_idx = 0;

		for (int i = M.rowdim() - 1; i >= 0; --i) 	//non pivot rows
		{
			if (non_pivot_rows_idxs[i] == MINUS_ONE)
				continue;

			rows_idxs_horizontal_bloc[curr_vect_in_bloc] = non_pivot_rows_idxs[i];
			curr_vect_in_bloc++;

			if(curr_vect_in_bloc == C.bloc_height ())
			{
				write_row_blocs_to_Left_Right_matrix(M, C, D, rows_idxs_horizontal_bloc, curr_vect_in_bloc, row_bloc_idx);

				row_bloc_idx++;
				curr_vect_in_bloc = 0;
			}

			if(destruct_original_matrix)
				M[non_pivot_rows_idxs[i]].free ();
		}

		//write the left rows of C and D
		if(curr_vect_in_bloc != 0)
			write_row_blocs_to_Left_Right_matrix(M, C, D, rows_idxs_horizontal_bloc, curr_vect_in_bloc, row_bloc_idx);
	}

};


#endif /* INDEXER_H_ */
