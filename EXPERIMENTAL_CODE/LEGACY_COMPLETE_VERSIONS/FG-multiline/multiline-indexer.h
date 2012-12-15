/*
 * multiline-index.h
 *
 *  Created on: 14 juin 2012
 *      Author: martani
 */

#ifndef MULTILINE_INDEX_H_
#define MULTILINE_INDEX_H_

//#include "../only-D/indexer.h"

using namespace LELA;

class MultiLineIndexer {
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
		assert(colSize >= 0);
		assert(rowSize >= 0);

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

	uint32 coldim, rowdim;

	//static const uint32 MINUS_ONE = 0xFFFFFFFF;	//stable for uint32!! lookout for other types
	static const uint32 MINUS_ONE = (uint32)-1 ;
	static const unsigned char MINUS_ONE_8 = (uint8)-1;// 0xFF;

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

	template <typename Element>
	long head(const MultiLineVector<Element>& v, const uint16 line_idx, Element& a, uint32& head_idx)
	{
		uint16 val=0;
		for(uint32 i=0; i<v.size(); ++i)
		{
			if((val = v.at_unchecked(line_idx, i)) != 0)
			{
				a = val;
				head_idx = i;
				return (long)v.IndexData[i];
			}
		}

		return -1;
	}

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
	template <class Element>
	void processMatrix(const SparseMatrix<Element>& M)
	{
		this->coldim = M.coldim ();
		this->rowdim = M.rowdim ();

		typename SparseMatrix<Element>::ConstRowIterator i_M;
		uint32 curr_row_idx = 0, entry;

		initArrays(this->rowdim, this->coldim);

		//for(uint32 i=0; i<this->coldim; ++i)
			//assert(pivot_rows_idxs_by_entry[i] == MINUS_ONE);

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

//				assert(piv_col_idx < this->Npiv);
				piv_col_idx++;
			}
			else
			{
				non_pivot_columns_map[i] = non_piv_col_idx;
				non_pivot_columns_rev_map[non_piv_col_idx] = i;

//				assert(non_piv_col_idx < (coldim - Npiv));
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

	//This is working only for the case where each multiline bloc has 2 lines
	template <class Element>
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
				h1 = head(*i_M, 0, h1_val, h1_idx);
				h2 = head(*i_M, 1, h2_val, h2_idx);

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

		assert(Npiv == piv_col_idx);
		assert(this->coldim - Npiv == non_piv_col_idx);

		for(uint32 i = 0; i < this->coldim; ++i){
			assert((pivot_columns_map[i] < Npiv) || (pivot_columns_map[i] == MINUS_ONE));
		}

		for(uint32 i = 0; i < this->coldim; ++i){
			assert((non_pivot_columns_map[i] < (this->coldim - Npiv))  || (non_pivot_columns_map[i] == MINUS_ONE));
		}

		_index_maps_constructed = true;
	}

	inline void push_back_element(MultiLineVector<uint16>& v_piv, MultiLineVector<uint16>& v_non_piv,
								const uint32 idx,
								const uint16 val1,
								const uint16 val2)
	{
		if(pivot_columns_map[idx] != MINUS_ONE)
		{
			v_piv.IndexData.push_back (pivot_columns_map[idx]);
			v_piv.ValuesData.push_back (val1);
			v_piv.ValuesData.push_back (val2);
		}
		else
		{
			v_non_piv.IndexData.push_back (non_pivot_columns_map[idx]);
			v_non_piv.ValuesData.push_back (val1);
			v_non_piv.ValuesData.push_back (val2);
		}
	}


	void constructSubMatrices(SparseMatrix<uint16>& M,
							  SparseMultilineMatrix<uint16>& A,
							  SparseMultilineMatrix<uint16>& B,
							  SparseMultilineMatrix<uint16>& C,
							  SparseMultilineMatrix<uint16>& D,
							  bool destruct_original)
	{
		//Rpiv <- the values in pivots_positions
		//Cpiv <- the indexes 0..pivots_size
		std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

		lela_check(M.coldim() == this->coldim);
		if(_index_maps_constructed == false)
		{
			report << "Indexes not constructed yet, call Indexer.processMatrix() first" << std::endl;
			throw "Index maps not constructed, call Indexer.processMatrix() first";
		}

		SparseMatrix<uint16>::Row::const_iterator it1, it2;
		uint32 curr_piv_AB = 0;

		uint32 row1=MINUS_ONE, row2=MINUS_ONE;

		//MATRIX A and B
		for (uint32 i = 0; i < M.coldim (); ++i) {
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)			//non pivot row
				continue;

			if(row1 == MINUS_ONE)
			{
				row1 = pivot_rows_idxs_by_entry[i];
				continue;
			}

			if(row2 != MINUS_ONE)
				throw "ERROR";
			else
				row2 = pivot_rows_idxs_by_entry[i];

			it1 = M[row1].begin ();
			it2 = M[row2].begin ();

			while(it1 != M[row1].end () && it2 != M[row2].end ())
			{
				if(it1->first < it2->first)
				{
					push_back_element(A[curr_piv_AB], B[curr_piv_AB], it1->first, it1->second, 0);
					++it1;
				}
				else if(it1->first > it2->first)
				{
					push_back_element(A[curr_piv_AB], B[curr_piv_AB], it2->first, 0, it2->second);
					++it2;
				}
				else	// it1->first == it2->first
				{
					push_back_element(A[curr_piv_AB], B[curr_piv_AB], it1->first, it1->second, it2->second);
					++it1;
					++it2;
				}
			}

			while(it1 != M[row1].end ())
			{
				push_back_element(A[curr_piv_AB], B[curr_piv_AB], it1->first, it1->second, 0);
				++it1;
			}

			while(it2 != M[row2].end ())
			{
				push_back_element(A[curr_piv_AB], B[curr_piv_AB], it2->first, 0, it2->second);
				++it2;
			}

			if(destruct_original)
			{
				M[row1].free ();
				M[row2].free ();
			}

			curr_piv_AB++;
			row1 = row2 = -1;
		}

		if(row1 != MINUS_ONE && row2 == MINUS_ONE){
			it1 = M[row1].begin ();

			while(it1 != M[row1].end ())
			{
				push_back_element(A[curr_piv_AB], B[curr_piv_AB], it1->first, it1->second, 0);
				++it1;
			}

			curr_piv_AB++;

			if(destruct_original)
			{
				M[row1].free ();
			}
		}

		//MATRIX C and D
		uint32 curr_piv_CD = 0;
		row1 = row2 = MINUS_ONE;
		for (uint32 i = 0; i < M.rowdim (); ++i) {
			if(non_pivot_rows_idxs[i] == MINUS_ONE)
				continue;

			if(row1==MINUS_ONE)
			{
				row1 = non_pivot_rows_idxs[i];
				continue;
			}

			if(row2 != MINUS_ONE)
				throw "ERROR";
			else
				row2 = non_pivot_rows_idxs[i];


			it1 = M[row1].begin ();
			it2 = M[row2].begin ();

			while(it1 != M[row1].end () && it2 != M[row2].end ())
			{
				if(it1->first < it2->first)
				{
					push_back_element(C[curr_piv_CD], D[curr_piv_CD], it1->first, it1->second, 0);
					++it1;
				}
				else if(it1->first > it2->first)
				{
					push_back_element(C[curr_piv_CD], D[curr_piv_CD], it2->first, 0, it2->second);
					++it2;
				}
				else	// it1->first == it2->first
				{
					push_back_element(C[curr_piv_CD], D[curr_piv_CD], it1->first, it1->second, it2->second);
					++it1;
					++it2;
				}
			}

			while(it1 != M[row1].end ())
			{
				push_back_element(C[curr_piv_CD], D[curr_piv_CD], it1->first, it1->second, 0);
				++it1;
			}

			while(it2 != M[row2].end ())
			{
				push_back_element(C[curr_piv_CD], D[curr_piv_CD], it2->first, 0, it2->second);
				++it2;
			}

			if(destruct_original)
			{
				M[row1].free ();
				M[row2].free ();
			}

			curr_piv_CD++;
			row1 = row2 = -1;
		}

		if(row1 != MINUS_ONE && row2 == MINUS_ONE){
			it1 = M[row1].begin ();

			while(it1 != M[row1].end ())
			{
				push_back_element(C[curr_piv_CD], D[curr_piv_CD], it1->first, it1->second, 0);
				++it1;
			}

			curr_piv_CD++;

			if(destruct_original)
			{
				M[row1].free ();
			}
		}

	}

	void constructSubMatrices(const SparseMultilineMatrix<uint16>& B,
							  const SparseMatrix<uint16>& D,
							  SparseMultilineMatrix<uint16>& B1,
							  SparseMultilineMatrix<uint16>& B2,
							  SparseMultilineMatrix<uint16>& D1,
							  SparseMultilineMatrix<uint16>& D2)
	{
		uint32 row1=MINUS_ONE, row2=MINUS_ONE;
		SparseMatrix<uint16>::Row::const_iterator it1, it2;
		uint32 curr_piv_D1D2=0;

		for (uint32 i = 0; i < D.coldim (); ++i) {
			if(pivot_rows_idxs_by_entry[i] == MINUS_ONE)			//non pivot row
				continue;

			if(row1 == MINUS_ONE)
			{
				row1 = pivot_rows_idxs_by_entry[i];
				continue;
			}

			if(row2 != MINUS_ONE)
				throw "ERROR";
			else
				row2 = pivot_rows_idxs_by_entry[i];

			it1 = D[row1].begin ();
			it2 = D[row2].begin ();

			while(it1 != D[row1].end () && it2 != D[row2].end ())
			{
				if(it1->first < it2->first)
				{
					push_back_element(D1[curr_piv_D1D2], D2[curr_piv_D1D2], it1->first, it1->second, 0);
					++it1;
				}
				else if(it1->first > it2->first)
				{
					push_back_element(D1[curr_piv_D1D2], D2[curr_piv_D1D2], it2->first, 0, it2->second);
					++it2;
				}
				else	// it1->first == it2->first
				{
					push_back_element(D1[curr_piv_D1D2], D2[curr_piv_D1D2], it1->first, it1->second, it2->second);
					++it1;
					++it2;
				}
			}

			while(it1 != D[row1].end ())
			{
				push_back_element(D1[curr_piv_D1D2], D2[curr_piv_D1D2], it1->first, it1->second, 0);
				++it1;
			}

			while(it2 != D[row2].end ())
			{
				push_back_element(D1[curr_piv_D1D2], D2[curr_piv_D1D2], it2->first, 0, it2->second);
				++it2;
			}

			curr_piv_D1D2++;
			row1 = row2 = -1;
		}

		if(row1 != MINUS_ONE && row2 == MINUS_ONE){
			it1 = D[row1].begin ();

			while(it1 != D[row1].end ())
			{
				push_back_element(D1[curr_piv_D1D2], D2[curr_piv_D1D2], it1->first, it1->second, 0);
				++it1;
			}

			curr_piv_D1D2++;
		}

		//B1 B2
		uint32 curr_piv_B1_B2 = 0;
		SparseMultilineMatrix<uint16>::ConstRowIterator i_B;
		uint32 idx;

		for(i_B = B.rowBegin (); i_B != B.rowEnd (); ++i_B)
		{
			for(uint32 j=0; j < i_B->size(); ++j)
			{
				idx = i_B->IndexData[j];

				if(pivot_columns_map[idx] != MINUS_ONE)
				{
					B1[curr_piv_B1_B2].push_back (pivot_columns_map[idx], i_B->at_unchecked(0, j), i_B->at_unchecked(1, j));
				}
				else
				{
					B2[curr_piv_B1_B2].push_back (non_pivot_columns_map[idx], i_B->at_unchecked(0, j), i_B->at_unchecked(1, j));
				}
			}

			++curr_piv_B1_B2;
		}
	}

	template <typename Multiline, typename SparseRow>
	inline void push_rowsAB_to_rowM(const Multiline& rowA, const Multiline& rowB, const uint16 line, SparseRow& rowM)
	{
		uint32 currA, currB, idxA, idxB;
		//rowA = &(A[this->pivot_rows_idxs_by_entry[i]/2]);
		//rowB = &(B[this->pivot_rows_idxs_by_entry[i]/2]);

		//nb_elts = rowA->size () + rowB->size ();
		//uint16 line = this->pivot_rows_idxs_by_entry[i] % 2;

		//rowM = &(M[piv]);
		//rowM->free ();
		//row_M->reserve (nb_elts);

		//it_A = rowA->begin ();
		//it_B = rowB->begin ();
		//while(it_A != rowA->end () && it_B != rowB->end ())

		currA = 0;
		currB = 0;

		while(currA < rowA.size () && currB < rowB.size ())
		{
			idxA = rowA.IndexData[currA];
			idxB = rowB.IndexData[currB];

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

	template <typename Matrix, typename Element>
	void reconstructMatrix(Matrix& M, const SparseMultilineMatrix<Element>& A, const SparseMultilineMatrix<Element>& B, const Matrix& D)
	{
		lela_check(A.rowdim() == B.rowdim());
		lela_check(A.rowdim() == A.coldim());
		lela_check(D.coldim() == B.coldim());

		typename SparseMultilineMatrix<Element>::ConstRow *rowA, *rowB;
		typename Matrix::ConstRow *rowD;
		typename Matrix::Row *row_M;

		typename Matrix::Row::const_iterator it_D;

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
				rowA = &(A[this->pivot_rows_idxs_by_entry[i]/2]);
				rowB = &(B[this->pivot_rows_idxs_by_entry[i]/2]);

				//nb_elts = rowA->size () + rowB->size ();
				uint16 line = this->pivot_rows_idxs_by_entry[i] % 2;

				row_M = &(M[piv]);
				row_M->free ();
				//row_M->reserve (nb_elts);

				push_rowsAB_to_rowM(*rowA, *rowB, line, *row_M);
			}

			++piv;
		}
	}


	template <typename Matrix, typename Element>
	void reconstructMatrix(Matrix& M, const SparseMultilineMatrix<Element>& A,
									  const SparseMultilineMatrix<Element>& B,
									  const SparseMultilineMatrix<Element>& D)
	{
		lela_check(A.rowdim() == B.rowdim());
		lela_check(A.rowdim() == A.coldim());
		lela_check(D.coldim() == B.coldim());

		typename SparseMultilineMatrix<Element>::ConstRow *rowA, *rowB, *rowD;
		typename Matrix::Row *row_M;

		typename Matrix::Row::const_iterator it_D;

		uint32 piv=0;
		uint16 line;

		for(uint32 i=0; i < this->coldim; i++)
		{
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			if(this->pivot_rows_idxs_by_entry[i] >= this->Npiv)		//D
			{
				rowD = &(D[(this->pivot_rows_idxs_by_entry[i] - this->Npiv)/2]);
				line = (this->pivot_rows_idxs_by_entry[i] - this->Npiv)%2;
				//nb_elts = rowD->size ();

				row_M = &(M[piv]);
				row_M->free ();
				//row_M->reserve (nb_elts);

				uint32 currD=0, idxD;
				while(currD < rowD->size ())
				{
					idxD = rowD->IndexData[currD];
					if(rowD->at_unchecked(line, currD) == 0)
					{
						++currD;
						continue;
					}

					row_M->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[idxD], rowD->at_unchecked(line, currD)));
					++currD;
				}
			}
			else															//A & B
			{
				rowA = &(A[this->pivot_rows_idxs_by_entry[i]/2]);
				rowB = &(B[this->pivot_rows_idxs_by_entry[i]/2]);

				//nb_elts = rowA->size () + rowB->size ();
				line = this->pivot_rows_idxs_by_entry[i] % 2;

				row_M = &(M[piv]);
				row_M->free ();
				//row_M->reserve (nb_elts);

				push_rowsAB_to_rowM(*rowA, *rowB, line, *row_M);
			}

			++piv;
		}
	}


	template <typename Matrix, typename Element>
	void reconstructMatrix(Matrix& M, const SparseMultilineMatrix<Element>& B, const SparseMultilineMatrix<Element>& D)
	{
		lela_check(D.coldim() == B.coldim());

		typename SparseMultilineMatrix<Element>::ConstRow *rowB, *rowD;
		typename Matrix::Row *row_M;

		typename Matrix::Row::const_iterator it_D;

		uint32 piv=0;
		uint16 line;

		for(uint32 i=0; i < this->coldim; i++)
		{
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			if(this->pivot_rows_idxs_by_entry[i] >= this->Npiv)		//D
			{
				rowD = &(D[(this->pivot_rows_idxs_by_entry[i] - this->Npiv)/2]);
				line = (this->pivot_rows_idxs_by_entry[i] - this->Npiv)%2;
				//nb_elts = rowD->size ();

				row_M = &(M[piv]);
				row_M->free ();
				//row_M->reserve (nb_elts);

				row_M->push_back (typename Matrix::Row::value_type(i, 1));

				uint32 currD=0, idxD;
				while(currD < rowD->size ())
				{
					idxD = rowD->IndexData[currD];
					if(rowD->at_unchecked(line, currD) == 0)
					{
						++currD;
						continue;
					}

					row_M->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[idxD], rowD->at_unchecked(line, currD)));
					++currD;
				}
			}
			else															//A & B
			{
				rowB = &(B[this->pivot_rows_idxs_by_entry[i]/2]);

				//nb_elts = rowA->size () + rowB->size ();
				line = this->pivot_rows_idxs_by_entry[i] % 2;

				row_M = &(M[piv]);
				row_M->free ();
				//row_M->reserve (nb_elts);

				row_M->push_back (typename Matrix::Row::value_type(i, 1));

				uint32 currB=0, idxB;
				while(currB < rowB->size ())
				{
					idxB = rowB->IndexData[currB];
					if(rowB->at_unchecked(line, currB) == 0)
					{
						++currB;
						continue;
					}

					row_M->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[idxB], rowB->at_unchecked(line, currB)));
					++currB;
				}
			}

			++piv;
		}
	}

	template <typename Matrix, typename Element>
	void reconstructMatrix(Matrix& M, const SparseMultilineMatrix<Element>& B)
	{
		typename SparseMultilineMatrix<Element>::ConstRow *rowB;
		typename Matrix::Row *row_M1, *row_M2;

		uint32 piv=0, j;
		//uint16 line;

		for(uint32 i=0; i < this->coldim; i++)
		{
			if(this->pivot_rows_idxs_by_entry[i] == MINUS_ONE)
				continue;

			rowB = &(B[this->pivot_rows_idxs_by_entry[i]/2]);

			j=i+1;
			while(j < this->coldim && this->pivot_rows_idxs_by_entry[j] == MINUS_ONE)
				++j;

			//LAST ONE
			if(j == this->coldim)
			{
				row_M1 = &(M[piv]);
				row_M1->free ();
				row_M1->push_back (typename Matrix::Row::value_type(i, 1));

				uint32 currB=0, idxB;
				while(currB < rowB->size ())
				{
					idxB = rowB->IndexData[currB];

					if(rowB->at_unchecked(0, currB) != 0)
					{
						row_M1->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[idxB], rowB->at_unchecked(0, currB)));
					}

					currB++;
				}

				break;
			}

			//assert(this->pivot_rows_idxs_by_entry[i]+1 == pivot_rows_idxs_by_entry[j]);

			//nb_elts = rowA->size () + rowB->size ();
			//line = this->pivot_rows_idxs_by_entry[i] % 2;

			row_M1 = &(M[piv]);
			row_M1->free ();
			row_M1->reserve (rowB->size());

			row_M2 = &(M[piv+1]);
			row_M2->free ();
			row_M2->reserve (rowB->size());

			row_M1->push_back (typename Matrix::Row::value_type(i, 1));
			row_M2->push_back (typename Matrix::Row::value_type(j, 1));

			uint32 currB=0, idxB;
			while(currB < rowB->size ())
			{
				idxB = rowB->IndexData[currB];

				if(rowB->at_unchecked(0, currB) != 0)
				{
					row_M1->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[idxB], rowB->at_unchecked(0, currB)));
				}

				if(rowB->at_unchecked(1, currB) != 0)
				{
					row_M2->push_back (typename Matrix::Row::value_type(non_pivot_columns_rev_map[idxB], rowB->at_unchecked(1, currB)));
				}

				++currB;
			}

			piv += 2;
			i = j;
		}
	}
















	///param only_rows is used to remap the new pivots only, no columns are considered;
	///this is basically used when the echelon form is wanted, and not the RREF
	void combineInnerIndexer(MultiLineIndexer& inner_idxr, bool only_rows = false)
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

	template <typename Elt>
	bool compare(Indexer<Elt>& idxr)
	{
		std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
		bool pass = true;

		if(this->coldim != idxr.coldim)
		{
			report << "DIFFER at coldim " << this->coldim << "    " << idxr.coldim << endl;
			return false;
		}
		if(this->rowdim != idxr.rowdim)
		{
			report << "DIFFER at rowdim " << this->rowdim << "    " << idxr.rowdim << endl;
			return false;
		}

		for(uint32 i=0; i<this->coldim; ++i)
		{
			if(pivot_rows_idxs_by_entry[i] != idxr.pivot_rows_idxs_by_entry[i])
			{
				report << "pivot_rows_idxs_by_entry DIFFER at index " << i << endl;
				pass = false;
			}

			if(pivot_columns_map[i] != idxr.pivot_columns_map[i])
			{
				report << "pivot_columns_map DIFFER at index " << i << endl;
				pass = false;
			}
			if(non_pivot_columns_map[i] != idxr.non_pivot_columns_map[i])
			{
				report << "non_pivot_columns_map DIFFER at index " << i << endl;
				pass = false;
			}

			if(pivot_columns_rev_map[i] != idxr.pivot_columns_rev_map[i])
			{
				report << "pivot_columns_map DIFFER at index " << i << endl;
				pass = false;
			}

			if(non_pivot_columns_rev_map[i] != idxr.non_pivot_columns_rev_map[i])
			{
				report << "non_pivot_columns_rev_map DIFFER at index " << i << endl;
				pass = false;
			}

			if(!pass)
				break;
		}

		for(uint32 i=0; i<this->rowdim; ++i)
		{
			if(non_pivot_rows_idxs[i] != idxr.non_pivot_rows_idxs[i])
			{
				report << "non_pivot_rows_idxs DIFFER at index " << i << endl;
				pass = false;

				if(!pass)
					break;
			}
		}

		return pass;
	}
};

#endif /* MULTILINE_INDEX_H_ */
