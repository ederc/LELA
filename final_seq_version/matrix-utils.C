/*
 * matrix-util.C
 *
 *  Created on: 4 juil. 2012
 *      Author: martani
 */


#ifndef MATRIX_UTILS_C_
#define MATRIX_UTILS_C_

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "matrix-utils.h"

using namespace LELA;
using namespace std;

template <typename Element, typename Index>
void MatrixUtils::copy(const SparseMatrix<Element>& A, SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B)
{
	lela_check(A.coldim () == B.coldim ());
	lela_check(A.rowdim () == B.rowdim ());

	if(A.rowdim () == 0)
	{
		for(uint32 j=0; j<B.rowBlocDim (); ++j)
			B[j].clear ();
		return;
	}

	switch(B.blocArrangement)
	{
		case SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::ArrangementTopDown_LeftRight:
			sparse_to_bloc_copyTopDowLeftRight(A, B);
		break;

		case SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::ArrangementTopDown_RightLeft:
			throw std::logic_error ("NotImplemented");
		break;

		case SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::ArrangementDownTop_LeftRight:
			sparse_to_bloc_copyDownTopLeftRight(A, B);
		break;
		case SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::ArrangementDownTop_RightLeft:
			sparse_to_bloc_copyDownTopRightLeft(A, B);
		break;
	}

	if(B.acceptRowsHybrid)
		transformSparseToHybridRows(B);

}


template <typename Element, typename Index>
void MatrixUtils::sparse_to_bloc_copyTopDowLeftRight(const SparseMatrix<Element>& A, SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B)
{
	//typename SparseMatrix<Element>::ConstRowIterator i_A, i_A__plus_1;
	typename SparseMatrix<Element>::Row::const_iterator it1, it2;

	uint32 row_bloc_idx, bloc_idx_in_row1, line_idx_in_bloc, elt_idx_in_line1;
	uint32 bloc_idx_in_row2, elt_idx_in_line2;
	uint32 i;

	for(i=0; i < A.rowdim(); ++i)
	{
		row_bloc_idx = i / B.bloc_height ();
		line_idx_in_bloc = i % B.bloc_height () / 2;

		//cout << "Line "<< i << " row_bloc_idx " << row_bloc_idx << " line_idx_in_bloc " << line_idx_in_bloc << endl;

		if(i % B.bloc_height () == 0)			//init blocs
		{
			B[row_bloc_idx].clear ();

			//cout << "after clear\n";

			uint32 nb_blocs_per_A_dim;

			if(B.isFilledWithEmptyBlocs())		//In case empty blocs should be appended
			{
				//cout << "B.isFilledWithEmptyBlocs()" << endl;
				nb_blocs_per_A_dim = (uint32) std::ceil((float)A.coldim () / B.bloc_width ());
				B.FirstBlocsColumIndexes[row_bloc_idx] = 0;
			}
			else
			{
				//search for the smallest column entry in the following B.bloc_height () rows
				uint32 entry = (uint32)-1;

				for(uint16 j=0; j<B.bloc_height () && j < A.rowdim () - i; ++j)
				{
					if(!(A[i+j].empty ()))
						entry = min(entry, A[i+j][0].first);
				}

				//cout << "entry is " << entry << endl;

				if(entry != (uint32)-1)			//number of blocs in this row
				{
					nb_blocs_per_A_dim = (uint32) std::ceil((float)A.coldim () / B.bloc_width ());
					nb_blocs_per_A_dim -= entry / B.bloc_width();

					//cout << "nb blocs here " << nb_blocs_per_A_dim << endl;

					B.FirstBlocsColumIndexes[row_bloc_idx] = (entry / B.bloc_width()) * B.bloc_width();
				}
				else
				{
					nb_blocs_per_A_dim = 0;
					B.FirstBlocsColumIndexes[row_bloc_idx] = 0;
				}
			}

			B[row_bloc_idx].resize (nb_blocs_per_A_dim);		//init memory of nb_blocs_per_A_dim blocs

			//TODO: initialize only blocs where there is actual data
			for(uint32 k=0; k<nb_blocs_per_A_dim; ++k)			//initiliaze blocs in a row
				B[row_bloc_idx][k].init(B.bloc_height(), B.bloc_width());
		}

		//cout << "cool\n";
		//if(i_A->empty ())
			//continue;

		it1 = A[i].begin (); //i_A->begin ();

		uint i_plus_1;
		if(i < A.rowdim() - 1)
		{
			i_plus_1 = i+1;
			it2 = A[i_plus_1].begin ();
		}
		else
		{
			i_plus_1 = i;
			it2 = A[i_plus_1].end ();
		}


		//cout << "ok \n";

		while(it1 !=  A[i].end () && it2 !=  A[i_plus_1].end ())
		{
			bloc_idx_in_row1 = (it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line1 = (it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			bloc_idx_in_row2 = (it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line2 = (it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			if(it1->first < it2->first)
			{
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line1);
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (it1->second);
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (0);

				++it1;
			}
			else if(it1->first > it2->first)
			{
				B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line2);
				B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].ValuesData.push_back (0);
				B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].ValuesData.push_back (it2->second);

				++it2;
			}
			else	// it1->first == it2->first
			{
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line1);
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (it1->second);
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (it2->second);
				++it1;
				++it2;
			}
		}

		//cout << "while 1\n";

		while(it1 !=  A[i].end ())
		{
			bloc_idx_in_row1 = (it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line1 = (it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			//cout << "bloc_idx_in_row1 " << bloc_idx_in_row1 << "\telt_idx_in_line1 " << elt_idx_in_line1 << endl;

			B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line1);
			B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (it1->second);
			B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (0);
			++it1;
		}

		//cout << "while 2\n";
		while(it2 !=  A[i_plus_1].end ())
		{
			bloc_idx_in_row2 = (it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line2 = (it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line2);
			B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].ValuesData.push_back (0);
			B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].ValuesData.push_back (it2->second);
			++it2;
		}

		++i;
		//cout << "while 3\n";
	}
}

template <typename Element, typename SparseBlocMatrix_>
void MatrixUtils::sparse_to_bloc_copyDownTopRightLeft(const SparseMatrix<Element>& A, SparseBlocMatrix_& B)
{
	typename SparseMatrix<Element>::Row::const_iterator it1, it2;

	uint32 row_bloc_idx, bloc_idx_in_row1, line_idx_in_bloc, elt_idx_in_line1;
	uint32 bloc_idx_in_row2, elt_idx_in_line2;
	uint32 i;

	const uint32 last_line_idx = A.rowdim ()-1;
	uint32 i_offset;

	std::vector<uint32> idx_stack_elt, idx_stack_bloc;
	std::vector<Element> v1_stack, v2_stack;

	i = A.rowdim();

	while(i != 0)
	{
		--i;

		i_offset = last_line_idx - i;

		row_bloc_idx = i_offset / B.bloc_height ();
		line_idx_in_bloc = i_offset % B.bloc_height ()/2;

		//cout << "Line "<< i << " row_bloc_idx " << row_bloc_idx << " line_idx_in_bloc " << line_idx_in_bloc << endl;

		if(i_offset % B.bloc_height () == 0)			//init blocs
		{
			B[row_bloc_idx].clear ();

			//cout << "after clear\n";

			uint32 nb_blocs_per_A_dim;

			if(B.isFilledWithEmptyBlocs())		//In case empty blocs should be appended
			{
				//cout << "B.isFilledWithEmptyBlocs()" << endl;
				nb_blocs_per_A_dim = (uint32) std::ceil((float)A.coldim () / B.bloc_width ());
				B.FirstBlocsColumIndexes[row_bloc_idx] = 0;
			}
			else
			{
				//search for the smallest/biggest column entry in the following B.bloc_height () rows
				uint32 entry = (uint32)-1;
				uint32 end = (uint32)0;

				int nb=0;
				for(int j=i; j>=0 && j > (int)(i - B.bloc_height ()); --j)
				{
					if(!(A[j].empty ()))
					{
						entry = min(entry, A[j].front ().first);
						end = max(end, A[j].back ().first);
					}
					nb++;
				}

				//cout << "entry is " << entry << " end is " << end << " nb " << nb <<endl;

				if(entry != (uint32)-1)			//number of blocs in this row
				{
					nb_blocs_per_A_dim = (uint32) std::ceil((float)(A.coldim () - entry) / B.bloc_width());
					nb_blocs_per_A_dim -= (A.coldim () - 1 - end) / B.bloc_width();

					B.FirstBlocsColumIndexes[row_bloc_idx] = (A.coldim () - 1 - end) / B.bloc_width() * B.bloc_width();
				}
				else
				{
					nb_blocs_per_A_dim = 0;
					B.FirstBlocsColumIndexes[row_bloc_idx] = 0;
				}
			}

			//cout << "nb blocs " << nb_blocs_per_A_dim << "\t B.FirstBlocsColumIndexes[row_bloc_idx] " << B.FirstBlocsColumIndexes[row_bloc_idx] << endl;

			B[row_bloc_idx].resize (nb_blocs_per_A_dim);		//init memory of nb_blocs_per_A_dim blocs

			for(uint32 k=0; k<nb_blocs_per_A_dim; ++k)			//initiliaze blocs in a row
			{
				B[row_bloc_idx][k].init(B.bloc_height(), B.bloc_width());
				//B[row_bloc_idx][k] = typename SparseBlocMatrix_::BlocType (B.bloc_height(), B.bloc_width());
			}
		}

		it1 = A[i].begin (); //i_A->begin ();

		uint i_plus_1;
		if(i > 0)
		{
			i_plus_1 = i-1;
			it2 = A[i_plus_1].begin ();
		}
		else
		{
			//cout << "LAST LINE \n";
			i_plus_1 = i;
			it2 = A[i_plus_1].end ();
		}


		//cout <<"ok\n";
		idx_stack_elt.clear();
		idx_stack_bloc.clear();
		v1_stack.clear();
		v2_stack.clear();

		while(it1 !=  A[i].end () && it2 !=  A[i_plus_1].end ())
		{
			bloc_idx_in_row1 = (A.coldim () - 1 - it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line1 = (A.coldim () - 1 - it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			bloc_idx_in_row2 = (A.coldim () - 1 - it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line2 = (A.coldim () - 1 - it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			if(it1->first < it2->first)
			{
				idx_stack_elt.push_back (elt_idx_in_line1);
				idx_stack_bloc.push_back(bloc_idx_in_row1);
				v1_stack.push_back (it1->second);
				v2_stack.push_back (0);

				++it1;
			}
			else if(it1->first > it2->first)
			{
				idx_stack_elt.push_back (elt_idx_in_line2);
				idx_stack_bloc.push_back(bloc_idx_in_row2);
				v1_stack.push_back (0);
				v2_stack.push_back (it2->second);

				++it2;
			}
			else	// it1->first == it2->first
			{
				idx_stack_elt.push_back (elt_idx_in_line1);
				idx_stack_bloc.push_back(bloc_idx_in_row1);
				v1_stack.push_back (it1->second);
				v2_stack.push_back (it2->second);
				++it1;
				++it2;
			}
		}

		//cout << "while 1\n";

		while(it1 !=  A[i].end ())
		{
			bloc_idx_in_row1 = (A.coldim () - 1 - it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line1 = (A.coldim () - 1 - it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			//cout << "bloc_idx_in_row1 " << bloc_idx_in_row1 << "\telt_idx_in_line1 " << elt_idx_in_line1 << endl;

			idx_stack_elt.push_back (elt_idx_in_line1);
			idx_stack_bloc.push_back(bloc_idx_in_row1);
			v1_stack.push_back (it1->second);
			v2_stack.push_back (0);
			++it1;
		}

		//cout << "while 2\n";
		while(it2 !=  A[i_plus_1].end ())
		{
			bloc_idx_in_row2 = (A.coldim () - 1 - it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line2 = (A.coldim () - 1 - it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			idx_stack_elt.push_back (elt_idx_in_line2);
			idx_stack_bloc.push_back(bloc_idx_in_row2);
			v1_stack.push_back (0);
			v2_stack.push_back (it2->second);
			++it2;
		}
		//cout << "done \n";

		i = i_plus_1;

		for(int k=idx_stack_elt.size()-1; k>=0; --k)
		{
			B[row_bloc_idx][idx_stack_bloc[k]][line_idx_in_bloc].IndexData.push_back (idx_stack_elt[k]);
			B[row_bloc_idx][idx_stack_bloc[k]][line_idx_in_bloc].ValuesData.push_back (v1_stack[k]);
			B[row_bloc_idx][idx_stack_bloc[k]][line_idx_in_bloc].ValuesData.push_back (v2_stack[k]);
		}

	}
}


template <typename Element, typename Index>
void MatrixUtils::sparse_to_bloc_copyDownTopLeftRight(const SparseMatrix<Element>& A, SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B)
{
	//typename SparseMatrix<Element>::ConstRowIterator i_A;
	typename SparseMatrix<Element>::Row::const_iterator it1, it2;

	uint32 row_bloc_idx, bloc_idx_in_row1, line_idx_in_bloc, elt_idx_in_line1;
	uint32 bloc_idx_in_row2, elt_idx_in_line2;
	uint32 i;

	const uint32 last_line_idx = A.rowdim ()-1;
	uint32 i_offset;

	i = A.rowdim();
	//i_A = A.rowEnd ();

	while(i != 0)	//TODO: caution! overflow
	{
		--i;
		//--i_A;

		i_offset = last_line_idx - i;

		row_bloc_idx = i_offset / B.bloc_height ();
		line_idx_in_bloc = i_offset % B.bloc_height ()/2;

		//cout << "Line "<< i << " row_bloc_idx " << row_bloc_idx << " line_idx_in_bloc " << line_idx_in_bloc << endl;

		if(i_offset % B.bloc_height () == 0)			//init blocs
		{
			B[row_bloc_idx].clear ();

			//cout << "after clear\n";

			uint32 nb_blocs_per_A_dim;

			if(B.isFilledWithEmptyBlocs())		//In case empty blocs should be appended
			{
				//cout << "B.isFilledWithEmptyBlocs()" << endl;
				nb_blocs_per_A_dim = (uint32) std::ceil((float)A.coldim () / B.bloc_width ());
				B.FirstBlocsColumIndexes[row_bloc_idx] = 0;
			}
			else
			{
				//search for the smallest/biggest column entry in the following B.bloc_height () rows
				uint32 entry = (uint32)-1;
				uint32 end = (uint32)0;

				for(int j=i; j>=0 && j > (int)(i - B.bloc_height ()); --j)
				{
					if(!(A[j].empty ()))
					{
						entry = min(entry, A[j][0].first);
						end = max(end, A[j].back ().first);
					}
				}

				//cout << "entry is " << entry << " end is " << end << endl;

				nb_blocs_per_A_dim = (uint32) std::ceil((float)A.coldim () / B.bloc_width ());

				if(entry != (uint32)-1)			//number of blocs in this row
				{
					nb_blocs_per_A_dim -= entry / B.bloc_width();
					nb_blocs_per_A_dim -= (A.coldim() - end) / B.bloc_width();
					B.FirstBlocsColumIndexes[row_bloc_idx] = entry / B.bloc_width() * B.bloc_width();
				}
				else
				{
					nb_blocs_per_A_dim = 0;
					B.FirstBlocsColumIndexes[row_bloc_idx] = 0;
				}
			}

			//cout << "nb blocs " << nb_blocs_per_A_dim << endl;

			B[row_bloc_idx].resize (nb_blocs_per_A_dim);		//init memory of nb_blocs_per_A_dim blocs

			//TODO: initialize only blocs where there is actual data
			for(uint32 k=0; k<nb_blocs_per_A_dim; ++k)			//initiliaze blocs in a row
				B[row_bloc_idx][k].init(B.bloc_height(), B.bloc_width());
		}

		it1 = A[i].begin (); //i_A->begin ();

		uint i_plus_1;
		if(i > 0)
		{
			i_plus_1 = i-1;
			it2 = A[i_plus_1].begin ();
		}
		else
		{
			//cout << "LAST LINE \n";
			i_plus_1 = i;
			it2 = A[i_plus_1].end ();
		}


		//cout << "ok \n";

		while(it1 !=  A[i].end () && it2 !=  A[i_plus_1].end ())
		{
			bloc_idx_in_row1 = (it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line1 = (it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			bloc_idx_in_row2 = (it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line2 = (it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			if(it1->first < it2->first)
			{
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line1);
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (it1->second);
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (0);

				++it1;
			}
			else if(it1->first > it2->first)
			{
				B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line2);
				B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].ValuesData.push_back (0);
				B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].ValuesData.push_back (it2->second);

				++it2;
			}
			else	// it1->first == it2->first
			{
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line1);
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (it1->second);
				B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (it2->second);
				++it1;
				++it2;
			}
		}

		//cout << "while 1\n";

		while(it1 !=  A[i].end ())
		{
			bloc_idx_in_row1 = (it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line1 = (it1->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			//cout << "bloc_idx_in_row1 " << bloc_idx_in_row1 << "\telt_idx_in_line1 " << elt_idx_in_line1 << endl;

			B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line1);
			B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (it1->second);
			B[row_bloc_idx][bloc_idx_in_row1][line_idx_in_bloc].ValuesData.push_back (0);
			++it1;
		}

		//cout << "while 2\n";
		while(it2 !=  A[i_plus_1].end ())
		{
			bloc_idx_in_row2 = (it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line2 = (it2->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].IndexData.push_back (elt_idx_in_line2);
			B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].ValuesData.push_back (0);
			B[row_bloc_idx][bloc_idx_in_row2][line_idx_in_bloc].ValuesData.push_back (it2->second);
			++it2;
		}
		//cout << "done \n";

		i = i_plus_1;
	}

}


template <typename Element, typename Index>
void MatrixUtils::copy(
		SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
		SparseMatrix<Element>& B, bool destruct_original)
{
	check_equal_or_raise_exception(A.coldim (), B.coldim ());
	check_equal_or_raise_exception(A.rowdim (), B.rowdim ());

	for(uint32 j=0; j<B.rowdim (); ++j)
		B[j].clear ();

	switch(A.blocArrangement)
	{
		case SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::ArrangementTopDown_LeftRight:
		{
			transformHybridToSparseRows(A);
			bloc_to_sparse_copyTopDowLeftRight(A, B);
		}
		break;

		case SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::ArrangementTopDown_RightLeft:
			throw std::logic_error ("NotImplemented");
		break;

		case SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::ArrangementDownTop_LeftRight:
			bloc_to_sparse_copyDownTopLeftRight(A, B, destruct_original);
		break;
		case SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::ArrangementDownTop_RightLeft:
		{
			transformHybridToSparseRows(A);
			bloc_to_sparse_copyDownTopRightLeft(A, B);
		}
		break;
	}
}

template <typename Element, typename Index>
void MatrixUtils::bloc_to_sparse_copyTopDowLeftRight(const SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A, SparseMatrix<Element>& B)
{
	typename SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::Row::const_iterator it_bloc;
	//typename SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::BlocType::Row::const_iterator it;

	uint32 curr_row_B = 0;

	for(uint32 i=0; i<A.rowBlocDim (); ++i)			//for each row of blocs in A
	{
		if(A[i].empty ())
		{
			curr_row_B += A.bloc_height();
			continue;
		}

		for(uint32 j=0;j<A[i].size (); ++j)		//for each bloc in the row
		{
			if(A[i][j].empty ())
				continue;

			uint32 bloc_idx = A.FirstBlocsColumIndexes[i] + A.bloc_width() * j;
			uint32 idx;
			Element val1, val2;

			for(uint16 k=0; k<A[i][j].bloc_height () && k*2<(B.rowdim() - curr_row_B); ++k)	//for each row in the bloc
			{
				if(A[i][j][k].empty ())
					continue;

				for(uint32 p=0; p<A[i][j][k].size (); ++p)
				//for(it = A[i][j][k].begin (); it != A[i][j][k].end (); ++it)			//for each element in this row
				{
					idx = A[i][j][k].IndexData[p];
					val1 = A[i][j][k].at_unchecked(0, p);
					val2 = A[i][j][k].at_unchecked(1, p);

					if(val1 != 0)
						B[curr_row_B+k*2].push_back(typename SparseMatrix<Element>::Row::value_type(bloc_idx + idx, val1));

					if(val2 != 0)
						B[curr_row_B+k*2+1].push_back(typename SparseMatrix<Element>::Row::value_type(bloc_idx + idx, val2));
				}
			}
		}

		curr_row_B += A.bloc_height();
	}
}


template <typename Element, typename Index>
void MatrixUtils::bloc_to_sparse_copyDownTopLeftRight(
		SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
		SparseMatrix<Element>& B,
		bool destruct_original)
{
	//typename SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::Row::const_iterator it_bloc;
	//typename SparseBlocMatrix<SparseMultilineBloc<Element, Index> >::BlocType::Row::const_iterator it;

	uint32 curr_row_B = B.rowdim() - 1;

	for(uint32 i=0; i<A.rowBlocDim (); ++i)			//for each row of blocs in A, starting at the last row in B
	{
		if(A[i].empty ())
		{
			curr_row_B -= A.bloc_height();
			continue;
		}

		if(A.isFilledWithEmptyBlocs())
				check_equal_or_raise_exception(A.FirstBlocsColumIndexes[i], 0);

		for(uint32 j=0;j<A[i].size (); ++j)		//for each bloc in the row
		{
			if(A[i][j].empty ())
				continue;

			uint32 bloc_idx = A.FirstBlocsColumIndexes[i] + (A.bloc_width() * j);

			uint32 idx;
			Element val1, val2;
			for(uint16 k=0; k<A[i][j].bloc_height () && k*2 <= curr_row_B; ++k)	//for each row in the bloc
			{
				if(A[i][j][k].empty ())
					continue;

				B[curr_row_B-k*2].reserve (A[i][j][k].size ());

				if(curr_row_B >= (uint32)k*2+1)
					B[curr_row_B-k*2-1].reserve (A[i][j][k].size ());

				if(A[i][j][k].is_sparse ())
				{
					for(uint32 p=0; p<A[i][j][k].size (); ++p)
					{
						idx = A[i][j][k].IndexData[p];
						val1 = A[i][j][k].at_unchecked(0, p);
						val2 = A[i][j][k].at_unchecked(1, p);

						if(val1 != 0)
							B[curr_row_B-k*2].push_back(typename SparseMatrix<Element>::Row::value_type(bloc_idx + idx, val1));

						if(val2 != 0)
							B[curr_row_B-k*2-1].push_back(typename SparseMatrix<Element>::Row::value_type(bloc_idx + idx, val2));
					}
				}
				else
				{
					for(uint32 p=0; p<A[i][j][k].size (); ++p)
					{
						val1 = A[i][j][k].at_unchecked(0, p);
						val2 = A[i][j][k].at_unchecked(1, p);

						if(val1 != 0)
							B[curr_row_B-k*2].push_back(typename SparseMatrix<Element>::Row::value_type(bloc_idx + p, val1));

						if(val2 != 0)
							B[curr_row_B-k*2-1].push_back(typename SparseMatrix<Element>::Row::value_type(bloc_idx + p, val2));
					}
				}

				if(destruct_original)
					A[i][j][k].free();
			}

			if(destruct_original)
					A[i][j].free();

		}

		curr_row_B -= A.bloc_height();
		//break;
	}

	if(destruct_original)
		A.free ();
}


template <typename Element, typename Index>
void MatrixUtils::bloc_to_sparse_copyDownTopRightLeft(const SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A, SparseMatrix<Element>& B)
{
	//cout << "IN COPY bloc size " << A.bloc_width() << ", " << A.bloc_height() << endl;

	uint32 curr_row_B = B.rowdim() - 1;

	for(uint32 i=0; i<A.rowBlocDim (); ++i)			//for each row of blocs in A, starting at the last row in B
	{
		if(A[i].empty ())
		{
			curr_row_B -= A.bloc_height();
			continue;
		}

		for(int j=A[i].size ()-1; j>=0; --j)		//for each bloc in the row starting from the end
		{
			if(A[i][j].empty ())
				continue;

			uint32 bloc_idx = A.coldim() - 1 - A.FirstBlocsColumIndexes[i] - A.bloc_width() * j;		//the colum index of the bloc starting from *THE LEFT*

			uint32 idx;
			Element val1, val2;
			for(uint16 k=0; k<A[i][j].bloc_height () && k*2 <= curr_row_B; ++k)	//for each row in the bloc
			{
				if(A[i][j][k].empty ())
					continue;

				for(int p=A[i][j][k].size ()-1; p>=0; --p)
				//for(it = A[i][j][k].begin (); it != A[i][j][k].end (); ++it)			//for each element in this row
				{
					idx = A[i][j][k].IndexData[p];
					val1 = A[i][j][k].at_unchecked(0, p);
					val2 = A[i][j][k].at_unchecked(1, p);

					if(val1 != 0)
						B[curr_row_B-k*2].push_back(typename SparseMatrix<Element>::Row::value_type(bloc_idx - idx, val1));

					if(val2 != 0)
						B[curr_row_B-k*2-1].push_back(typename SparseMatrix<Element>::Row::value_type(bloc_idx - idx, val2));
				}
			}
		}

		curr_row_B -= A.bloc_height();
	}
}


template <typename Element, typename Index>
void MatrixUtils::transformSparseToHybridRows(SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A)
{
	MultiLineVector<Element, Index> tmp;

	for(uint32 i=0; i<A.rowBlocDim(); ++i)
	{
		for(uint32 j=0; j<A[i].size (); ++j)
		{
			for(uint16 k=0; k<A[i][j].bloc_height (); ++k)
			{
				if(A[i][j][k].empty ())
					continue;

				if((float)A[i][j][k].size () / (float)A.bloc_width() < A.get_HYBRID_REPRESENTATION_THRESHOLD ())
					continue;

				tmp.clear();
				uint32 idx=0;

				for(uint32 p=0; p<A.bloc_width(); ++p) // p<A[i][j][k].size (); ++p)
				{
					if(idx < A[i][j][k].size () && A[i][j][k].IndexData[idx] == p)
					{
						tmp.ValuesData.push_back(A[i][j][k].at_unchecked(0, idx));
						tmp.ValuesData.push_back(A[i][j][k].at_unchecked(1, idx));
						++idx;
					}
					else
					{
						tmp.ValuesData.push_back(0);
						tmp.ValuesData.push_back(0);
					}
				}

				A[i][j][k].swap(tmp);
			}
		}
	}
}

template <typename Element, typename Index>
void MatrixUtils::transformHybridToSparseRows(SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A)
{
	for(uint32 i=0; i<A.rowBlocDim(); ++i)
	{
		for(uint32 j=0; j<A[i].size (); ++j)
		{
			for(uint16 k=0; k<A[i][j].bloc_height (); ++k)
			{
				if(A[i][j][k].is_sparse ())
					continue;

				MultiLineVector<Element, Index> tmp;
				Element e1, e2;

				for(uint32 p=0; p<A.bloc_width(); ++p) // p<A[i][j][k].size (); ++p)
				{
					e1 = A[i][j][k].at_unchecked(0, p);
					e2 = A[i][j][k].at_unchecked(1, p);

					if(e1 != 0 || e2 != 0)
					{
						tmp.ValuesData.push_back(e1);
						tmp.ValuesData.push_back(e2);
						tmp.IndexData.push_back(p);
					}
				}

				A[i][j][k].swap(tmp);
			}
		}
	}
}


template <typename Element, typename Index>
void MatrixUtils::transformHybridMatrixToSparseMatrix(SparseMultilineMatrix<Element, Index>& A)
{
	MultiLineVector<Element, Index> tmp;
	Element e1, e2;

	for(uint32 i=0; i<A.multiline_rowdim(); ++i)
	{
		if(A[i].is_sparse ())
			continue;

		tmp.clear ();
		
		if(!A[i].is_sparse(A.coldim ()))
		{
			for(uint32 p=0; p<A[i].size (); ++p) // p<A[i][j][k].size (); ++p)
			{
				e1 = A[i].at_unchecked(0, p);
				e2 = A[i].at_unchecked(1, p);

				if(e1 != 0 || e2 != 0)
				{
					tmp.ValuesData.push_back(e1);
					tmp.ValuesData.push_back(e2);
					tmp.IndexData.push_back(p);
				}
			}

			A[i].swap(tmp);
		}
	}
}


template <typename Element, typename Index>
bool MatrixUtils::equal(const SparseMultilineMatrix<Element, Index>& A, const SparseMultilineMatrix<Element, Index>& B)
{
	if(A.multiline_rowdim () != B.multiline_rowdim ())
		return false;

	for(uint32 i=0; i<A.multiline_rowdim(); ++i)
	{
		if(!A[i].equal(B[i], A.coldim()))
		{
			std::cout << "NOT EQUAL LINE " << i << std::endl;
			return false;
		}
	}

	return true;
}

template <typename Element, typename Index>
void MatrixUtils::copyBlocToMultilineMatrix_DownTop_RightLeft(const SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& inMatrix,
		SparseMultilineMatrix<Element>& outMatrix)
{
	outMatrix = SparseMultilineMatrix<Element> (inMatrix.rowdim (), inMatrix.coldim ());

	int curr_row_base = inMatrix.rowdim() / NB_ROWS_PER_MULTILINE + inMatrix.rowdim() % NB_ROWS_PER_MULTILINE - 1;

	cout << "MatrixUtils::copyBlocToMultilineMatrix_DownTop_RightLeft" << endl;

	for(uint32 i=0; i<inMatrix.rowBlocDim (); ++i)			//for each row of blocs in A, starting at the last row in B
	{
		if(inMatrix[i].empty ())
		{
			curr_row_base -= inMatrix.bloc_height () / NB_ROWS_PER_MULTILINE;
			continue;
		}

		for(int j=inMatrix[i].size ()-1; j>=0; --j) //for each bloc in the row
		{
			if(inMatrix[i][j].empty ())
				continue;

			//check_equal_or_raise_exception(inMatrix.FirstBlocsColumIndexes[i], 0);
			//uint32 bloc_idx = inMatrix.FirstBlocsColumIndexes[i] + (inMatrix.bloc_width () * j);
			uint32 bloc_idx = inMatrix.coldim() - 1 - inMatrix.FirstBlocsColumIndexes[i] - inMatrix.bloc_width() * j;

			//cout << "bloc_idx " << bloc_idx << endl;

			uint32 idx;
			Element val1, val2;
			for (uint16 k = 0;
					k < inMatrix[i][j].bloc_height(); ++k) //for each row in the bloc
			{
				if(inMatrix[i][j][k].empty ())
					continue;

				if(inMatrix[i][j][k].is_sparse ())
				{
					for (int p = inMatrix[i][j][k].size() - 1; p >= 0; --p)
					{
						idx = inMatrix[i][j][k].IndexData[p];
						val1 = inMatrix[i][j][k].at_unchecked(0, p);
						val2 = inMatrix[i][j][k].at_unchecked(1, p);

						outMatrix[curr_row_base-k].IndexData.push_back(bloc_idx - idx);
						outMatrix[curr_row_base-k].ValuesData.push_back(val2);
						outMatrix[curr_row_base-k].ValuesData.push_back(val1);
					}
				}
				else
				{
					for (int p = inMatrix.bloc_width () - 1; p >= 0; --p)
					{
						val1 = inMatrix[i][j][k].at_unchecked(0, p);
						val2 = inMatrix[i][j][k].at_unchecked(1, p);

						if(val1 != 0 || val2 != 0)
						{
							outMatrix[curr_row_base-k].IndexData.push_back(bloc_idx - p);
							outMatrix[curr_row_base-k].ValuesData.push_back(val2);
							outMatrix[curr_row_base-k].ValuesData.push_back(val1);
						}
					}
				}
			}
		}

		curr_row_base -= inMatrix.bloc_height () / NB_ROWS_PER_MULTILINE;
	}
}











template <typename Element, typename Index>
void MatrixUtils::dumpMatrixAsPbmImage(SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A, const char *outputFileName)
{
	SparseMatrix<Element> tmp (A.rowdim(), A.coldim());
	copy(A, tmp);
	dumpMatrixAsPbmImage(tmp, outputFileName);
}


template <typename Element>
void MatrixUtils::dumpMatrixAsPbmImage(const SparseMatrix<Element>& A, const char *outputFileName)
{
	typename SparseMatrix<Element>::ConstRowIterator i_A = A.rowBegin ();
	uint32 j, col, m, row_size;
	char buffer[512];
	unsigned char output_byte = 0;
	m = A.coldim ();

	FILE *outStream = fopen(outputFileName, "wb");

	//magic PBM header
#ifdef __LP64__	//64 bit machine
	sprintf(buffer, "P4\n# matrix size(%lu, %lu)\n%lu %lu\n", A.rowdim (), A.coldim (), A.coldim (), A.rowdim ());
#else			//32 bit machine
	sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", A.rowdim (), A.coldim (), A.coldim (), A.rowdim ());
#endif

	fwrite(buffer, sizeof(char), strlen(buffer), outStream);


	while(i_A != A.rowEnd())
	{
		row_size = i_A->size ();

		j=0;
		for(col=0; col<m; ++col){

			if(j<row_size && (*i_A)[j].first == col)
			{
				output_byte |= (1 << (7 - (col%8)));
				j++;
			}
			else
			{
				output_byte &= ~(1 << (7 - (col%8)));
			}

			if(col%8 == 7) //flush byte every 8 cols
			{
				fwrite(&output_byte, sizeof(unsigned char), 1, outStream);
				output_byte = 0;
			}
		}

		if(col%8 != 0)
			fwrite(&output_byte, sizeof(unsigned char), 1, outStream);

		fflush(outStream);

		++i_A;
	}

	fclose(outStream);
}


template<typename Element, typename Index>
void MatrixUtils::dumpMatrixAsPbmImage(const SparseMultilineMatrix<Element, Index>& A, const char *outputFileName)
{
	typename SparseMultilineMatrix<Element>::ConstRowIterator i_A = A.rowBegin();
	uint32 j, col, m, row_size;
	char buffer[512];
	unsigned char output_byte = 0;
	m = A.coldim();

	FILE *outStream = fopen(outputFileName, "wb");

	//magic PBM header
#ifdef __LP64__	//64 bit machine
	sprintf(buffer, "P4\n# matrix size(%lu, %lu)\n%lu %lu\n", A.rowdim (), A.coldim (), A.coldim (), A.rowdim ());
#else			//32 bit machine
	sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", A.rowdim(),
			A.coldim(), A.coldim(), A.rowdim());
#endif

	fwrite(buffer, sizeof(char), strlen(buffer), outStream);

	while (i_A != A.rowEnd()) //for each multiline
	{
		row_size = i_A->size();

		for (uint16 i = 0; i < A.nb_lines_per_bloc(); ++i) //for each line in the multiline
		{
			j = 0;
			for (col = 0; col < m; ++col)
			{
				if (j < row_size && i_A->IndexData[j] == col)
				{
					if (i_A->at(i, j) != 0)
						output_byte |= (1 << (7 - (col % 8)));
					else
						output_byte &= ~(1 << (7 - (col % 8)));

					++j;
				}
				else
				{
					output_byte &= ~(1 << (7 - (col % 8)));
				}

				if (col % 8 == 7) //flush byte every 8 cols
				{
					fwrite(&output_byte, sizeof(unsigned char), 1, outStream);
					output_byte = 0;
				}
			}

			if (col % 8 != 0)
				fwrite(&output_byte, sizeof(unsigned char), 1, outStream);

			fflush(outStream);
		}

		++i_A;
	}

	fclose(outStream);
}

template <typename Element, typename Index>
bool MatrixUtils::equal(const Modular<Element>& R, SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A, const SparseMatrix<Element>& B)
{
	SparseMatrix<Element> tmp (A.rowdim(), A.coldim());
	copy(A, tmp);

	Context<Modular<Element> > ctx (R);

	return BLAS3::equal(ctx, B, tmp);
}


template <typename Element, typename Index>
bool MatrixUtils::equal(const Modular<Element>& R, const SparseMatrix<Element>& A, SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& B)
{
	return equal(R, B, A);
}



void MatrixUtils::show_mem_usage(std::string msg)
{
	std::string unit = "KB"; // KB, MB
	double vm, rss;
	process_mem_usage(vm, rss);
	if(vm > 1024)
	{
		vm = vm / 1024.0;
		rss = rss / 1024.0;
		unit = "MB";
	}

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "[[[" << msg << "]]]\t\t" << " Memory (RSS: " << rss << unit << "; VM: " << vm << unit << ")" << std::endl;
}

uint32 MatrixUtils::loadF4Modulus(const char *fileName)
{
	uint32 mod;

    std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	FILE *f = fopen(fileName, "r");
	if (f == NULL)
	{
		report << "Can't open " << fileName << std::endl;
		throw std::runtime_error ("Can't open file");
	}

	fseek(f, 2 * sizeof(uint32), SEEK_SET);
	if(fread(&mod, sizeof(uint32),     1,f) != 1)
	{
		report << "Error while reading file " << fileName << std::endl;
		return 0;
	}

	assert(mod >= 2);

	fclose(f);
	return mod;
}


///Reads the matrix row by row from the file, does not load the whole file to memory. more efficient than dump_matrix.c
///Caller must free memory once the matrix is not longer needed
template<class Ring>
void MatrixUtils::loadF4Matrix(const Ring &R, SparseMatrix<typename Ring::Element>& A,
		const char *fileName)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	// Code adapted from C version of dump_matrix.c

	FILE *f = fopen(fileName, "r");
	if (f == NULL)
	{
		commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "Can't open " << fileName << std::endl;
		throw std::runtime_error ("Can't open file");
	}

	unsigned short int *nz, *onz;
	unsigned int *pos, *opos;
	unsigned int *sz, *osz;
	unsigned int n;
	unsigned int m;
	unsigned int mod;
	unsigned long long nb;

	if (fread(&n, sizeof(unsigned int), 1, f) != 1)
		throw std::runtime_error ("Error while reading file");
	if (fread(&m, sizeof(unsigned int), 1, f) != 1)
		throw std::runtime_error ("Error while reading file");
	if (fread(&mod, sizeof(unsigned int), 1, f) != 1)
		throw std::runtime_error ("Error while reading file");
	if (fread(&nb, sizeof(unsigned long long), 1, f) != 1)
		throw std::runtime_error ("Error while reading file");

	report << "loading file " << fileName << std::endl;
	report << n << " x " << m << " matrix" << std::endl;
	report << "mod " << mod << std::endl;
	{
		double Nz = (double) (n) * (double) (m);
		Nz = (double) (nb) / Nz;
		Nz *= 100.0;
		report << "Nb of Nz elements " << nb << " (density " << Nz << ")"
				<< std::endl;
	}

	A = SparseMatrix<typename Ring::Element> (n, m);
	typename SparseMatrix<typename Ring::Element>::RowIterator i_A;

	uint32 i;
	onz = nz = (unsigned short int*) malloc(nb * sizeof(unsigned short int));
	opos = pos = (unsigned int*) malloc(nb * sizeof(unsigned int));
	osz = sz = (unsigned int*) malloc(n * sizeof(unsigned int));

	if (fread(nz, sizeof(short int), nb, f) != nb)
		throw std::runtime_error ("Error while reading file");

	if (fread(pos, sizeof(unsigned int), nb, f) != nb)
		throw std::runtime_error ("Error while reading file");

	if (fread(sz, sizeof(unsigned int), n, f) != n)
		throw std::runtime_error ("Error while reading file");

	for (i_A = A.rowBegin(), i = 0; i < n; i++, ++i_A)
	{
		const unsigned int szi = sz[i];
		unsigned int j;

		i_A->reserve(szi);

		for (j = 0; j < szi; j++)
		{
			i_A->push_back(
					typename Vector<Ring>::Sparse::value_type(pos[j],
							typename Ring::Element()));
			R.init(i_A->back().second, nz[j]);

		}

		nz += szi;
		pos += szi;
	}

	//free memory
	free(onz);
	free(opos);
	free(osz);
	fclose(f);
}




template<class Ring>
void MatrixUtils::loadF4Matrix__low_memory(const Ring &R,
		SparseMatrix<typename Ring::Element>& A, const char *fileName)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL,
			INTERNAL_DESCRIPTION);
	// Code adapted from C version of dump_matrix.c

	FILE *f = fopen(fileName, "r");
	if (f == NULL)
	{
		commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "Can't open " << fileName << std::endl;
		throw std::runtime_error ("Can't open file");
	}

	uint16 *nz;
	uint32 *pos;
	uint32 sz;
	uint32 n;
	uint32 m;
	uint32 mod;
	uint64 nb;

	if (fread(&n, sizeof(unsigned int), 1, f) != 1)
		throw std::runtime_error ("Error while reading file");
	if (fread(&m, sizeof(unsigned int), 1, f) != 1)
		throw std::runtime_error ("Error while reading file");
	if (fread(&mod, sizeof(unsigned int), 1, f) != 1)
		throw std::runtime_error ("Error while reading file");
	if (fread(&nb, sizeof(unsigned long long), 1, f) != 1)
		throw std::runtime_error ("Error while reading file");

	report << n << " x " << m << " matrix" << std::endl;
	report << "mod " << mod << std::endl;
	{
		double Nz = (double) (n) * (double) (m);
		Nz = (double) (nb) / Nz;
		Nz *= 100.0;
		report << "Nb of Nz elements " << nb << " (density " << Nz << ")"
				<< std::endl;
	}

	A = SparseMatrix<typename Ring::Element>(n, m);
	typename SparseMatrix<typename Ring::Element>::RowIterator i_A;


	uint32 i;
	nz = new unsigned short int[m]; //has a size of at most a full row of the matrix
	pos = new uint32[m];

	//save the arrays original pointers
	uint16 *oNz = nz;
	uint32 *oPos = pos;

	uint32 header_size = sizeof(uint32) * 3 + sizeof(uint64); //size of n, m, mod and nb in the header of the file
	uint64 row_sizes_offset, row_values_offset, row_positions_offset; //cursors in the file

	//row sizes if positioned after the values and the positions of the elements in the file
	row_sizes_offset = nb * sizeof(uint16) + nb * sizeof(uint32) + header_size;
	row_values_offset = header_size;
	row_positions_offset = nb * sizeof(uint16) + header_size;

	for (i_A = A.rowBegin(), i = 0; i < n; i++, ++i_A)
	{
		//get the size of the current row
		fseek(f, row_sizes_offset, SEEK_SET);
		if (fread(&sz, sizeof(uint32), 1, f) != 1)
			throw "Error while reading file";

		row_sizes_offset += sizeof(uint32);

		//assert(sz <= m);
		//number of elements in a row at max equal to the size of a row in the matrix

		//read sz elements from the values part of the file
		fseek(f, row_values_offset, SEEK_SET);
		if (fread(nz, sizeof(uint16), sz, f) != sz)
			throw "Error while reading file";

		row_values_offset += sz * sizeof(uint16);

		//read sz elements from the posistions part of the file
		fseek(f, row_positions_offset, SEEK_SET);
		if (fread(pos, sizeof(uint32), sz, f) != sz)
			throw "Error while reading file";

		row_positions_offset += sz * sizeof(uint32);

		i_A->reserve(sz);
		for (uint32 j = 0; j < sz; j++)
		{
			i_A->push_back(
					typename Vector<Ring>::Sparse::value_type(pos[j],
							typename Ring::Element()));
			R.init(i_A->back().second, nz[j]);
			//assert(pos[j] < m);
		}
	}

	//free memory
	delete[] oNz;
	delete[] oPos;
	fclose(f);
}



template<class Ring>
void MatrixUtils::loadF4Matrix__low_memory_syscall_no_checks(const Ring &R,
		SparseMatrix<typename Ring::Element>& A, const char *fileName)
{
	std::ostream &report = commentator.report(Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << std::endl << "Loading file: " << fileName << std::endl;

	// Code adapted from C version of dump_matrix.c
	struct stat fileStat;
	int f = open(fileName, O_RDONLY);

	if (f < 0)
	{
		commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "Can't open " << fileName << std::endl;
		throw std::runtime_error ("Can't open file");
		return;
	}

	uint16 *nz;
	uint32 *pos;
	uint32 sz;
	uint32 n;
	uint32 m;
	uint32 mod;
	uint64 nb;

	if (read(f, &n, sizeof(unsigned int)) != sizeof(unsigned int))
		throw std::runtime_error ("Error while reading file");
	if (read(f, &m, sizeof(unsigned int)) != sizeof(unsigned int))
		throw std::runtime_error ("Error while reading file");
	if (read(f, &mod, sizeof(unsigned int)) != sizeof(unsigned int))
		throw std::runtime_error ("Error while reading file");
	if (read(f, &nb, sizeof(unsigned long long)) != sizeof(unsigned long long))
		throw std::runtime_error ("Error while reading file");

	stat(fileName, &fileStat);

	report << n << " x " << m << " matrix" << std::endl;
	report << "mod " << mod << std::endl;
	{
		double Nz = (double) (n) * (double) (m);
		Nz = (double) (nb) / Nz;
		Nz *= 100.0;
		report << "Nb of Nz elements " << nb << " (density " << Nz << "% ) -\tsize "
			   << fileStat.st_size / 1024 / 1024 << " MB" << std::endl;
	}

	A = SparseMatrix<typename Ring::Element>(n, m);
	typename SparseMatrix<typename Ring::Element>::RowIterator i_A;


	uint32 i;
	nz = new unsigned short int[m]; //has a size of at most a full row of the matrix
	pos = new uint32[m];

	//save the arrays original pointers
	uint16 *oNz = nz;
	uint32 *oPos = pos;

	uint32 header_size = sizeof(uint32) * 3 + sizeof(uint64); //size of n, m, mod and nb in the header of the file
	uint64 row_sizes_offset, row_values_offset, row_positions_offset; //cursors in the file

	//row sizes if positioned after the values and the positions of the elements in the file
	row_sizes_offset = nb * sizeof(uint16) + nb * sizeof(uint32) + header_size;
	row_values_offset = header_size;
	row_positions_offset = nb * sizeof(uint16) + header_size;

	int ret;

	for (i_A = A.rowBegin(), i = 0; i < n; i++, ++i_A)
	{
		//get the size of the current row
		lseek(f, row_sizes_offset, SEEK_SET);
		ret = read(f, &sz, sizeof(uint32));
//		if (ret != sizeof(uint32))
//			throw "Error while reading file";

		row_sizes_offset += sizeof(uint32);

		//assert(sz <= m);
		//number of elements in a row at max equal to the size of a row in the matrix

		//read sz elements from the values part of the file
		lseek(f, row_values_offset, SEEK_SET);
		ret = read(f, nz, sizeof(uint16) * sz);
//		if (ret != sizeof(uint16) * sz)
//			throw "Error while reading file";

		row_values_offset += sz * sizeof(uint16);

		//read sz elements from the posistions part of the file
		lseek(f, row_positions_offset, SEEK_SET);
		ret = read(f, pos, sizeof(uint32) * sz);
//		if (ret != sizeof(uint32) * sz)
//			throw "Error while reading file";

		row_positions_offset += sz * sizeof(uint32);

		i_A->reserve(sz);
		for (uint32 j = 0; j < sz; j++)
		{
			i_A->push_back(
					typename Vector<Ring>::Sparse::value_type(pos[j], nz[j]));

			//TODO: If using Mosular<>.init, it could take 50 times slower to load a matrix from file
//							typename Ring::Element()));
//			R.init(i_A->back().second, nz[j]);
			//assert(pos[j] < m);
		}
	}

	ret++;

	//free memory
	delete[] oNz;
	delete[] oPos;
	close(f);
}


template <typename Element>
std::pair<uint64, double> MatrixUtils::getMatrixSizeAndDensity(const SparseMatrix<Element>& A, bool exact)
{
	uint64 nb_elts = 0;

    typename SparseMatrix<Element>::ConstRowIterator i_A = A.rowBegin ();

	while(i_A != A.rowEnd ())
	{
			nb_elts += i_A->size ();
			++i_A;
	}

	double Nz = (double)(A.rowdim ())*(double)(A.coldim ());
	Nz = (double)(nb_elts)/Nz;
	Nz *= 100.0;

	return std::pair<uint64, double>(nb_elts, Nz);
}

template <typename Element>
std::pair<uint64, double> MatrixUtils::getMatrixSizeAndDensity(const SparseMultilineMatrix<Element>& A, bool exact)
{
	uint64 nb_elts = 0;

    for(uint32 i=0; i<A.multiline_rowdim(); ++i)
	{
			nb_elts += A[i].size ();
	}

	double Nz = (double)(A.rowdim ())*(double)(A.coldim ());
	Nz = (double)(nb_elts)/Nz;
	Nz *= 100.0;

	return std::pair<uint64, double>(nb_elts, Nz);
}

//Returns an **approximation** of the density
template <typename Element, typename Index>
std::pair<uint64, double> MatrixUtils::getMatrixSizeAndDensity(const SparseBlocMatrix<SparseMultilineBloc<Element, Index> >& A,
		bool exact)
{
	uint64 nb_elts = 0;

	if(!exact)
		for(uint32 i=0; i<A.rowBlocDim(); ++i)
		{
			for(uint32 j=0; j<A[i].size (); ++j)
			{
				for(uint16 k=0; k<A[i][j].bloc_height (); ++k)
				{
					nb_elts += A[i][j][k].size ();
				}
			}
		}
	else
		for(uint32 i=0; i<A.rowBlocDim(); ++i)
		{
			for(uint32 j=0; j<A[i].size (); ++j)
			{
				for(uint16 k=0; k<A[i][j].bloc_height (); ++k)
				{
					Element e1, e2;

					for(uint32 p=0; p<A[i][j][k].size (); ++p) // p<A[i][j][k].size (); ++p)
					{
						e1 = A[i][j][k].at_unchecked(0, p);
						e2 = A[i][j][k].at_unchecked(1, p);

						if(e1 != 0)
							nb_elts++;

						if(e2 != 0)
							nb_elts++;
					}
				}
			}
		}

	double Nz = (double)(A.rowdim ())*(double)(A.coldim ());
	Nz = (double)(nb_elts)/Nz;
	Nz *= 100.0;

	return std::pair<uint64, double>(nb_elts, Nz);
}


//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void MatrixUtils::process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
			   >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
			   >> utime >> stime >> cutime >> cstime >> priority >> nice
			   >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

#endif	/* MATRIX_UTILS_C_ */
