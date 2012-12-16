 
 
/*
 * test-FG-multiline.C
 *
 *  Created on: 14 juin 2012
 *      Author: martani
 */

#include <iostream>

#include "FG-types.h"
#include "../../new-faugere-lachartre/matrix-util.h"
#include "../../new-faugere-lachartre/matrix-op.h"
#include "../../new-faugere-lachartre/indexer.h"
#include "multiline-indexer.h"

#include "../../../util/support.h"

using namespace LELA;
using namespace std;

const int NB_LINES_PER_BLOC = 8;

#define START_UNROLL_CODE				\
	while(x<xl)					\
	{						\
		idx = v.IndexData[x];			\
		val1 = v.at_unchecked(0, x);		\
		val2 = v.at_unchecked(1, x);

#define MIDDLE_UNROLL_CODE				\
	++x;						\
	}						\
							\
	for(x=xl; x<sz; x+=STEP)			\
	{

#define MIDDLE_UNROLL_CODE2				\
		for(uint8 t=0; t<STEP; ++t)		\
		{					\
			idx = v.IndexData[x+t];		\
							\
			val1 = v.at_unchecked(0, x+t);	\
			val2 = v.at_unchecked(1, x+t);

#define END_UNROLL_CODE					\
		}					\
	}

template <typename Element>
void dumpMatrixAsPbmImage(const SparseMultilineMatrix<Element>& A, const char *outputFileName)
{
	typename SparseMultilineMatrix<Element>::ConstRowIterator i_A = A.rowBegin ();
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

	while(i_A != A.rowEnd())					//for each multiline
	{
		row_size = i_A->size ();

		for(uint16 i=0; i < A.nb_lines_per_bloc (); ++i)		//for each line in the multiline
		{
			j=0;
			for(col=0; col<m; ++col){
				if(j<row_size && i_A->IndexData[j] == col)
				{
					if(i_A->at (i, j) != 0)
						output_byte |= (1 << (7 - (col%8)));
					else
						output_byte &= ~(1 << (7 - (col%8)));

					++j;
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
		}

		++i_A;
	}

	fclose(outStream);
}

template <typename Matrix, typename Element>
bool equal(Modular<Element>& R, SparseMultilineMatrix<Element>& A, Matrix& B)
{
	SparseVector<Element> v;
	typename SparseMultilineMatrix<Element>::ConstRowIterator i_A = A.rowBegin ();
	typename Matrix::ConstRowIterator i_B = B.rowBegin ();

	Context<Modular<Element> > ctx (R);

	if(A.rowdim() != B.rowdim ())
		return false;

	uint32 line=0;

	while(i_A != A.rowEnd ())			//for each multiline of A
	{
		for(uint16 i=0; i<A.nb_lines_per_bloc (); ++i)		//for each line in the multiline
		{
			//cout << "Line " << t << " size: " << i_A->size () << endl << "[";
			for(uint32 j=0; j < i_A->size (); ++j)
			{
				if(i_A->at(i, j) != 0)
				{
					v.push_back(typename Matrix::Row::value_type(i_A->IndexData[j], i_A->at(i, j)));
				}
			}

			if(v.empty () && i_B == B.rowEnd ())
			{
				continue;
			}

			++line;
			if(!BLAS1::equal(ctx, *i_B, v))
			{
				BLAS1::write(ctx, cout, v);
				cout << endl << endl;
				BLAS1::write(ctx, cout, *i_B);
				cout << endl << endl;

				cout << "Diff on line: " << line;
				cout << endl << endl;


				for(uint32 j=0; j < i_B->size (); ++j)
				{
					if((*i_B)[j].first != v[j].first)
					{
						cout << "DIFF ON Index: " << (*i_B)[j].first << "      " << v[j].first<< endl;
						return false;
					}
					else
					if((*i_B)[j].second != v[j].second)
					{
						cout << "DIFF ON Value: " << (*i_B)[j].second << "      " << v[j].second
								<< " At index: " << (*i_B)[j].first << endl;
						return false;
					}

				}


				return false;
			}
			else
				;//cout << "Line " << line << " OK " << endl;

			v.clear ();
			++i_B;
		}

		++i_A;
	}

	return true;
}


template <typename Matrix, typename Element>
bool equal_reverse(Modular<Element>& R, SparseMultilineMatrix<Element>& A, Matrix& B)
{
	SparseVector<Element> v;
	typename SparseMultilineMatrix<Element>::ConstRowIterator i_A = A.rowEnd ();
	typename Matrix::ConstRowIterator i_B = B.rowEnd ();

	Context<Modular<Element> > ctx (R);

	if(A.rowdim() != B.rowdim ())
		return false;

	uint32 line=A.rowdim() - 1;
	int last_null_rows = A.rowdim() % A.nb_lines_per_bloc();
	if(last_null_rows != 0)
		last_null_rows = A.nb_lines_per_bloc() - last_null_rows;

	--i_A;
	--i_B;
	int nb_elts=0;

	do
	{
		for(int i=A.nb_lines_per_bloc ()-1; i>=0; --i)		//for each line in the multiline
		{
			nb_elts=0;
			//cout << "Line " << t << " size: " << i_A->size () << endl << "[";
			for(uint32 j=0; j < i_A->size (); ++j)
			{
				if(i_A->at(i, j) != 0)
				{
					v.push_back(typename Matrix::Row::value_type(i_A->IndexData[j], i_A->at(i, j)));
					//++nb_elts;
				}
			}

			//if(nb_elts == 0)
				//cout << "ADDED " << nb_elts << endl;

			if(last_null_rows>0)		//skip null rows at the end
			{
				last_null_rows--;
				v.clear ();
				continue;
			}

			//if(v.empty () && i_B == B.rowEnd ())
				//return true;

			if(!BLAS1::equal(ctx, *i_B, v))
			{
				cout << "MULTI_LINE" << endl;
				BLAS1::write(ctx, cout, v);
				cout << endl << endl;
				cout << "SPARSE VECTOR" << endl;
				BLAS1::write(ctx, cout, *i_B);
				cout << endl << endl;

				cout << "Diff on line: " << line;
				cout << endl << endl;
				return false;
			}
			else
				;//cout << "Line " << line << " OK " << endl;

			--line;
			v.clear ();
			--i_B;
		}

		--i_A;
	} while(i_A != A.rowBegin ());			//for each multiline of A

	return true;
}

inline void razArray64(uint64 *arrays[], uint16 nb_arrays, uint32 arrSize)
{
	for(uint16 j=0; j<nb_arrays; ++j)
		memset(arrays[j], 0, arrSize*sizeof(uint64));
}

//Copy the lines of a multiline bloc into arrays
void copyMultilineToDenseArrays64(const MultiLineVector<uint16>& v, uint64 *arrays[], uint16 nb_arrays)
{
	register uint32 idx;
	for(uint32 i=0; i<v.size (); ++i)
	{
		idx = v.IndexData[i];

		for(uint16 j=0; j<nb_arrays; ++j)
			arrays[j][idx] = (uint64)v.at_unchecked(j, i);
	}
}

TIMER_DECLARE_(PREPARE);
TIMER_DECLARE_(COMPUTE);

//Copy a list of dense arrays to a multiline vector
template <typename Ring>
void copyDenseArraysToMultilineVector64(const Ring& R, uint64 *arrays[], uint16 nb_arrays,
		uint32 size,
		MultiLineVector<uint16>& v,
		bool reduce = true)
{
	MultiLineVector<uint16> tmp (v.nb_lines());
	typename Ring::Element e;

	bool index_added = false;
	uint16 nb_pending_zeros = 0;

	if(reduce)
		for (uint32 i = 0; i < size; ++i){

			index_added = false;
			nb_pending_zeros = 0;

			for(uint16 j=0; j<nb_arrays; ++j)
				if((arrays[j][i] % R._modulus != 0))
				{
					if(!index_added)
					{
						tmp.IndexData.push_back(i);

						for(uint16 l=0; l<nb_pending_zeros; ++l)
							tmp.ValuesData.push_back(0);

						index_added = true;
					}

					ModularTraits<typename Ring::Element>::reduce (e, arrays[j][i], R._modulus);
					tmp.ValuesData.push_back(e);
				}
				else
					if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
						nb_pending_zeros++;
					else
						tmp.ValuesData.push_back(0);
		}
	else
		for (uint32 i = 0; i < size; ++i){
			index_added = false;
			nb_pending_zeros = 0;

			for(uint16 j=0; j<nb_arrays; ++j)
				if((arrays[j][i] != 0))
				{
					assert(arrays[j][i] % R._modulus != 0);
					if(!index_added)
					{
						tmp.IndexData.push_back(i);

						for(uint16 l=0; l<nb_pending_zeros; ++l)
							tmp.ValuesData.push_back(0);

						index_added = true;
					}

					//ModularTraits<typename Ring::Element>::reduce (e, arrays[j][i], R._modulus);
					tmp.ValuesData.push_back(arrays[j][i]);
				}
				else
					if(!index_added)	//if element is zero in this column; and no non zero is found yet, count that it should be appended
						nb_pending_zeros++;
					else
						tmp.ValuesData.push_back(0);
		}

	v.swap(tmp);
}

inline void axpy2(uint16 *Av[],
		const uint16 nb_av,
		const MultiLineVector<uint16>& v,
		uint64 *arrays[])
{

	TIMER_START_(PREPARE);
	uint32 non_zero_Av_idx[NB_LINES_PER_BLOC][NB_LINES_PER_BLOC],
		   nb_non_zero_Av[NB_LINES_PER_BLOC],
		   nb_non_zeo_lines=0;

	uint64 non_zero_lines_dx[NB_LINES_PER_BLOC];

	bool line_added = false;

	for(uint16 i=0; i<NB_LINES_PER_BLOC; ++i)
		nb_non_zero_Av[i]=0;

	for(uint16 i=0; i<NB_LINES_PER_BLOC; ++i)
	{
		line_added = false;
		for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
		{
			if(Av[i][j] != 0)
			{
				non_zero_Av_idx[i][nb_non_zero_Av[i]] = j;
				nb_non_zero_Av[i]++;

				if(!line_added)
				{
					non_zero_lines_dx[nb_non_zeo_lines] = i;
					nb_non_zeo_lines++;
					line_added=true;
				}
			}
		}
	}
	TIMER_STOP_(PREPARE);

	/*cout << "VAL" << endl;
	for(uint16 i=0; i<nb_av; ++i)
	{
		cout << "[";
		for(uint16 j=0;j<NB_LINES_PER_BLOC; ++j)
			cout << Av[i][j] << ", ";
		cout << "]\n";
	}*/

	const int STEP=32;
	const uint8 xl = v.size () % STEP;
	const uint32 sz = v.size();
	register uint32 x=0;


	register uint32 idx;
	register uint16 val[NB_LINES_PER_BLOC];
	register unsigned char line_idx;
	register unsigned char av_idx;

	TIMER_START_(COMPUTE);
	/*for(uint32 i=0; i<v.size (); ++i)
	{
		idx = v.IndexData[i];

		for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
			val[j] = v.at(j, i);

		for(uint16 j=0; j<nb_non_zeo_lines; ++j)
		{
			line_idx = non_zero_lines_dx[j];

			for(uint16 k=0; k<nb_non_zero_Av[line_idx]; ++k)
			{
				av_idx = non_zero_Av_idx[line_idx][k];

				arrays[line_idx][idx] += (uint32)Av[line_idx][av_idx] * val[av_idx];
			}
		}
			//for(uint16 k=0; k<NB_LINES_PER_BLOC; ++k)
				//arrays[non_zero_lines_dx[j]][idx] += (uint32)(Av[non_zero_lines_dx[j]][k]) * val[k];
	}*/


	/*for(uint32 i=0; i<v.size (); ++i)
	{
		idx = v.IndexData[i];

		for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
			val[j] = v.at(j, i);

		for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
			for(uint16 k=0; k<NB_LINES_PER_BLOC; ++k)
				arrays[j][idx] += (uint32)(Av[j][k]) * val[k];
	}*/

	while(x<xl)
	{
		idx = v.IndexData[x];

		for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
			val[j] = v.at_unchecked(j, x);

#pragma loop unroll
		//for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
			//arrays[j][idx] += (uint32)Av_col1[j] * val;
		for(uint16 j=0; j<nb_non_zeo_lines; ++j)
		{
			line_idx = non_zero_lines_dx[j];

			/*for(uint16 k=0; k<nb_non_zero_Av[line_idx]; ++k)
			{
				av_idx = non_zero_Av_idx[line_idx][k];
				arrays[line_idx][idx] += (uint32)Av[line_idx][av_idx] * val[av_idx];
			}*/
			for(uint16 k=0; k<NB_LINES_PER_BLOC; ++k)
			{
				arrays[line_idx][idx] += (uint32)Av[line_idx][k] * val[k];
			}
		}

		++x;
	}

	for(x=xl; x<sz; x+=STEP)
	{

#pragma loop unroll
		for(uint8 t=0; t<STEP; ++t)
		{
			idx = v.IndexData[x+t];

			for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
				val[j] = v.at_unchecked(j, x+t);

			for(uint16 j=0; j<nb_non_zeo_lines; ++j)
			{
				line_idx = non_zero_lines_dx[j];

				/*for(uint16 k=0; k<nb_non_zero_Av[line_idx]; ++k)
				{
					av_idx = non_zero_Av_idx[line_idx][k];

					arrays[line_idx][idx] += (uint32)Av[line_idx][av_idx] * val[av_idx];
				}*/
				for(uint16 k=0; k<NB_LINES_PER_BLOC; ++k)
				{
					arrays[line_idx][idx] += (uint32)Av[line_idx][k] * val[k];
				}
			}
		}
	}

	TIMER_STOP_(COMPUTE);
}

inline void axpy(const uint16 Av_col1[],
		  const MultiLineVector<uint16>& v,
		  const uint16 line,
		  uint64 **arrays)
{
	lela_check(line < v.nb_lines());
	const uint16 STEP=32;
	const uint32 sz = v.size ();

	uint32 non_zero_Av[NB_LINES_PER_BLOC], nb_non_zero_Av=0;
	uint64 *non_zero_arrays[NB_LINES_PER_BLOC];

	for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
		if(Av_col1[j] != 0)
		{
			non_zero_Av[nb_non_zero_Av] = (uint32)Av_col1[j];
			non_zero_arrays[nb_non_zero_Av] = arrays[j];
			nb_non_zero_Av++;
		}

	/*for(uint32 i=0; i<sz; ++i)
	{
		idx = v.IndexData[i];
		val = v.at_unchecked(line, i);

#pragma loop unroll
		for(uint16 j=0; j<nb_non_zero_Av; ++j)
			non_zero_arrays[j][idx] += (uint32)non_zero_Av[j] * val;
	}*/

	const uint8 xl = v.size () % STEP;
	register uint32 x=0;
	register uint32 idx;
	register uint16 val;
	while(x<xl)
	{
		idx = v.IndexData[x];
		val = v.at_unchecked(line, x);

#pragma loop unroll
		//for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
			//arrays[j][idx] += (uint32)Av_col1[j] * val;
		for(uint16 j=0; j<nb_non_zero_Av; ++j)
			non_zero_arrays[j][idx] += (uint32)non_zero_Av[j] * val;

		++x;
	}

	for(x=xl; x<sz; x+=STEP)
	{

#pragma loop unroll
		for(uint8 t=0; t<STEP; ++t)
		{
			idx = v.IndexData[x+t];
			val = v.at_unchecked(line, x+t);

			for(uint16 j=0; j<nb_non_zero_Av; ++j)
				non_zero_arrays[j][idx] += (uint32)non_zero_Av[j] * val;
		}
	}

	/*uint8 xl = v.size () % STEP;
	register uint32 x=0;
	if(av1_col1 != 0 && av2_col1 != 0)
	{
		while(x<xl)
		{
			idx = v.IndexData[x];
			val1 = v.at_unchecked(line, x);
			arr1[idx] += av1_col1_32 * val1;
			arr2[idx] += av2_col1_32 * val1;

			++x;
		}

		for(x=xl; x<sz; x+=STEP)
		{

	#pragma loop unroll
			for(uint8 t=0; t<STEP; ++t)
			{
				idx = v.IndexData[x+t];
				val1 = v.at_unchecked(line, x+t);
				arr1[idx] += av1_col1_32 * val1;
				arr2[idx] += av2_col1_32 * val1;
			}
		}
	}
	else if(av1_col1 != 0)
	{
		while(x<xl)
		{
			idx = v.IndexData[x];
			val1 = v.at_unchecked(line, x);
			arr1[idx] += av1_col1_32 * val1;

			++x;
		}

		for(x=xl; x<sz; x+=STEP)
		{

	#pragma loop unroll
			for(uint8 t=0; t<STEP; ++t)
			{
				idx = v.IndexData[x+t];
				val1 = v.at_unchecked(line, x+t);
				arr1[idx] += av1_col1_32 * val1;
			}
		}
	} else 	//av2_col1 != 0
	{
		while(x<xl)
		{
			idx = v.IndexData[x];
			val1 = v.at_unchecked(line, x);
			arr2[idx] += av2_col1_32 * val1;

			++x;
		}

		for(x=xl; x<sz; x+=STEP)
		{

	#pragma loop unroll
			for(uint8 t=0; t<STEP; ++t)
			{
				idx = v.IndexData[x+t];
				val1 = v.at_unchecked(line, x+t);
				arr2[idx] += av2_col1_32 * val1;
			}
		}
	}*/

}

template <typename Ring>
void reducePivotsByPivots(const Ring& R, const SparseMultilineMatrix<uint16>& A, SparseMultilineMatrix<uint16>& B)
{
	lela_check(A.coldim () == B.coldim ());
	lela_check(A.rowdim () == A.coldim ());

	//typedef Modular<uint16> Ring;
	typedef SparseMultilineMatrix<uint16> Matrix;


	typename Matrix::ConstRowIterator i_A;
	typename Matrix::RowIterator i_B;

	uint32 B_coldim = B.coldim ();
	uint64 *tmpDenseArray[NB_LINES_PER_BLOC];

	for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
		tmpDenseArray[j] = new uint64[B_coldim];

	uint16 *Vals[NB_LINES_PER_BLOC];
	for(uint16 j=0; j<NB_LINES_PER_BLOC; ++j)
		Vals[j] = new uint16[NB_LINES_PER_BLOC];

	TIMER_DECLARE_(RazArrayTimer);
	TIMER_DECLARE_(CopySparseVectorToDenseArrayTimer);
	TIMER_DECLARE_(CopyDenseArrayToSparseVectorTimer);
	TIMER_DECLARE_(AxpyTimer);
	TIMER_DECLARE_(Axpy2Timer);
	TIMER_DECLARE_(AxpyTimerOuter);
	TIMER_DECLARE_(AxpyTimerLeft);

#ifdef SHOW_PROGRESS
	uint32 i=A.rowdim () / NB_LINES_PER_BLOC;
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "In spec Modular<uint16> Multiline" << std::endl;
#endif

	i_A = A.rowEnd ();
	i_B = B.rowEnd ();

	while(i_A != A.rowBegin ()) {		//for each multiline
		--i_A;
		--i_B;

		TIMER_START_(RazArrayTimer);
			 razArray64(tmpDenseArray, NB_LINES_PER_BLOC, B_coldim);
		TIMER_STOP_(RazArrayTimer);

		TIMER_START_(CopySparseVectorToDenseArrayTimer);
			copyMultilineToDenseArrays64(*i_B, tmpDenseArray, NB_LINES_PER_BLOC);
		TIMER_STOP_(CopySparseVectorToDenseArrayTimer);

		//line = A.nb_lines_per_bloc() - 1;		//last line
		//do {		//for each line
#ifdef SHOW_PROGRESS
		--i;
        report << "                                                                    \r";
        report << "\t" << i << std::ends;
#endif

		typename Ring::Element Av_col1[NB_LINES_PER_BLOC];
		uint32 Ap1;
		uint32 Ap[NB_LINES_PER_BLOC];
		int nb_consecutive = 0;

		TIMER_START_(AxpyTimerOuter);
		for(uint32 j=NB_LINES_PER_BLOC; j < i_A->size (); ++j) 		//skip first two elements
		{
			Ap1 = i_A->IndexData[j];
			Ap[0] = Ap1;

			for(int a=0; a<NB_LINES_PER_BLOC; ++a)
				memset(Vals[a], 0, NB_LINES_PER_BLOC*sizeof(uint16));

			for(uint16 k=0; k<NB_LINES_PER_BLOC; ++k)
			{
				R.copy(Av_col1[k], i_A->at_unchecked(k, j));

				if(Av_col1[k] != 0)
					R.negin(Av_col1[k]);

				Vals[k][0] = Av_col1[k];
			}

			nb_consecutive=0;
			int idx=-1;


			if(Ap1 % NB_LINES_PER_BLOC == 0 && j < i_A->size() - NB_LINES_PER_BLOC && i_A->size()>NB_LINES_PER_BLOC)
			{
				for(uint16 l=1; l<NB_LINES_PER_BLOC; ++l)
				{
					Ap[l] = i_A->IndexData[j+l];
					if((Ap[l] - Ap[0]) >= NB_LINES_PER_BLOC)
						break;

					++nb_consecutive;

					idx = Ap[l] - Ap[0];
					assert(idx >= 0 && idx < NB_LINES_PER_BLOC);

					for(uint16 k=0; k<NB_LINES_PER_BLOC; ++k)
					{
						R.copy(Vals[k][idx], i_A->at_unchecked(k, j+l));

						if(Vals[k][idx] != 0)
							R.negin(Vals[k][idx]);
					}
				}

				//cout << "NB CONS " << nb_consecutive << endl;
				j += nb_consecutive;

				if(nb_consecutive > 0)
				{
					//cout << "AXPY 2 "<< endl;
					TIMER_START_(Axpy2Timer);
						axpy2(Vals,
								nb_consecutive+1,
								B[Ap1/A.nb_lines_per_bloc()],
								tmpDenseArray);
					TIMER_STOP_(Axpy2Timer);
				}
				else
				{
					//cout << "AXPY1" << endl;
					TIMER_START_(AxpyTimer);
						axpy(Av_col1,
							 B[Ap1/A.nb_lines_per_bloc()],
							 Ap1%A.nb_lines_per_bloc(),
							 tmpDenseArray);
					TIMER_STOP_(AxpyTimer);
				}
			}
			else
			{
			TIMER_START_(AxpyTimer);
				axpy(Av_col1,
				 	 B[Ap1/A.nb_lines_per_bloc()],
					 Ap1%A.nb_lines_per_bloc(),
					 tmpDenseArray);
			TIMER_STOP_(AxpyTimer);
			}
		}
		TIMER_STOP_(AxpyTimerOuter);

		TIMER_START_(AxpyTimerLeft);
		typename Ring::Element Av;
		bool reduce = true;

		if(i_A->size () > 1)			//reduce lines within the same multiline
		{
			uint16 last_non_null_line = i_A->size () < NB_LINES_PER_BLOC ? i_A->size () : NB_LINES_PER_BLOC;

			for(uint32 j=last_non_null_line-1; j>0; --j) 		//fetch the elements from 1 to nb_lines_per_bloc
			{													//starting from the last line
				Ap1 = i_A->IndexData[j];

				for(uint32 t=0; t < B_coldim; ++t)
					tmpDenseArray[j][t] %= R._modulus;	//Make sure product in next loop doesn't overflow on line j in multiline

				//For all lines situated above line j
				for(uint16 k=0; k<j; ++k)
				{
					R.copy(Av, i_A->at_unchecked(k, j));		//all elements at index (j, >j) are null

					if(Av != 0)
					{
						R.negin(Av);

						for(uint32 t=0; t < B_coldim; ++t)
							tmpDenseArray[k][t] += (uint32)Av * tmpDenseArray[j][t];
					}
				}
			}

			for(uint32 t=0; t < B_coldim; ++t)
					tmpDenseArray[0][t] %= R._modulus;
			reduce = false;
		}
		TIMER_STOP_(AxpyTimerLeft);



		TIMER_START_(CopyDenseArrayToSparseVectorTimer);
			copyDenseArraysToMultilineVector64(R, tmpDenseArray, NB_LINES_PER_BLOC, B_coldim, *i_B, reduce);
		TIMER_STOP_(CopyDenseArrayToSparseVectorTimer);


		//} while(line != 0);		//end for each line
	}		//for each multiline

#ifdef SHOW_PROGRESS
        report << "\r                                                                    \n";
#endif

	TIMER_REPORT_(RazArrayTimer);
	TIMER_REPORT_(CopySparseVectorToDenseArrayTimer);
	TIMER_REPORT_(CopyDenseArrayToSparseVectorTimer);
	TIMER_REPORT_(AxpyTimer);
	TIMER_REPORT_(Axpy2Timer);
	TIMER_REPORT_(AxpyTimerOuter);
	TIMER_REPORT_(AxpyTimerLeft);
	TIMER_REPORT_(PREPARE);
	TIMER_REPORT_(COMPUTE);
}









bool testFaugereLachartre(const char *file_name)
{
	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtil::loadF4Modulus(file_name);
	Modular<modulus_type> R (modulus);
	Context<Ring> ctx (R);

	const uint16 nb_lines_per_bloc = NB_LINES_PER_BLOC;
	MultiLineIndexer multilineIndexer;
	Indexer<uint32> simpleIndexer;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);


	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("Loading matrix");

	commentator.start("[MultiLineIndexer] constructing indexes");
		multilineIndexer.processMatrix(A);
	commentator.stop(MSG_DONE);

	commentator.start("[SimpleIndexer] constructing indexes");
		simpleIndexer.processMatrix(A);
	commentator.stop(MSG_DONE);
	//MatrixUtil::show_mem_usage("[Indexer] constructing indexes");

	report << "Pivots found: " << multilineIndexer.Npiv << endl << endl;

	SparseMultilineMatrix<Ring::Element>	sub_A (multilineIndexer.Npiv, multilineIndexer.Npiv, nb_lines_per_bloc),
						sub_B (multilineIndexer.Npiv, A.coldim () - multilineIndexer.Npiv, nb_lines_per_bloc),
						sub_C (A.rowdim () - multilineIndexer.Npiv, multilineIndexer.Npiv, nb_lines_per_bloc),
						sub_D (A.rowdim () - multilineIndexer.Npiv, A.coldim () - multilineIndexer.Npiv, nb_lines_per_bloc);

	commentator.start("[multilineIndexer] constructing sub matrices");
		multilineIndexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, false);
	commentator.stop(MSG_DONE);
	report << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////
	SparseMatrix<Ring::Element>		sub_A_ (simpleIndexer.Npiv, simpleIndexer.Npiv),
									sub_B_ (simpleIndexer.Npiv, A.coldim () - simpleIndexer.Npiv),
									sub_C_ (A.rowdim () - simpleIndexer.Npiv, simpleIndexer.Npiv),
									sub_D_ (A.rowdim () - simpleIndexer.Npiv, A.coldim () - simpleIndexer.Npiv);

	commentator.start("[SimpleIndexer] constructing sub matrices");
		simpleIndexer.constructSubMatrices(A, sub_A_, sub_B_, sub_C_, sub_D_, false);
	commentator.stop(MSG_DONE);

////////////////////////////////////////////////////////////////////////////////////////////////////

	if(equal(R, sub_A, sub_A_))
		cout << "sub_A EQUAL" << endl;
	else
		cout << "sub_A NOT EQUAL" << endl;

	if(equal_reverse(R, sub_A, sub_A_))
		cout << "sub_A EQUAL" << endl;
	else
		cout << "sub_A NOT EQUAL" << endl;

	if(equal_reverse(R, sub_B, sub_B_))
		cout << "sub_B EQUAL" << endl;
	else
		cout << "sub_B NOT EQUAL" << endl;

	if(equal(R, sub_B, sub_B_))
		cout << "sub_B EQUAL" << endl;
	else
		cout << "sub_B NOT EQUAL" << endl;

	dumpMatrixAsPbmImage(sub_A, "sub_A_multiline.pbm");
	dumpMatrixAsPbmImage(sub_B, "sub_B_multiline.pbm");

	//MatrixUtil::dumpMatrixAsPbmImage(sub_A_, "sub_A_simple.pbm");

	commentator.start("B = A^-1 x B [multiline]");
		reducePivotsByPivots(R, sub_A, sub_B);
	commentator.stop(MSG_DONE);

	commentator.start("B = A^-1 x B [simple]");
		MatrixOp::reducePivotsByPivots(R, sub_A_, sub_B_, false, 8);
	commentator.stop(MSG_DONE);

	if(equal_reverse(R, sub_B, sub_B_))
		cout << "EQUAL" << endl;
	else
		cout << "NOT EQUAL" << endl;

	if(equal(R, sub_B, sub_B_))
		cout << "EQUAL" << endl;
	else
		cout << "NOT EQUAL" << endl;


	//dumpMatrixAsPbmImage(sub_B, "sub_B_multiline.pbm");
	//MatrixUtil::dumpMatrixAsPbmImage(sub_B_, "sub_B_simple.pbm");

	return true;
}










int main (int argc, char **argv)
{
	char *fileName = NULL;

	bool pass = true;

	static Argument args[] = {
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &fileName},
		{ '\0' }
	};

	parseArguments (argc, argv, args, "", 0);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	cout << "BLOCS OF " << NB_LINES_PER_BLOC << " LINES" << endl;
	commentator.start ("Faugère-Lachartre Bloc Version", "FaugereLachartre");

	pass = testFaugereLachartre(fileName);

	commentator.stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}


