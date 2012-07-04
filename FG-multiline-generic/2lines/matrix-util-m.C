/*
 * matrix-util.C
 *
 *  Created on: 19 juin 2012
 *      Author: martani
 */

#ifndef MATRIX_UTIL_MULTILINE_
#define MATRIX_UTIL_MULTILINE_

#include <iostream>
#include "FG-types.h"

#include "lela/ring/modular.h"

using namespace std;
using namespace LELA;

namespace NS
{

template<typename Element>
static void dumpMatrixAsPbmImage(const SparseMultilineMatrix<Element>& A, const char *outputFileName)
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

template <typename Matrix, typename Element>
static bool equal(Modular<Element>& R, SparseMultilineMatrix<Element>& A, Matrix& B)
{
	SparseVector<Element> v;
	typename SparseMultilineMatrix<Element>::ConstRowIterator i_A = A.rowBegin();
	typename Matrix::ConstRowIterator i_B = B.rowBegin();

	Context<Modular<Element> > ctx(R);

	if (A.rowdim() != B.rowdim())
		return false;

	uint32 line = 0;

	while (i_A != A.rowEnd()) //for each multiline of A
	{
		for (uint16 i = 0; i < A.nb_lines_per_bloc(); ++i) //for each line in the multiline
		{
			//cout << "Line " << t << " size: " << i_A->size () << endl << "[";
			for (uint32 j = 0; j < i_A->size(); ++j)
			{
				if (i_A->at(i, j) != 0)
				{
					v.push_back(
							typename Matrix::Row::value_type(i_A->IndexData[j],
									i_A->at(i, j)));
				}
			}

			if (v.empty() && i_B == B.rowEnd())
				return true;

			++line;
			if (!BLAS1::equal(ctx, *i_B, v))
			{
				BLAS1::write(ctx, cout, v);
				cout << endl << endl;
				BLAS1::write(ctx, cout, *i_B);
				cout << endl << endl;

				cout << "Diff on line: " << line;
				cout << endl << endl;

				for (uint32 j = 0; j < i_B->size(); ++j)
				{
					if ((*i_B)[j].first != v[j].first)
					{
						cout << "DIFF ON Index: " << (*i_B)[j].first << "      "
								<< v[j].first << endl;
						return false;
					}
					else if ((*i_B)[j].second != v[j].second)
					{
						cout << "DIFF ON Value: " << (*i_B)[j].second
								<< "      " << v[j].second << " At index: "
								<< (*i_B)[j].first << endl;
						return false;
					}

				}

				return false;
			}
			else
				; //cout << "Line " << line << " OK " << endl;

			v.clear();
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





}




#endif
