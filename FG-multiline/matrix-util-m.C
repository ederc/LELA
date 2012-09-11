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

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"

using namespace std;
using namespace LELA;

namespace NS
{

template<class Ring>
static void loadF4Matrix__low_memory_syscall_no_checks(const Ring &R,
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

	
template<typename Element>
static void write(const MultiLineVector<Element>& v, uint32 size)
{
	cout << "[";
	uint32 j=0;
	for(uint32 i=0; i<size; ++i)
	{
		if(j<v.size() && v.IndexData[j] == i)
		{
			cout << " " << v.at(0, j);
			++j;
		}
		else
			cout << " 0";
	}
	cout << "]" << endl;

	j=0;
	cout << "[";
	for(uint32 i=0; i<size; ++i)
	{
		if(j<v.size() && v.IndexData[j] == i)
		{
			cout << " " << v.at(1, j);
			++j;
		}
		else
			cout << " 0";
	}
	cout << "]" << endl;
}

template<typename Element>
static void write(const SparseMultilineMatrix<Element>& A)
{
	typename SparseMultilineMatrix<Element>::ConstRowIterator i_A = A.rowBegin();

	cout << "Matrix " << endl;
	while(i_A != A.rowEnd())
	{
		write(*i_A, A.coldim());
		++i_A;
	}
	cout << endl;
}


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
static bool equal(Modular<Element>& R, const SparseMultilineMatrix<Element>& A, const Matrix& B, bool show_diff)
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
				if(show_diff)
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
bool equal_reverse(Modular<Element>& R, const SparseMultilineMatrix<Element>& A, const Matrix& B, bool show_diff)
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
				if(show_diff)
				{
					cout << "MULTI_LINE" << endl;
					BLAS1::write(ctx, cout, v);
					cout << endl << endl;
					cout << "SPARSE VECTOR" << endl;
					BLAS1::write(ctx, cout, *i_B);
					cout << endl << endl;

					cout << "Diff on line: " << line;
					cout << endl << endl;
				}
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
