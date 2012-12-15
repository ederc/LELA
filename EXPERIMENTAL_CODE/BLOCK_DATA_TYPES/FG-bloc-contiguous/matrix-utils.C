/*
 * matrix-util.C
 *
 *  Created on: 4 juil. 2012
 *      Author: martani
 */


#ifndef MATRIX_UTILS_C_
#define MATRIX_UTILS_C_

#include "matrix-utils.h"

using namespace LELA;
using namespace std;

template <typename Element, typename Index>
void MatrixUtils::copy(const SparseMatrix<Element>& A, SparseBlocMatrix<ContiguousBloc<Element, Index> >& B)
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
		case SparseBlocMatrix<ContiguousBloc<Element, Index> >::ArrangementTopDown_LeftRight:
			sparse_to_bloc_copyTopDowLeftRight(A, B);
		break;

		case SparseBlocMatrix<ContiguousBloc<Element, Index> >::ArrangementTopDown_RightLeft:
			throw std::logic_error ("NotImplemented");
		break;

		case SparseBlocMatrix<ContiguousBloc<Element, Index> >::ArrangementDownTop_LeftRight:
			sparse_to_bloc_copyDownTopLeftRight(A, B);
		break;
		case SparseBlocMatrix<ContiguousBloc<Element, Index> >::ArrangementDownTop_RightLeft:
			sparse_to_bloc_copyDownTopRightLeft(A, B);
		break;
	}
}


template <typename Element, typename Index>
void MatrixUtils::sparse_to_bloc_copyTopDowLeftRight(const SparseMatrix<Element>& A, SparseBlocMatrix<ContiguousBloc<Element, Index> >& B)
{
	typename SparseMatrix<Element>::ConstRowIterator i_A;
	typename SparseMatrix<Element>::Row::const_iterator it;

	uint32 row_bloc_idx, bloc_idx_in_row, line_idx_in_bloc, elt_idx_in_line;
	uint32 i;

	for(i_A = A.rowBegin (), i=0; i_A != A.rowEnd (); ++i_A, ++i)
	{
		row_bloc_idx = i / B.bloc_height ();
		line_idx_in_bloc = i % B.bloc_height ();

		if(i % B.bloc_height () == 0)			//init blocs
		{
			B[row_bloc_idx].clear ();
			uint32 nb_blocs_per_A_dim;

			if(B.isFilledWithEmptyBlocs())		//In case empty blocs should be appended
			{
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

				if(entry != (uint32)-1)			//number of blocs in this row
				{
					nb_blocs_per_A_dim = (uint32) std::ceil((float)A.coldim () / B.bloc_width ());
					nb_blocs_per_A_dim -= entry / B.bloc_width();
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

		if(i_A->empty ())
			continue;

		uint32 old_bloc_idx = (i_A->front ().first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
		uint32 nb_elts_added = 0;

		for(it = i_A->begin (); it != i_A->end (); ++it)
		{
			bloc_idx_in_row = (it->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line = (it->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			if(bloc_idx_in_row != old_bloc_idx)
			{
				B[row_bloc_idx][old_bloc_idx].set_row_size(line_idx_in_bloc, nb_elts_added);
				nb_elts_added = 0;
				old_bloc_idx = bloc_idx_in_row;
			}

			B[row_bloc_idx][bloc_idx_in_row].push_back_value(it->second);
			B[row_bloc_idx][bloc_idx_in_row].push_back_pos(elt_idx_in_line);

			++nb_elts_added;
		}

		if(bloc_idx_in_row == old_bloc_idx)
		{
			B[row_bloc_idx][old_bloc_idx].set_row_size(line_idx_in_bloc, nb_elts_added);
		}
	}

	//fix offsets
	for(uint32 i=0; i<B.rowBlocDim(); ++i)
	{
		for(uint32 j=0;j<B[i].size (); ++j)		//for each bloc in the row
		{
			B[i][j].construct_rows_starting_offset();
		}
	}
}

template <typename Element, typename SparseBlocMatrix_>
void MatrixUtils::sparse_to_bloc_copyDownTopRightLeft(const SparseMatrix<Element>& A, SparseBlocMatrix_& B)
{
	typename SparseMatrix<Element>::ConstRowIterator i_A;
	typename SparseMatrix<Element>::Row::const_iterator it;

	uint32 row_bloc_idx, bloc_idx_in_row, line_idx_in_bloc, elt_idx_in_line;
	uint32 i;

	const uint32 last_line_idx = A.rowdim ()-1;
	uint32 i_offset;

	i = A.rowdim();
	i_A = A.rowEnd ();

	while(i_A != A.rowBegin ())
	{
		--i;
		--i_A;

		i_offset = last_line_idx - i;

		row_bloc_idx = i_offset / B.bloc_height ();
		line_idx_in_bloc = i_offset % B.bloc_height ();

		if(i_offset % B.bloc_height () == 0)			//init blocs
		{
			B[row_bloc_idx].clear ();
			uint32 nb_blocs_per_A_dim;

			if(B.isFilledWithEmptyBlocs())		//In case empty blocs should be appended
			{
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

			B[row_bloc_idx].resize (nb_blocs_per_A_dim);		//init memory of nb_blocs_per_A_dim blocs

			for(uint32 k=0; k<nb_blocs_per_A_dim; ++k)			//initiliaze blocs in a row
			{
				B[row_bloc_idx][k].init(B.bloc_height(), B.bloc_width());
			}
		}

		if(i_A->empty ())
			continue;

		uint32 old_bloc_idx = (A.coldim () - 1 - i_A->front ().first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
		uint32 nb_elts_added = 0;

		it = i_A->end ();
		while(it != i_A->begin ())
		{
			--it;
			bloc_idx_in_row = (A.coldim () - 1 - it->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line = (A.coldim () - 1 - it->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			if(bloc_idx_in_row != old_bloc_idx)
			{
				B[row_bloc_idx][old_bloc_idx].set_row_size(line_idx_in_bloc, nb_elts_added);
				nb_elts_added = 0;
				old_bloc_idx = bloc_idx_in_row;
			}

			B[row_bloc_idx][bloc_idx_in_row].push_back_value(it->second);
			B[row_bloc_idx][bloc_idx_in_row].push_back_pos(elt_idx_in_line);

			++nb_elts_added;
		}

		if(bloc_idx_in_row == old_bloc_idx)
		{
			B[row_bloc_idx][old_bloc_idx].set_row_size(line_idx_in_bloc, nb_elts_added);
		}
	}

	//fix offsets
	for(uint32 i=0; i<B.rowBlocDim(); ++i)
	{
		for(uint32 j=0;j<B[i].size (); ++j)		//for each bloc in the row
		{
			B[i][j].construct_rows_starting_offset();
		}
	}
}


template <typename Element, typename Index>
void MatrixUtils::sparse_to_bloc_copyDownTopLeftRight(const SparseMatrix<Element>& A, SparseBlocMatrix<ContiguousBloc<Element, Index> >& B)
{
	typename SparseMatrix<Element>::ConstRowIterator i_A;
	typename SparseMatrix<Element>::Row::const_iterator it;

	uint32 row_bloc_idx, bloc_idx_in_row, line_idx_in_bloc, elt_idx_in_line;
	uint32 i;

	const uint32 last_line_idx = A.rowdim ()-1;
	uint32 i_offset;

	i = A.rowdim();
	i_A = A.rowEnd ();

	while(i_A != A.rowBegin ())
	{
		--i;
		--i_A;

		i_offset = last_line_idx - i;

		row_bloc_idx = i_offset / B.bloc_height ();
		line_idx_in_bloc = i_offset % B.bloc_height ();

		if(i_offset % B.bloc_height () == 0)			//init blocs
		{
			B[row_bloc_idx].clear ();
			uint32 nb_blocs_per_A_dim;

			if(B.isFilledWithEmptyBlocs())		//In case empty blocs should be appended
			{
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
			
			B[row_bloc_idx].resize (nb_blocs_per_A_dim);		//init memory of nb_blocs_per_A_dim blocs

			//TODO: initialize only blocs where there is actual data
			for(uint32 k=0; k<nb_blocs_per_A_dim; ++k)			//initiliaze blocs in a row
				B[row_bloc_idx][k].init(B.bloc_height(), B.bloc_width());
		}

		if(i_A->empty ())
			continue;

		uint32 old_bloc_idx = (i_A->front ().first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
		uint32 nb_elts_added = 0;

		for(it = i_A->begin (); it != i_A->end (); ++it)
		{
			bloc_idx_in_row = (it->first - B.FirstBlocsColumIndexes[row_bloc_idx]) / B.bloc_width();
			elt_idx_in_line = (it->first - B.FirstBlocsColumIndexes[row_bloc_idx]) % B.bloc_width();

			if(bloc_idx_in_row != old_bloc_idx)
			{
				B[row_bloc_idx][old_bloc_idx].set_row_size(line_idx_in_bloc, nb_elts_added);
				nb_elts_added = 0;
				old_bloc_idx = bloc_idx_in_row;
			}

			B[row_bloc_idx][bloc_idx_in_row].push_back_value(it->second);
			B[row_bloc_idx][bloc_idx_in_row].push_back_pos(elt_idx_in_line);

			++nb_elts_added;
		}

		if(bloc_idx_in_row == old_bloc_idx)
		{
			B[row_bloc_idx][old_bloc_idx].set_row_size(line_idx_in_bloc, nb_elts_added);
		}
	}

	//fix offsets
	for(uint32 i=0; i<B.rowBlocDim(); ++i)
	{
		for(uint32 j=0;j<B[i].size (); ++j)		//for each bloc in the row
		{
			B[i][j].construct_rows_starting_offset();
		}
	}

}


template <typename Element, typename Index>
void MatrixUtils::copy(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, SparseMatrix<Element>& B)
{
	check_equal_or_raise_exception(A.coldim (), B.coldim ());
	check_equal_or_raise_exception(A.rowdim (), B.rowdim ());

	for(uint32 j=0; j<B.rowdim (); ++j)
		B[j].clear ();

	switch(A.blocArrangement)
	{
		case SparseBlocMatrix<ContiguousBloc<Element, Index> >::ArrangementTopDown_LeftRight:
			bloc_to_sparse_copyTopDowLeftRight(A, B);
		break;

		case SparseBlocMatrix<ContiguousBloc<Element, Index> >::ArrangementTopDown_RightLeft:
			throw std::logic_error ("NotImplemented");
		break;

		case SparseBlocMatrix<ContiguousBloc<Element, Index> >::ArrangementDownTop_LeftRight:
			bloc_to_sparse_copyDownTopLeftRight(A, B);
		break;
		case SparseBlocMatrix<ContiguousBloc<Element, Index> >::ArrangementDownTop_RightLeft:
			bloc_to_sparse_copyDownTopRightLeft(A, B);
		break;
	}
}

template <typename Element, typename Index>
void MatrixUtils::bloc_to_sparse_copyTopDowLeftRight(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, SparseMatrix<Element>& B)
{
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

			for(uint16 k=0; k<A.bloc_height () && k<(B.rowdim() - curr_row_B); ++k)	//for each row in the bloc
			{
				for (uint32 l = 0; l < A[i][j].get_row_size(k); ++l)
				{
					B[curr_row_B + k].push_back(
							typename SparseMatrix<Element>::Row::value_type(
									bloc_idx + A[i][j].pos_in_row(k, l),
									A[i][j].value_in_row(k, l)));
				}
			}
		}

		curr_row_B += A.bloc_height();
	}
}


template <typename Element, typename Index>
void MatrixUtils::bloc_to_sparse_copyDownTopLeftRight(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, SparseMatrix<Element>& B)
{
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
			{
				continue;
			}

			uint32 bloc_idx = A.FirstBlocsColumIndexes[i] + (A.bloc_width() * j);

			for(uint16 k=0; k<A.bloc_height () && k <= curr_row_B; ++k)	//for each row in the bloc
			{
				for (uint32 l = 0; l < A[i][j].get_row_size(k); ++l)
				{
					B[curr_row_B-k].push_back(
							typename SparseMatrix<Element>::Row::value_type(
									bloc_idx + A[i][j].pos_in_row(k, l),
									A[i][j].value_in_row(k, l)));
				}
			}
		}

		curr_row_B -= A.bloc_height();
	}
}


template <typename Element, typename Index>
void MatrixUtils::bloc_to_sparse_copyDownTopRightLeft(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, SparseMatrix<Element>& B)
{
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

			uint32 bloc_idx = B.coldim() - 1 - A.FirstBlocsColumIndexes[i] - A.bloc_width() * j;		//the colum index of the bloc starting from *THE LEFT*

			for(uint16 k=0; k<A.bloc_height () && k <= curr_row_B; ++k)	//for each row in the bloc
			{
				for (int l = A[i][j].get_row_size(k)-1; l >= 0; --l)
				{
					B[curr_row_B-k].push_back(
							typename SparseMatrix<Element>::Row::value_type(
									(uint32)(bloc_idx - A[i][j].pos_in_row(k, l)),
									A[i][j].value_in_row(k, l)));
				}
			}
		}

		curr_row_B -= A.bloc_height();
	}
}

template <typename Element, typename Index>
void MatrixUtils::dumpMatrixAsPbmImage(const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, const char *outputFileName)
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

template <typename Element, typename Index>
bool MatrixUtils::equal(const Modular<Element>& R, const SparseBlocMatrix<ContiguousBloc<Element, Index> >& A, const SparseMatrix<Element>& B)
{
	SparseMatrix<Element> tmp (A.rowdim(), A.coldim());
	copy(A, tmp);

	Context<Modular<Element> > ctx (R);

	return BLAS3::equal(ctx, B, tmp);
}


template <typename Element, typename Index>
bool MatrixUtils::equal(const Modular<Element>& R, const SparseMatrix<Element>& A, const SparseBlocMatrix<ContiguousBloc<Element, Index> >& B)
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

		assert(sz <= m);
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
			assert(pos[j] < m);
		}
	}

	//free memory
	delete[] oNz;
	delete[] oPos;
	fclose(f);
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
