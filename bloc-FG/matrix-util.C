/*
 * matrix-util.C
 *
 *  Created on: 12 juin 2012
 *      Author: martani
 */

#ifndef MATRIX_UTIL_C_
#define MATRIX_UTIL_C_

#include "matrix-util.h"
#include "bloc-types.h"

void MatrixUtil::show_mem_usage(std::string msg)
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

uint32 MatrixUtil::loadF4Modulus(const char *fileName)
{
	uint32 mod;

    std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);

	FILE *f = fopen(fileName, "r");
	if (f == NULL)
	{
		report << "Can't open " << fileName << std::endl;
		exit(1);
	}

	fseek(f, 2 * sizeof(uint32), SEEK_SET);
	assert(fread(&mod, sizeof(uint32),     1,f) == 1);
	assert(mod >= 2);

	fclose(f);
	return mod;
}

//reads the matrix row by row from the file, does not load the whole file to memory. more efficient than dump_matrix.c
template <class Ring>
SparseMatrix<typename Ring::Element> MatrixUtil::loadF4Matrix(const Ring &R, const char *fileName)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	// Code adapted from C version of dump_matrix.c

	FILE *f = fopen(fileName, "r");
	if (f == NULL)
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                        << "Can't open " << fileName << std::endl;
		throw "Can't open file";
	}

	uint16 *nz;
	uint32       *pos;
	uint32        sz;
	uint32 n;
	uint32 m;
	uint32  mod;
	uint64  nb;

	assert(fread(&n, sizeof(uint32),       1,f) == 1);
	assert(fread(&m, sizeof(uint32),       1,f) == 1);
	assert(fread(&mod, sizeof(uint32),     1,f) == 1);
	assert(fread(&nb, sizeof(uint64),1,f) == 1);

	report << n << " x " << m << " matrix" << std::endl;
	report << "mod " << mod << std::endl;
	{
		double Nz=(double)(n)*(double)(m);
		Nz=(double)(nb)/Nz;
		Nz*=100.0;
		report << "Nb of Nz elements " << nb << " (density " << Nz << ")" << std::endl;
	}

	SparseMatrix<typename Ring::Element> A (n, m);
	typename SparseMatrix<typename Ring::Element>::RowIterator i_A;

	uint32 i;
	nz = new unsigned short int [m];	//has a size of at most a full row of the matrix
	pos = new uint32 [m];

	//save the arrays original pointers
	uint16 *oNz = nz;
	uint32 *oPos = pos;

	uint32 header_size = sizeof(uint32) * 3 + sizeof(uint64);	//size of n, m, mod and nb in the header of the file
	uint64 row_sizes_offset, row_values_offset, row_positions_offset;		//cursors in the file

	//row sizes if positioned after the values and the positions of the elements in the file
	row_sizes_offset = nb*sizeof(uint16) + nb*sizeof(uint32) + header_size;
	row_values_offset = header_size;
	row_positions_offset = nb*sizeof(uint16) + header_size;

	for(i_A = A.rowBegin (), i=0; i<n; i++, ++i_A){
		//get the size of the current row
		fseek(f, row_sizes_offset, SEEK_SET);
		assert(fread(&sz, sizeof(uint32), 1, f) == 1 );
		row_sizes_offset += sizeof(uint32);

		assert(sz <= m);		//number of elements in a row at max equal to the size of a row in the matrix

		//read sz elements from the values part of the file
		fseek(f, row_values_offset, SEEK_SET);
		assert(fread(nz, sizeof(uint16), sz, f) == sz );
		row_values_offset += sz*sizeof(uint16);

		//read sz elements from the posistions part of the file
		fseek(f, row_positions_offset, SEEK_SET);
		assert(fread(pos, sizeof(uint32), sz, f) == sz );
		row_positions_offset += sz*sizeof(uint32);

		i_A->reserve (sz);
		for(uint32 j=0; j<sz; j++)
		{
			i_A->push_back (typename Vector<Ring>::Sparse::value_type (pos[j], typename Ring::Element ()));
			R.init(i_A->back ().second, nz[j]);
		}
	}

	//free memory
	delete[] oNz; delete[] oPos;
	fclose(f);

	return A;
}

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void MatrixUtil::process_mem_usage(double& vm_usage, double& resident_set)
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

//dumps the matrix content as a PBM image. null elements are represented with a white pixel
//other elements are represented with a black pixel

template <typename Matrix>
void MatrixUtil::dumpMatrixAsPbmImage(const Matrix& A, const char *outputFileName)
{
	typename Matrix::ConstRowIterator i_A = A.rowBegin ();
	uint32 j, col, m, row_size;
	char buffer[512];
	unsigned char output_byte = 0;
	m = A.coldim ();

	FILE *outStream = fopen(outputFileName, "wb");

	//magic PBM header
	sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", A.rowdim (), A.coldim (), A.coldim (), A.rowdim ());
	fwrite(buffer, sizeof(char), strlen(buffer), outStream);


	while(i_A != A.rowEnd())
	{
		row_size = (*i_A).size ();

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

//dumps the matrix content as a PBM image. null elements are represented with a white pixel
//other elements are represented with a black pixel

template <typename Element, typename Index>
void dumpMatrixAsPbmImage(const FG_Types::SparseBlocMatrix<uint16, uint16>& A, const char *outputFileName)
{
	typedef FG_Types::SparseBlocMatrix<uint16, uint16> _SparseBlocMatrix;

	typename _SparseBlocMatrix::ConstRowIterator i_A = A.rowBegin ();
	//typename _SparseBlocMatrix::RowOfBlocs::

	uint32 j, col, m, row_size;
	char buffer[512];
	unsigned char output_byte = 0;
	m = A.coldim ();

	FILE *outStream = fopen(outputFileName, "wb");

	//magic PBM header
	sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", A.rowdim (), A.coldim (), A.coldim (), A.rowdim ());
	fwrite(buffer, sizeof(char), strlen(buffer), outStream);


	while(i_A != A.rowEnd())
	{
		row_size = (*i_A).size ();

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

#endif /* MATRIX_UTIL_H_ */


