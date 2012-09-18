/*
 * matrix-util.C
 * Copyright 2012 Martani Fayssal (LIP6 / UPMC University Paris06)
 *
 *  Created on: 30 mai 2012
 *      Author: martani (LIP6 / UPMC University Paris06)
 */


#ifndef MATRIX_UTIL_C_
#define MATRIX_UTIL_C_

#include <assert.h>
//#include <boost/functional/hash.hpp>
 #include <unistd.h>
#include <sys/resource.h>

#include "matrix-util.h"

using namespace LELA;

/*void MatrixUtil::show_mem_usage(std::string msg)
{
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "\t[[[ " << msg << " ]]]\t\t" << " Memory : " << ru.ru_maxrss / 1024.0 << "MB" << std::endl;
	//printf(" Using: %.2fMB of memory\n", );
}*/

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
	if(fread(&mod, sizeof(uint32),     1,f) != 1)
	{
		report << "Error while reading file " << fileName << std::endl;
		throw "Error while reading file";
	}

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

	if(fread(&n, sizeof(uint32),       1,f) != 1)
		throw "Error while reading file";
	if(fread(&m, sizeof(uint32),       1,f) != 1)
		throw "Error while reading file";
	if(fread(&mod, sizeof(uint32),     1,f) != 1)
		throw "Error while reading file";
	if(fread(&nb, sizeof(uint64),1,f) != 1)
		throw "Error while reading file";

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
		if(fread(&sz, sizeof(uint32), 1, f) != 1 )
			throw "Error while reading file";

		row_sizes_offset += sizeof(uint32);
                                
		assert(sz <= m);		//number of elements in a row at max equal to the size of a row in the matrix

		//read sz elements from the values part of the file
		fseek(f, row_values_offset, SEEK_SET);
		if(fread(nz, sizeof(uint16), sz, f) != sz )
			throw "Error while reading file";

		row_values_offset += sz*sizeof(uint16);

		//read sz elements from the posistions part of the file
		fseek(f, row_positions_offset, SEEK_SET);
		if(fread(pos, sizeof(uint32), sz, f) != sz )
			throw "Error while reading file";

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

//Expects a SparseMatrix
template <class Ring, typename Matrix>
void MatrixUtil::writeF4MatrixToFile(const Ring &R, const char *fileName, const Matrix& A)
{
	unsigned int n;
	unsigned int m;
	unsigned int  mod;
	unsigned long long  nb;

	n = A.rowdim ();
	m = A.coldim ();
	mod = R._modulus;

	uint32 i;
	nb = 0;
	for (i = 0; i < n; ++i) {
		nb += A[i].size ();
	}

	FILE *f = fopen(fileName, "w");
	if (f == NULL)
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) 
                        << "Can't open " << fileName << std::endl;
		throw "Can't open file";	//TODO: handle properly with exception classes
	}

	fwrite(&n, sizeof(unsigned int), 1, f);
	fwrite(&m, sizeof(unsigned int), 1, f);
	fwrite(&mod, sizeof(unsigned int), 1, f);
	fwrite(&nb, sizeof(unsigned long long), 1, f);

	typename Matrix::Row::const_iterator it;
	//write values
	for (i = 0; i < n; ++i) {
		it = A[i].begin ();
		while(it != A[i].end ())
		{
			fwrite(&(it->second), sizeof(unsigned short int), 1, f);
			++it;
		}
	}

	//write positions
	for (i = 0; i < n; ++i) {
		it = A[i].begin ();
		while(it != A[i].end ())
		{
			fwrite(&(it->first), sizeof(unsigned int), 1, f);
			++it;
		}
	}

	size_t sz;
	//write rows' sizes
	for (i = 0; i < n; ++i) {
		sz = A[i].size ();
		fwrite(&sz, sizeof(unsigned int), 1, f);
	}

	fclose(f);
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

//given a path to a file, returns prefix|parent directory|__|file|suffix	where "|" denotes concatenation
std::string MatrixUtil::getOutputFileNameWithExtension(const char *inputfile, const char* prefix, const char *suffix)
{
	std::string f (inputfile), pre (prefix), suff (suffix);
	std::string file, dir;

	size_t found;
	found = f.find_last_of("/");

	file = f.substr(found+1);
	dir = f.substr(0, found);

	found = dir.find_last_of("/");
	dir = dir.substr(found+1);

	std::stringstream ss;
	ss << prefix << dir << "__" << file << suffix;

	return ss.str ();
}

template <typename Matrix, typename Ring>
bool MatrixUtil::verifyMatrixRowsAreUnitary(Ring& R, const Matrix& A, typename Ring::Element& det)
{
	typename Matrix::RowIterator i_A = A.rowBegin ();
	R.init(det, 1);

	while(i_A != A.rowEnd ())
	{
		if(!i_A->empty ())
		{
			R.mulin(det, i_A->front ().second);
		}
		++i_A;
	}

	if(det != 1)
			return false;
		else
			return true;
}

template <typename Matrix, typename Ring>
bool MatrixUtil::verifyEveryRowEntryIsNoGreaterThanThePreceedingRows(Ring& R, const Matrix& A)
{
	typename Ring::Element e;
	Context<Ring> ctx (R);

	uint32 h, last_h = BLAS1::head(ctx, e, A[0]);;
	for (uint32 i = 1; i < A.rowdim (); ++i) {
		h = BLAS1::head(ctx, e, A[i]);
		if(h  > last_h)
		{
			return false;
		}

		last_h = h;
	}

	return true;
}


/* ACTIVATE this if you have the boost library */
//verifies the equality of matrices by the hash of their rows (if rows are not in the same order, the matrices are still equal)
/*template <typename Matrix, typename Ring>
bool MatrixUtil::matrixEqualUsingHash(Ring& R, const Matrix& A, const Matrix& B)
{
	assert(A.coldim () == B.coldim ());
	assert(A.rowdim () == B.rowdim ());

	std::size_t res1 = 0, res2 = 0;
	typename Ring::Element arr[A.coldim ()];

	typename Matrix::Row::const_iterator it;

	for (uint32 i = 0; i < A.rowdim (); ++i) {
		memset(arr, 0, A.coldim () * sizeof(typename Ring::Element));

		for (it = A[i].begin(); it != A[i].end (); ++it) {
			arr[it->first] = it->second;
		}

		res1 ^= hasharray(arr, A.coldim ());
	}

	for (uint32 i = 0; i < B.rowdim (); ++i) {
		memset(arr, 0, A.coldim () * sizeof(typename Ring::Element));

		for (it = B[i].begin(); it != B[i].end (); ++it) {
			arr[it->first] = it->second;
		}

		res2 ^= hasharray(arr, A.coldim ());
	}

	return res1 == res2;
}*/

/* ACTIVATE this if you have the boost library */
/*template <typename T>
std::size_t MatrixUtil::hasharray(const T arr[], int N)
{
	 return boost::hash_range(arr, arr+N);
}*/

template <typename Matrix>
void MatrixUtil::invertMatrixRows(Matrix& A)
{
	uint32 i, rowdim = A.rowdim ()-1;
	for (i = 0; i <  (rowdim+1)/2; ++i) {
		std::swap(A[i], A[rowdim - i]);
	}
}


// For each r in A: r <- entry(r)^-1 * r
template <typename Ring, typename Matrix>
void MatrixUtil::makeRowsUnitary(const Ring& R, Matrix& A)
{
	typename Matrix::RowIterator i_A;
	typename Matrix::Row::iterator it;
	typename Ring::Element inv;

	Context<Ring> ctx (R);

	for (i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A)
	{
		if(i_A->empty ())
			continue;

		it = i_A->begin ();
		if(R.inv(inv, it->second) != true)		//should be invertible
			throw "Non invertible value";

		BLAS1::scal(ctx, inv, *i_A);
	}
}

template <typename Ring>
SparseMatrix<typename Ring::Element> MatrixUtil::generateIDMatrix(const Ring R, size_t size)
{
	SparseMatrix<typename Ring::Element> ID (size, size);
	for (size_t i = 0;  i < size; ++ i) {
		ID[i].push_back (typename Vector<Ring>::Sparse::value_type (i, 1));
	}

	return ID;
}

template <typename Matrix>
void MatrixUtil::freeMatrixMemory(Matrix& A)
{
	typename Matrix::RowIterator i_M = A.rowBegin ();

	while(i_M != A.rowEnd ())
	{
		(*i_M).free ();
		++i_M;
	}
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


template <typename Matrix>
std::pair<uint64, double> MatrixUtil::getMatrixSizeAndDensity(const Matrix& A)
{
	uint64 nb_elts = 0;
        
        typename Matrix::ConstRowIterator i_A = A.rowBegin ();
        
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

#endif // MATRIX_UTIL_C_
