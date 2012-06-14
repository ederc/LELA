/*
 * test-FG-multiline.C
 *
 *  Created on: 14 juin 2012
 *      Author: martani
 */

#include <iostream>

#include "FG-types.h"
#include "../new-faugere-lachartre/matrix-util.h"
#include "../new-faugere-lachartre/indexer.h"
#include "multiline-indexer.h"

#include "../util/support.h"

using namespace LELA;
using namespace std;

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
	sprintf(buffer, "P4\n# matrix size(%lu, %lu)\n%lu %lu\n", A.rowdim (), A.coldim (), A.coldim (), A.rowdim ());
	fwrite(buffer, sizeof(char), strlen(buffer), outStream);

	while(i_A != A.rowEnd())
	{
		row_size = i_A->size ();

		for(uint16 i=0; i < A.nb_lines_per_bloc (); ++i)
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
		}

		if(col%8 != 0)
			fwrite(&output_byte, sizeof(unsigned char), 1, outStream);

		fflush(outStream);

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

	uint32 t=0;

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

			if(!BLAS1::equal(ctx, *i_B, v))
			{
				uint32 t=0;
				for(uint32 j=0; j < v.back ().first+1; ++j)
				{
					if(j != v[t].first)
						cout << "0, ";
					else
					{
						cout << v[t].second << ", ";
						++t;
					}
				}
				cout << endl << endl;
				BLAS1::write(ctx, cout, *i_B);
				cout << endl << endl;

				for(uint32 j=0; j < i_B->size (); ++j)
				{
					cout << (*i_B)[j].first << ", ";
				}

				cout << endl << endl;

				for(uint32 j=0; j < v.size (); ++j)
				{
					cout << v[j].first << ", ";
				}

				cout << endl << endl;



				for(uint32 j=0; j < i_B->size (); ++j)
				{
					if((*i_B)[j].first != v[j].first)
					{
						cout << "DIFF ON " << (*i_B)[j].first << "      " << v[j].first;
						return false;
					}
				}


				return false;
			}

			v.clear ();
			++i_B;
			++t;
		}

		++i_A;
	}

	return true;
}

bool testFaugereLachartre(const char *file_name)
{
	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtil::loadF4Modulus(file_name);
	Modular<modulus_type> R (modulus);
	Context<Ring> ctx (R);

	MultiLineIndexer indexer;
	Indexer<uint32> simpleIndexer;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);


	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("Loading matrix");

	commentator.start("[Indexer] constructing indexes");
		indexer.processMatrix(A);
		simpleIndexer.processMatrix(A);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("[Indexer] constructing indexes");

	report << "Pivots found: " << indexer.Npiv << endl << endl;

	SparseMultilineMatrix<Ring::Element>	sub_A (indexer.Npiv, indexer.Npiv),
											sub_B (indexer.Npiv, A.coldim () - indexer.Npiv),
											sub_C (A.rowdim () - indexer.Npiv, indexer.Npiv),
											sub_D (A.rowdim () - indexer.Npiv, A.coldim () - indexer.Npiv);

	commentator.start("[Indexer] constructing sub matrices");
		indexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, false);
	commentator.stop(MSG_DONE);

////////////////////////////////////////////////////////////////////////////////////////////////////
	SparseMatrix<Ring::Element>		sub_A_ (indexer.Npiv, indexer.Npiv),
									sub_B_ (indexer.Npiv, A.coldim () - indexer.Npiv),
									sub_C_ (A.rowdim () - indexer.Npiv, indexer.Npiv),
									sub_D_ (A.rowdim () - indexer.Npiv, A.coldim () - indexer.Npiv);

	commentator.start("[Indexer2] constructing sub matrices");
		simpleIndexer.constructSubMatrices(A, sub_A_, sub_B_, sub_C_, sub_D_, false);
	commentator.stop(MSG_DONE);

////////////////////////////////////////////////////////////////////////////////////////////////////

	if(equal(R, sub_A, sub_A_))
		cout << "EQUAL" << endl;
	else
		cout << "NOT EQUAL" << endl;

	dumpMatrixAsPbmImage(sub_A, "sub_A.pbm");

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

	commentator.start ("Faugère-Lachartre Bloc Version", "FaugereLachartre");

	pass = testFaugereLachartre(fileName);

	commentator.stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}


