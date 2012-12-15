/*
 * test-splicer-indexer.C
 *
 *  Created on: 30 mai 2012
 *      Author: martani (LIP6 / UPMC University Paris06)
 *
 */


#include <iostream>
#include <ctime>
#include <cmath>

#include "matrix-util.h"
#include "indexer.h"
#include "LELA-FG-util.C"

#include "../util/support.h"

#include "lela/util/commentator.h"
#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"

using namespace LELA;
using namespace std;

void testLELASplicerVsIndexer(const char *file_name)
{
	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtil::loadF4Modulus(file_name);
	Ring R (modulus);
	Context<Ring> ctx (R);
	Indexer<uint32> indexer;

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
        MatrixUtil::invertMatrixRows(A);	//accomodate format to LELA's splicer
	commentator.stop(MSG_DONE);
        
        commentator.start("Constructing submatrices using Splicer", "USING LELA SPLICER");
                LELA_GF_UTIL::spliceMatrix(R, A);
        commentator.stop(MSG_DONE, "USING LELA SPLICER");
        
        report << endl;
        commentator.start("Constructing submatrices using Indexer", "USING INDEXER");
                commentator.start("Constructing indexes");
                        indexer.processMatrix(A);
                commentator.stop(MSG_DONE);

                SparseMatrix<Ring::Element>    sub_A (indexer.Npiv, indexer.Npiv),
                                                        sub_B (indexer.Npiv, A.coldim () - indexer.Npiv),
                                                        sub_C (A.rowdim () - indexer.Npiv, indexer.Npiv),
                                                        sub_D (A.rowdim () - indexer.Npiv, A.coldim () - indexer.Npiv);

                commentator.start("Constructing sub matrices");
                        indexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, true);
                commentator.stop(MSG_DONE);
        
        commentator.stop(MSG_DONE, "USING INDEXER");
}


int main (int argc, char **argv)
{
	char *fileName = NULL;
	static Argument args[] = {
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &fileName},
		{ '\0' }
	};

	parseArguments (argc, argv, args, "", 0);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	testLELASplicerVsIndexer(fileName);
	
	return 0;
}
