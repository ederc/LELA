/*
 * test-structured-gauss.C
 * Copyright 2012 Martani Fayssal (LIP6 / UPMC University Paris06)
 *
 *  Created on: 24 mai 2012
 *      Author: martani  (LIP6 / UPMC University Paris06)
 * 
 * Tests the structured Gaussian elmination method
 * run with; ./test-structured-gauss - -f path/to/matrix/file
 */

#include "matrix-util.h"
#include "structured-gauss-lib.h"
#include "../util/support.h"

#include "lela/util/commentator.h"
#include "lela/ring/modular.h"
#include "lela/matrix/sparse.h"

using namespace LELA;
using namespace std;

int main (int argc, char **argv)
{
	char *file_name = NULL;

	static Argument args[] = {
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &file_name},
		{ '\0' }
	};

	parseArguments (argc, argv, args, "", 0);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator.getMessageClass (PROGRESS_REPORT).setMaxDepth (3);

	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;
	size_t rank;

	modulus_type modulus = MatrixUtil::loadF4Modulus(file_name);
	Ring R (modulus);

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	
	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
	commentator.stop(MSG_DONE);

	
	/*commentator.start("Writing matrix to file");
		report << "\t" << getOutputFileNameWithExtension(file_name, strdup("pbm/"), strdup("-orig.pbm")).c_str() << endl;
		dumpMatrixAsPbmImage(A, getOutputFileNameWithExtension(file_name, strdup("pbm/"), strdup("-orig.pbm")).c_str());
	commentator.stop(MSG_DONE);
	report << endl;*/

	commentator.start("Structured Gauss _ rref", "STRUCTURED GAUSS");
		rank = StructuredGauss::echelonize_reduced(R, A);
		report << "Rank " << rank << endl;
	commentator.stop(MSG_DONE, "STRUCTURED GAUSS");

	//std::pair<uint64, double> p = MatrixUtil::getMatrixSizeAndDensity(A);
	//report << endl;
	//report << "Nb of Nz elements: " << p.first << " (density " << p.second << ")" << std::endl;

	/*commentator.start("Writing matrix to file");
		report << "\t" << getOutputFileNameWithExtension(file_name, strdup("pbm/"), strdup("-structured_rref.pbm")).c_str() << endl;
		dumpMatrixAsPbmImage(A, getOutputFileNameWithExtension(file_name, strdup("pbm/"), strdup("-structured_rref.pbm")).c_str());
	commentator.stop(MSG_DONE);*/


	return 0;
}

