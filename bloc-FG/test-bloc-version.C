#include <iostream>

#include "bloc-types.h"
#include "matrix-util.h"
#include "indexer.h"
#include "../util/support.h"

using namespace LELA;
using namespace std;
using namespace FG_Types;




bool testFaugereLachartre(const char *file_name)
{
	typedef uint16 modulus_type;
	typedef Modular<modulus_type> Ring;

	modulus_type modulus = MatrixUtil::loadF4Modulus(file_name);
	Modular<modulus_type> R (modulus);
	Context<Ring> ctx (R);

	Indexer<uint32> indexer;


	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);


	commentator.start("Loading matrix");
		SparseMatrix<Ring::Element> A = MatrixUtil::loadF4Matrix(R, file_name);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("Loading matrix");

	commentator.start("[Indexer] constructing indexes");
		indexer.processMatrix(A);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("[Indexer] constructing indexes");

	report << "Pivots found: " << indexer.Npiv << endl << endl;

	SparseBlocMatrix<Ring::Element, uint16>	sub_A (indexer.Npiv, indexer.Npiv),
									sub_C (A.rowdim () - indexer.Npiv, indexer.Npiv);

	HybridBlocMatrix<Ring::Element, uint16> sub_B (indexer.Npiv, A.coldim () - indexer.Npiv),
									sub_D (A.rowdim () - indexer.Npiv, A.coldim () - indexer.Npiv);


	commentator.start("[Indexer] constructing sub matrices");
		indexer.constructSubMatrices(A, sub_A, sub_B, sub_C, sub_D, true);
	commentator.stop(MSG_DONE);
	MatrixUtil::show_mem_usage("[Indexer] constructing sub matrices");

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

