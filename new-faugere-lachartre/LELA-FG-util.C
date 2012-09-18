/*
 * LELA-FG-util.C
 * Copyright 2012 Martani Fayssal (LIP6 / UPMC University Paris06)
 *
 *  Created on: 25 mai 2012
 *      Author: martani (LIP6 / UPMC University Paris06)
 * 
 * This exposes the exact functionality in ../lela/algorithms/faugere-lachartre
 * The algorithm is broken into parts to return the intermediate computation results
 */

#ifndef _LELA_FG_util_
#define _LELA_FG_util_

#include "lela/util/commentator.h"
#include "lela/algorithms/faugere-lachartre.h"
#include "lela/matrix/sparse.h"
#include "lela/matrix/dense.h"

using namespace LELA;


namespace LELA_GF_UTIL
{
	template <class Ring, class Matrix1, class Matrix2, class Matrix3>
	class MatrixGrid1
	{
		const Ring &R;
		Matrix1 &X;
		Matrix2 &A;
		Matrix3 &B, &C, &D;

		template <class Container>
		void moveRowSpecialised (const Block &horiz_block, int row, const Container &vert_blocks, VectorRepresentationTypes::Generic)
		{
			typename Container::const_iterator vert_block;
			typename Matrix1::ConstRowIterator v_X = X.rowBegin () + (horiz_block.sourceIndex () + row);
			typename Matrix1::ConstRow::const_iterator i = v_X->begin ();

			if (horiz_block.dest () == 0) {
				typename Matrix2::RowIterator v_A = A.rowBegin () + (horiz_block.destIndex () + row);
				typename Matrix3::RowIterator v_B = B.rowBegin () + (horiz_block.destIndex () + row);

				for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
					if (vert_block->dest () == 0)
						i = Splicer::moveBlock (R, *v_A, i, v_X->end (), vert_block->sourceIndex (), vert_block->destIndex (), vert_block->size (),
									typename VectorTraits<Ring, typename Matrix1::ConstRow>::RepresentationType ());
					else
						i = Splicer::moveBlock (R, *v_B, i, v_X->end (), vert_block->sourceIndex (), vert_block->destIndex (), vert_block->size (),
									typename VectorTraits<Ring, typename Matrix1::ConstRow>::RepresentationType ());
				}
			} else {
				typename Matrix3::RowIterator v_C = C.rowBegin () + (horiz_block.destIndex () + row);
				typename Matrix3::RowIterator v_D = D.rowBegin () + (horiz_block.destIndex () + row);

				for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
					if (vert_block->dest () == 0)
						i = Splicer::moveBlock (R, *v_C, i, v_X->end (), vert_block->sourceIndex (), vert_block->destIndex (), vert_block->size (),
									typename VectorTraits<Ring, typename Matrix1::ConstRow>::RepresentationType ());
					else
						i = Splicer::moveBlock (R, *v_D, i, v_X->end (), vert_block->sourceIndex (), vert_block->destIndex (), vert_block->size (),
									typename VectorTraits<Ring, typename Matrix1::ConstRow>::RepresentationType ());
				}
			}
		}

		template <class Container>
		void moveRowSpecialised (const Block &horiz_block, int row, const Container &vert_blocks, VectorRepresentationTypes::Dense01)
		{
			typename Container::const_iterator vert_block;
			typename Matrix1::ConstRowIterator v_X = X.rowBegin () + (horiz_block.sourceIndex () + row);

			if (horiz_block.dest () == 0) {
				typename Matrix2::RowIterator v_A = A.rowBegin () + (horiz_block.destIndex () + row);
				typename Matrix3::RowIterator v_B = B.rowBegin () + (horiz_block.destIndex () + row);

				for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
					typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType v (*v_X, vert_block->sourceIndex (), vert_block->sourceIndex () + vert_block->size ());

					if (vert_block->dest () == 0)
						Splicer::moveBitBlockDense (R, *v_A, v, vert_block->sourceIndex (), vert_block->destIndex ());
					else
						Splicer::moveBitBlockDense (R, *v_B, v, vert_block->sourceIndex (), vert_block->destIndex ());
				}
			} else {
				typename Matrix3::RowIterator v_C = C.rowBegin () + (horiz_block.destIndex () + row);
				typename Matrix3::RowIterator v_D = D.rowBegin () + (horiz_block.destIndex () + row);

				for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
					typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType v (*v_X, vert_block->sourceIndex (), vert_block->sourceIndex () + vert_block->size ());

					if (vert_block->dest () == 0)
						Splicer::moveBitBlockDense (R, *v_C, v, vert_block->sourceIndex (), vert_block->destIndex ());
					else
						Splicer::moveBitBlockDense (R, *v_D, v, vert_block->sourceIndex (), vert_block->destIndex ());
				}
			}
		}

		template <class Container>
		void moveRowSpecialised (const Block &horiz_block, int row, const Container &vert_blocks, VectorRepresentationTypes::Hybrid01)
		{
			typename Container::const_iterator vert_block;
			typename Matrix1::ConstRowIterator v_X = X.rowBegin () + (horiz_block.sourceIndex () + row);
			typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType::const_iterator i;

			typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType vp (*v_X, 0, vert_blocks.back ().sourceIndex () + vert_blocks.back ().size ());
			i = vp.begin ();

			if (horiz_block.dest () == 0) {
				typename Matrix2::RowIterator v_A = A.rowBegin () + (horiz_block.destIndex () + row);
				typename Matrix3::RowIterator v_B = B.rowBegin () + (horiz_block.destIndex () + row);

				for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
					typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType w (*v_X, i, vert_block->sourceIndex (), vert_block->sourceIndex () + vert_block->size ());

					if (vert_block->dest () == 0)
						i = Splicer::moveBitBlockHybrid (R, *v_A, w, vert_block->destIndex ());
					else
						i = Splicer::moveBitBlockHybrid (R, *v_B, w, vert_block->destIndex ());
				}
			} else {
				typename Matrix3::RowIterator v_C = C.rowBegin () + (horiz_block.destIndex () + row);
				typename Matrix3::RowIterator v_D = D.rowBegin () + (horiz_block.destIndex () + row);

				for (vert_block = vert_blocks.begin (); vert_block != vert_blocks.end (); ++vert_block) {
					typename VectorTraits<Ring, typename Matrix1::ConstRow>::ConstSubvectorType w (*v_X, i, vert_block->sourceIndex (), vert_block->sourceIndex () + vert_block->size ());

					if (vert_block->dest () == 0)
						i = Splicer::moveBitBlockHybrid (R, *v_C, w, vert_block->destIndex ());
					else
						i = Splicer::moveBitBlockHybrid (R, *v_D, w, vert_block->destIndex ());
				}
			}
		}

	public:
		typedef GridTypeRowOptimised GridType;

		MatrixGrid1 (const Ring &__R, Matrix1 &__X, Matrix2 &__A, Matrix3 &__B, Matrix3 &__C, Matrix3 &__D)
			: R (__R), X (__X), A (__A), B (__B), C (__C), D (__D)
			{}

		template <class Container>
		void moveRow (const Block &horiz_block, int row, const Container &vert_blocks)
			{ moveRowSpecialised (horiz_block, row, vert_blocks, typename VectorTraits<Ring, typename Matrix1::Row>::RepresentationType ()); }
	};

	template <class Ring, class Matrix>
	class MatrixGrid2
	{
		const Ring &R;
		Matrix &B, &B1, &B2;

	public:
		typedef GridTypeNormal GridType;

		MatrixGrid2 (const Ring &__R, Matrix &__B, Matrix &__B1, Matrix &__B2)
			: R (__R), B (__B), B1 (__B1), B2 (__B2)
			{}

		void operator () (const Block &horiz_block, const Block &vert_block)
		{
			if (horiz_block.dest () == 0) {
				if (vert_block.dest () == 0)
					Splicer::copyBlock (R, B, B1, horiz_block, vert_block);
				else
					Splicer::copyBlock (R, B, B2, horiz_block, vert_block);
			}
		}
	};

	template <typename Ring, typename Matrix>
	void setup_splicer (Ring& R, Splicer &splicer, Splicer &reconst_splicer, const Matrix &A, size_t &num_pivot_rows, typename Ring::Element &det)
	{
		commentator.start ("Finding pivot-rows", __FUNCTION__);

		typename Matrix::ConstRowIterator i_A;
		int last_col = -1, col, first_col_in_block = 0, height = 0, dest_col_tr = 0, dest_col_res = 0;
		int row = 0, first_row_in_block = 0, dest_row_tr = 0, dest_row_res = 0;
		bool last_was_same_col = false;
		typename Ring::Element a;
		Context<Ring> ctx (R);

		num_pivot_rows = 0;

		for (i_A = A.rowBegin (); i_A != A.rowEnd (); ++i_A, ++row) {
			col = BLAS1::head (ctx, a, *i_A);

			if (col == -1) {
				if (!last_was_same_col) {
					splicer.addHorizontalBlock (Block (0, 0, first_row_in_block, dest_row_tr, row - first_row_in_block));
					dest_row_tr += row - first_row_in_block;
					first_row_in_block = row;
				}
				break;
			}
			else if (col == last_col) {
				if (!last_was_same_col) {
					splicer.addHorizontalBlock (Block (0, 0, first_row_in_block, dest_row_tr, row - first_row_in_block));
					dest_row_tr += row - first_row_in_block;
					first_row_in_block = row;
					last_was_same_col = true;
				}
			}
			else if (col < last_col)
				throw "WrongMatrixForm (row)";
			else {
				if (last_was_same_col) {
					splicer.addHorizontalBlock (Block (0, 1, first_row_in_block, dest_row_res, row - first_row_in_block));
					dest_row_res += row - first_row_in_block;
					first_row_in_block = row;
					last_was_same_col = false;
				}

				if (col == last_col + 1) {
					last_col = col;
					++height;
					++num_pivot_rows;
				} else {
					if (height > 0) {
						Block newblock (0, 0, dest_col_tr, first_col_in_block, height);

						splicer.addVerticalBlock (Block (0, 0, first_col_in_block, dest_col_tr, height));
						reconst_splicer.addVerticalBlock (newblock);
						reconst_splicer.addHorizontalBlock (newblock);
					}

					if (col - last_col > 1) {
						Block newblock (1, 0, dest_col_res, first_col_in_block + height, col - last_col - 1);

						splicer.addVerticalBlock (Block (0, 1, first_col_in_block + height, dest_col_res, col - last_col - 1));
						reconst_splicer.addVerticalBlock (newblock);
						reconst_splicer.addHorizontalBlock (newblock);
					}

					dest_col_tr += height;
					dest_col_res += col - last_col - 1;
					height = 1;
					first_col_in_block = last_col = col;
					++num_pivot_rows;
				}
			}

			ctx.F.mulin (det, a);
		}

		if (height > 0) {
			Block newblock (0, 0, dest_col_tr, first_col_in_block, height);

			splicer.addVerticalBlock (Block (0, 0, first_col_in_block, dest_col_tr, height));
			reconst_splicer.addVerticalBlock (newblock);
			reconst_splicer.addHorizontalBlock (newblock);
		}

		if (first_col_in_block + height < (int) A.coldim ()) {
			Block newblock (1, 0, dest_col_res, first_col_in_block + height, A.coldim () - first_col_in_block - height);

			splicer.addVerticalBlock (Block (0, 1, first_col_in_block + height, dest_col_res, A.coldim () - first_col_in_block - height));
			reconst_splicer.addVerticalBlock (newblock);
			reconst_splicer.addHorizontalBlock (newblock);
		}

		if (first_row_in_block < (int) A.rowdim ()) {
			if (row < (int) A.rowdim () || last_was_same_col)
				splicer.addHorizontalBlock (Block (0, 1, first_row_in_block, dest_row_res, A.rowdim () - first_row_in_block));
			else
				splicer.addHorizontalBlock (Block (0, 0, first_row_in_block, dest_row_tr, A.rowdim () - first_row_in_block));
		}

		commentator.stop (MSG_DONE, NULL, __FUNCTION__);
	}

	template <typename Matrix, typename Ring>
	static bool checkSubMatricesOfM0SameAsLELASplicer(Ring& R, const Matrix& X, const Matrix& sub_A, const Matrix& sub_B, const Matrix& sub_C, const Matrix& sub_D)
	{
		commentator.start ("LELA-FG-util", __FUNCTION__);

		bool pass = true;
		typedef SparseMatrix<typename Ring::Element> _Sparse;
		typedef DenseMatrix<typename Ring::Element> _Dense;
		Context<Ring> ctx (R);

		std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_WARNING);

		Splicer X_splicer, X_reconst_splicer;

		size_t num_pivot_rows;
		typename Ring::Element det;
		R.copy (det, R.one ());

		setup_splicer (R, X_splicer, X_reconst_splicer, X, num_pivot_rows, det);
		//rank = num_pivot_rows;

		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Found " << num_pivot_rows << " pivots" << std::endl;
		reportUI << "Splicer:" << std::endl << X_splicer << std::endl;

		_Sparse A (num_pivot_rows, num_pivot_rows);
		_Dense B (num_pivot_rows, X.coldim () - num_pivot_rows);
		_Dense C (X.rowdim () - num_pivot_rows, num_pivot_rows);
		_Dense D (X.rowdim () - num_pivot_rows, X.coldim () - num_pivot_rows);

		X_splicer.splice (MatrixGrid1<Ring, const Matrix, _Sparse, _Dense > (R, X, A, B, C, D));

		std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING);

		if(BLAS3::equal(ctx, sub_A, A))
			report << "<<A CORRECT>>\t";
		else
		{
			report << "<<<<<<<<<<A NOT CORRECT>>>>>>>>" << std::endl;
			pass = false;
		}

		if(BLAS3::equal(ctx, sub_B, B))
			report << "<<B CORRECT>>\t";
		else
		{
			report << "<<<<<<<<<<B NOT CORRECT>>>>>>>>" << std::endl;
			pass = false;
		}

		if(BLAS3::equal(ctx, sub_C, C))
			report << "<<C CORRECT>>\t";
		else
		{
			report << "<<<<<<<<<<C NOT CORRECT>>>>>>>>" << std::endl;
			pass = false;
		}

		if(BLAS3::equal(ctx, sub_D, D))
			report << "<<D CORRECT>>\t";
		else
		{
			reportUI << "<<<<<<<<<<D NOT CORRECT>>>>>>>>" << std::endl;
			pass = false;
		}
		reportUI << std::endl;

		commentator.stop(pass ? MSG_PASSED : MSG_FAILED);
		return pass;
	}

	template <typename Matrix, typename Ring>
	static bool checkSubMatricesOfDSameAsLELASplicer(Ring& R, const Matrix& orig_B, const Matrix& orig_D,
											const Matrix& sub_B1, const Matrix& sub_B2,
											const Matrix& sub_D1, const Matrix& sub_D2)
	{
		commentator.start ("LELA-FG-util", __FUNCTION__);

		Splicer D_splicer, D_reconst_splicer;
		size_t num_pivot_rows;
		typename Ring::Element det;
		Context<Ring> ctx (R);

		R.copy (det, R.one ());

		DenseMatrix<typename Ring::Element> D (orig_D.rowdim (), orig_D.coldim ()),
						    B (orig_B.rowdim (), orig_B.coldim ());
		BLAS3::copy(ctx, orig_B, B);
		BLAS3::copy(ctx, orig_D, D);

		setup_splicer (R, D_splicer, D_reconst_splicer, D, num_pivot_rows, det);

		std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_WARNING);

		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "(In D) found " << num_pivot_rows << " pivots" << std::endl;
		reportUI << "Splicer:" << std::endl << D_splicer << std::endl;

		DenseMatrix<typename Ring::Element> B1 (B.rowdim (), num_pivot_rows);
		DenseMatrix<typename Ring::Element> B2 (B.rowdim (), D.coldim () - num_pivot_rows);
		DenseMatrix<typename Ring::Element> D1 (num_pivot_rows, num_pivot_rows);
		DenseMatrix<typename Ring::Element> D2 (num_pivot_rows, D.coldim () - num_pivot_rows);

		Splicer B_splicer (D_splicer);
		B_splicer.clearHorizontalBlocks ();
		B_splicer.addHorizontalBlock (Block (0, 0, 0, 0, B.rowdim ()));

		B_splicer.splice (MatrixGrid2<Ring, DenseMatrix<typename Ring::Element> > (R, B, B1, B2));
		D_splicer.splice (MatrixGrid2<Ring, DenseMatrix<typename Ring::Element> > (R, D, D1, D2));

		bool pass = true;

		std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING);

		if(BLAS3::equal(ctx, sub_B1, B1))
			report << "<<B1 CORRECT>>\t";
		else
		{
			report << "<<<<<<<<<<B1 NOT CORRECT>>>>>>>>" << std::endl;
			pass = false;
		}

		if(BLAS3::equal(ctx, sub_B2, B2))
			report << "<<B2 CORRECT>>\t";
		else
		{
			report << "<<<<<<<<<<B2 NOT CORRECT>>>>>>>>" << std::endl;
			pass = false;
		}

		if(BLAS3::equal(ctx, sub_D1, D1))
			report << "<<D1 CORRECT>>\t";
		else
		{
			report << "<<<<<<<<<<D1 NOT CORRECT>>>>>>>>" << std::endl;
			pass = false;
		}

		if(BLAS3::equal(ctx, sub_D2, D2))
			report << "<<D2 CORRECT>>\t";
		else
		{
			report << "<<<<<<<<<<D2 NOT CORRECT>>>>>>>>" << std::endl;
			pass = false;
		}
		report << std::endl;

		commentator.stop(pass ? MSG_PASSED : MSG_FAILED);
		return pass;
	}

	template <typename Matrix, typename Ring>
        static void spliceMatrix(Ring& R, const Matrix& X)
        {
                commentator.start ("LELA-FG-util - Splicer setup", __FUNCTION__);

                typedef SparseMatrix<typename Ring::Element> _Sparse;
                typedef DenseMatrix<typename Ring::Element> _Dense;
                Context<Ring> ctx (R);
                std::ostream &reportUI = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_WARNING);

                Splicer X_splicer, X_reconst_splicer;

                size_t num_pivot_rows;
                typename Ring::Element det;
                R.copy (det, R.one ());
                
                LELA_GF_UTIL::setup_splicer (R, X_splicer, X_reconst_splicer, X, num_pivot_rows, det);
                //rank = num_pivot_rows;

                commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                        << "Found " << num_pivot_rows << " pivots" << std::endl;
                reportUI << "Splicer:" << std::endl << X_splicer << std::endl;

                _Sparse A (num_pivot_rows, num_pivot_rows);
                _Dense B (num_pivot_rows, X.coldim () - num_pivot_rows);
                _Dense C (X.rowdim () - num_pivot_rows, num_pivot_rows);
                _Dense D (X.rowdim () - num_pivot_rows, X.coldim () - num_pivot_rows);

                X_splicer.splice (MatrixGrid1<Ring, const Matrix, _Sparse, _Dense > (R, X, A, B, C, D));       
                
                commentator.stop(MSG_DONE);
        }
}


#endif  //#define _LELA-FG-util_
