#include <linbox/field/modular.h>
#include <linbox/blackbox/dense.h>
#include <linbox/solutions/det.h>

using namespace LinBox;

main()
{
    typedef Modular<double> Field;
    typedef DenseMatrix<Field> Matrix;
    Field F(65521);
    Matrix A; A.read("datafile");
    Field::Element d;
    det(d, A, blasElimination);
    F.write(cout << "the determinant is ", d) << endl;
}

