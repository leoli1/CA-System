'''
Created on 28.03.2019

@author: Leonard
'''
import unittest
import AlgebraicStructures as AS
import Matrix,Symbol


class TestMatrix(unittest.TestCase):


    def test_Matrix(self):
        Z5 = AS.IntegerResidueClassField(5)
        m = Matrix.Matrix(2,2,Z5,[[2,1],[0,2]])
        n = Matrix.Matrix(2,2,Z5,[[3,1],[0,3]])
        self.assertEqual(str(m), "[2|1]\n[0|2]")
        self.assertEqual(m.invert(),n)
        polys = AS.PolynomialDomain(Z5,Symbol.Symbol('t'))
        self.assertEqual(m.charPoly(), polys([4,1,1]))
        self.assertEqual(m.determinant(), Z5.elementFromValue(4))
        m2 = Matrix.Matrix(2,2,Z5,[[1,3],[2,6]])
        with self.assertRaises(Matrix.MatrixInversionError):
            m2.invert()
        self.assertEqual(m2.determinant(), Z5(0))
        
        k = Matrix.getUnitMatrix(2, Z5)
        self.assertTrue(k.isDiagonalizable())
        (p,d) = k.getDiagonalizerMatrix()
        self.assertEqual(p, d)
        self.assertEqual(p, k)
        
        m = Matrix.Matrix(2,2,Z5,[[1,1],[1,1]])
        self.assertTrue(m.isDiagonalizable())
        (p,d) = m.getDiagonalizerMatrix()
        self.assertEqual(p,Matrix.Matrix(2,2,Z5,[[4,1],[1,1]]))
        
        Q = AS.F_Q
        m = Matrix.Matrix(3,3,Q,[[1,1,0],[2,0,2],[-1,1,1]])
        n = Matrix.Matrix(3,3,Q,[[2,1,-2],[4,-1,2],[-2,2,2]]).scalarMultiplication(Q(1,6))
        self.assertEqual(n, m.invert())
        self.assertEqual(n.invert(), m)
        
        self.assertEqual(m.determinant(),Q(-6))
        polys = AS.PolynomialDomain(Q,Symbol.Symbol("t"))
        self.assertEqual(m.charPoly(), polys([6,-3,-2,1]))
        
        
        m = Matrix.Matrix(3,3,Q,[[1,1,0],[0,0,1],[0,0,0]])
        self.assertEqual(m.nullity(), 1)
        
        
        


if __name__ == "__main__":
    unittest.main()