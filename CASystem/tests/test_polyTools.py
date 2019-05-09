'''
Created on 30.03.2019

@author: Leonard
'''
import unittest
import AlgebraicStructures as AS
import polyTools,Symbol


class TestPolyTools(unittest.TestCase):


    def test_rootsOverFiniteFields(self):
        Z7 = AS.IntegerResidueClassField(7)
        Z7x = AS.PolynomialDomain(Z7,Symbol.Symbol('x'))
        poly = Z7x([1,0,1])
        self.assertEqual(polyTools.rootsOverFiniteField(poly), [])
        poly = Z7x([-1,0,1])
        self.assertEqual(polyTools.rootsOverFiniteField(poly), [[Z7(6),1],[Z7(1),1]])
        
        Z5 = AS.IntegerResidueClassRing(5)
        Z5x,x = AS.polyRing(Z5,'x')
        f = x*(x+4)**3*(x**4+4)**4
        #print(polyTools.SFF(f))
        
        f = 3+2*x+3*x**2+3*x**3+3*x**4+3*x**5
        print(polyTools.getBerlekampMatrix(f,5))
        print(polyTools.BerlekampFactorization(f/3))
        
        Z3x,x = AS.polyRing(AS.IntegerResidueClassField(3),'x')
        
        f = x**11+2*x**9+2*x**8+x**6+x**5+2*x**3+2*x**2+1
        #print(polyTools.SFF(f))

if __name__ == "__main__":
    unittest.main()