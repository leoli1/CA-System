'''
Created on 30.03.2019

@author: Leonard
'''
import unittest
import AlgebraicStructures as AS
import polyTools,Symbol


class TestPolyTools(unittest.TestCase):


    def test_rootsOverFiniteFields(self):
        Z7 = AS.ResidueClassField(7)
        Z7x = AS.PolynomialDomain(Z7,Symbol.Symbol('x'))
        poly = Z7x([1,0,1])
        self.assertEqual(polyTools.rootsOverFiniteField(poly), None)
        poly = Z7x([-1,0,1])
        self.assertEqual(polyTools.rootsOverFiniteField(poly), [[Z7(6),1],[Z7(1),1]])


if __name__ == "__main__":
    unittest.main()