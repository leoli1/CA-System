'''
Created on 27.03.2019

@author: Leonard
'''
import unittest
import AlgebraicStructures as AS
import BaseAlgebraicStructures as BAS
import Symbol

class TestAlgebraicStructures(unittest.TestCase):
    def test_IntegerResidueClass(self):
        Z5 = AS.IntegerResidueClassRing(5)
        a = Z5(4)
        b = Z5(3)
        print(a/b)
    def test_Rationals(self):
        Q = AS.F_Q
        self.assertIsInstance(Q, BAS.QuotientField)
        a1 = Q(12,3)
        a2 = Q.elementFromValue((12,3),simplify=False)
        b1 = Q(5,-3)
        b2 = Q.elementFromValue((5,-3),simplify=False)
        self.assertEqual(str(a2), "12/3")
        self.assertEqual(a1, a2)
        self.assertEqual(a1+b1, b2+a2)
        self.assertEqual(a1+b1,Q(7,3))
    def test_ResidueClass(self):
        r1 = AS.IntegerResidueClassRing(5)
        self.assertIsInstance(r1, AS.IntegerResidueClassField)
        r2 = AS.IntegerResidueClassRing(6)
        self.assertNotIsInstance(r2, AS.IntegerResidueClassField)
        
        a = r1.elementFromValue(3)
        b = r1.elementFromValue(4)
        c = r1.mulInverse(a)
        d = r1.elementFromValue(2)
        self.assertEqual(c, d)
        self.assertEqual(a+b, c)
        
    def test_PolynomialDomain(self):
        domain = AS.IntegerResidueClassField(7)
        polydomain = AS.PolynomialDomain(domain, Symbol.Symbol('x'))
        poly1 = polydomain.elementFromValue([1,1])
        poly2 = polydomain.elementFromValue([0,2,3,0,7])
        poly3 = polydomain.elementFromValue([2])
        self.assertEqual(poly1.degree, 1)
        self.assertEqual(poly2.degree, 2)
        self.assertEqual(str(poly2), "3x^2+2x")
        self.assertEqual(str(poly1+poly2), "3x^2+3x+1")
        self.assertFalse(poly1.isUnit())
        self.assertTrue(poly3.isUnit())
        self.assertEqual(poly3.inverse(), polydomain([4]))
        p3 = poly2*poly1
        p4 = polydomain([0,2,5,3])
        self.assertEqual(p3, p4)
        self.assertNotEqual(poly2*poly3, poly2)
        
        Q = AS.F_Q
        polydomain = AS.PolynomialDomain(Q,Symbol.Symbol('x'))
        p1 = polydomain([(3,2),1,(4,4)])
        p2 = polydomain([(9,4),3,4,2,1])
        self.assertEqual(str(p1), "x^2+x+3/2")
        self.assertEqual(p1**2, p2)
        
        p1 = polydomain([2,2,1])
        p2 = polydomain([1,2])
        self.assertEqual(str(polydomain.divisionWithRemainder(p1, p2)), "(1/2x+3/4, 5/4)")
        self.assertEqual(p1.evaluate(Q(2)), Q(10))
        
        rational = polydomain.getQuotientField()(p1,p2)
        self.assertEqual(p1/p2, rational)
        
        

if __name__ == "__main__":
    unittest.main()