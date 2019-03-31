'''
Created on 28.03.2019

@author: Leonard
'''

import BaseAlgebraicStructures as BAS

class Integers(BAS.EuclideanDomain):
    __instance = None
    def __new__(cls):
        if Integers.__instance is None:
            Integers.__instance = object.__new__(cls)
            super(Integers,Integers.__instance).__init__()
            Integers.__instance._zero = 0
            Integers.__instance._one = 1
        return Integers.__instance
        
    def getQuotientField(self):
        if self._quotientField!=None:
            return self._quotientField
        return Rationals()
    
    # ring methods
    def add(self, a, b):
        return a+b
    def mul(self, a, b):
        return a*b
        
    def addInverse(self, a):
        return -a
    
    def mulInverse(self, a):
        if self.isUnit(a):
            return a
        return super(PolynomialDomain,self).mulInverse(a)
    
    # euclidean methods
    def euclidFunction(self, a):
        return abs(a)
    def divisionWithRemainder(self, a, b):
        return (a//b,a%b)
    def getAssociateRepresentant(self, a):
        return abs(a),abs(a)/a
    
    def elementFromValue(self, value):
        return value
    
    def getElementType(self):
        return int
    def isUnit(self, a):
        return a in [-1,1]

class Rationals(BAS.QuotientField):
    __instance = None
    def __new__(cls):
        if Rationals.__instance==None:
            Rationals.__instance = object.__new__(cls)
            super(Rationals,Rationals.__instance).__init__(Integers())
        return Rationals.__instance
    def __init__(self):
        pass
    
    def getElementType(self):
        return DomainElement.QuotientFieldElementRational

class IntegerIdeal(BAS.Ideal):
    def __init__(self, gen):
        super(IntegerIdeal,self).__init__(Integers(),[gen])
        self.gen = gen
        
    def isMaximal(self):
        return NTU.isPrime(self.gen)
    
    
class IntegerResidueClassRing(BAS.QuotientRing):
    def __new__(cls, order, _makeFieldOpt=True):
        if type(order)!=int or order<2:
            raise Exception()
        if NTU.isPrime(order) and _makeFieldOpt:
            obj = IntegerResidueClassField(order)
            return obj
        
        obj = object.__new__(cls)
        obj._init(order)
        super(IntegerResidueClassRing,obj).__init__(Integers(),IntegerIdeal(order))
        return obj
    
    def getSimpleRepresentant(self,val):
        return val%self.order
    
    def getElementType(self):
        return DomainElement.IntegerResidueEquivalenceClass
    def __init__(self, *args):
        pass
    def _init(self, order):
        self.order = order
        self.__mulInverses={}
        
class IntegerResidueClassField(IntegerResidueClassRing, BAS.Field):
    def __new__(cls, order):
        if not NTU.isPrime(order):
            raise Exception("no field")
        obj = object.__new__(cls) 
        super(IntegerResidueClassField, obj)._init(order)
        super(IntegerResidueClassRing,obj).__init__(Integers(),IntegerIdeal(order))
        return obj
    def __init__(self, *args):
        pass
    def isUnit(self, a):
        return True if a!=self.zero else False
    def isField(self):
        return True


class SquareMatricesRing(BAS.Ring):
    def __init__(self, size, domain):
        self.size = size
        self.basedomain = domain
        self._zero = DomainElement.SquareMatrix(size, domain)
        self._one = Matrix.getUnitMatrix(size, domain)
        self.__mulInverses={}
    
    def add(self, a, b):
        return a+b
    def mul(self, a,b):
        return a*b
    def addInverse(self, a):
        return -a
    def mulInverse(self, a):
        if not self.isUnit(a):
            raise Exception("no unit")
        if a in self.__mulInverses:
            return self.__mulInverses[a]
        inverse = a.invert()
        self.__mulInverses[a] = inverse
        return inverse
    
    def isUnit(self, a):
        if a in self.__mulInverses:
            return True
        return a.hasFullRank()
        #return self.basedomain.isUnit(a.determinant())
    
    def getElementType(self):
        return DomainElement.SquareMatrix
    
    def __eq__(self, other):
        if type(other)!=SquareMatricesRing:
            return False
        return self.size==other.size and self.basedomain==other.basedomain

class PolynomialDomain(BAS.EuclideanDomain):
    """
    univariate polynomial domain
    """
    def __init__(self, domain, symbol):
        super(PolynomialDomain,self).__init__()
        self.basedomain = domain
        self.symbol = symbol
        
        self._zero = self([self.basedomain.zero])
        self._one = self([self.basedomain.one])
        
    def getQuotientField(self):
        if self._quotientField!=None:
            return self._quotientField
        return RationalFunctionsDomain(self)
        
    def add(self, a,b):
        coeffs = []
        for i in range(max(a.degree,b.degree)+1):
            coeffs.append(a[i]+b[i])
        return self(coeffs)#DomainElement.Polynomial(self, self.basedomain, coeffs)
    def mul(self, a,b):
        coeffs = []
        for i in range(a.degree+b.degree+1):
            c = sum((a[k]*b[i-k] for k in range(i+1)), self.basedomain.zero)
            coeffs.append(c)
        return self(coeffs)
    def addInverse(self, a):
        coeffs = []
        for i in range(a.degree+1):
            coeffs.append(self.basedomain.addInverse(a[i]))
        return DomainElement.Polynomial(self, self.basedomain, coeffs)
    def mulInverse(self, a):
        if self.isUnit(a):
            return self([self.basedomain.mulInverse(a[0])])
        return super(PolynomialDomain,self).mulInverse(a)
    
    def __eq__(self, other):
        if type(other)!=PolynomialDomain:
            return False
        return self.basedomain==other.basedomain and self.symbol==other.symbol
        
    def isUnit(self, a):
        return a.degree==0 and self.basedomain.isUnit(a[0])
    
    #euclid stuff
    def euclidFunction(self, a):
        return a.degree
    def divisionWithRemainder(self, a, b):
        degA,degB = a.degree, b.degree
        if degA<degB:
            return (self.zero,a)
        power = degA-degB
        coeff = a[degA]/b[degB]
        mon = self([self.basedomain.zero]*power+[coeff])
        newA = a-mon*b
        if newA==self.zero:
            return (mon,self.zero)
        else:
            (q,r) = self.divisionWithRemainder(newA, b)
            return (mon+q,r)
        
    def getAssociateRepresentant(self, a):
        if a==self.zero:
            return self.zero
        f = self.mulInverse(self([a.lcoeff]))
        return a*f,f
    
    
    def getElementType(self):
        return DomainElement.Polynomial
    
    def elementFromValue(self, value):
        if type(value)==self.getElementType():
            return value
        return self.getElementType().fromValue(value, self, self.basedomain)
    
class RationalFunctionsDomain(BAS.QuotientField):
    def __init__(self, basering):
        super(RationalFunctionsDomain, self).__init__(basering)
        if not isinstance(basering, PolynomialDomain):
            raise Exception()
        basering._rationalsDomain = self
        
        
'''class DomainAdjunctionExtension(object):
    """
    Field Extensions of type F(x) where F is a ring/field, x is a new symbol and either transcendental or algebraic over F
    """
    def __init__(self, baseField, variable, algebraic=False,minimalPolynomial=None):
        self.baseField = baseField
        self.variable = variable
        self.algebraic = algebraic
        self.minimalPolynomial = minimalPolynomial'''
        
        
import DomainElement,Matrix
import NTheoryUtils as NTU
#F_Q = Rationals()
#F_Q._setupZeroOne()
Int = Integers()
F_Q = Rationals()

FIELDS = [F_Q]