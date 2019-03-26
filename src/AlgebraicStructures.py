'''
Created on 19.01.2019

@author: Leonard
'''


class Ring(object):
    zero=None
    one=None
    
    def getZero(self):
        if self.zero==None:
            raise NotImplementedError()
        return self.zero
    def getOne(self):
        if self.one==None:
            raise NotImplementedError()
        return self.one
    
    def add(self, a, b):
        raise NotImplementedError()
    def mul(self, a, b):
        raise NotImplementedError()
    def addInverse(self, a):
        raise NotImplementedError()
    def mulInverse(self,a):
        raise NotImplementedError()
    
    def isUnit(self, a):
        return NotImplementedError()
    
    def sub(self, a,b):
        return self.add(a,self.addInverse(b))
    def div(self, a,b):
        return self.mul(a,self.mulInverse(b))
    
    def getElementType(self):
        raise NotImplementedError()
    def elementFromValue(self, value):
        if type(value)==self.getElementType():
            return value
        return self.getElementType().fromValue(value)
    
class Field(Ring):
    def isUnit(self, a):
        return a!=self.getZero()
    

class ResidueClassRing(Ring):
    def __init__(self, order):
        if type(order)!=int or order<2:
            raise Exception()
        self.order = order
        self.__mulInverses={}
        self.zero = Number.ModularInteger(0,self)
        self.one = Number.ModularInteger(1,self)
        
    def add(self, a, b):
        return Number.ModularInteger(a.val+b.val, self)
    def mul(self, a, b):
        return Number.ModularInteger(a.val*b.val, self)
    def addInverse(self, a):
        return Number.ModularInteger(self.order-a-1, self)
    def mulInverse(self, a):
        if not self.isUnit(a):
            raise Exception("cant invert non-unit")
        ## euclidean alg
        if a in self.__mulInverses:
            return self.__mulInverses[a]
        else:
            return Number.ModularInteger(NTU.euclid(a,self.order)[0], self)

    def isUnit(self, a):
        return NTU.gcd(a.val,self.order)==1
    def getElementType(self):
        return Number.ModularInteger
    def elementFromValue(self, value):
        if type(value)==self.getElementType():
            return value
        return self.getElementType().fromValue(value, self)
    
class ResidueClassField(ResidueClassRing, Field):
    def __init__(self, order):
        if not NTU.isPrime(order):
            raise Exception("no field")
        super(ResidueClassField, self).__init__(order)
        
class Rationals(Field):

    def __init__(self):
        pass
    def _setupZeroOne(self):
        self.zero = Number.RationalNumber(0,1)
        self.one = Number.RationalNumber(1,1)
    def add(self, a, b):
        #a,b = self._makeRat(a),self._makeRat(b)
        return Number.RationalNumber(a.num*b.denom+a.denom*b.num,a.denom*b.denom)
    def mul(self,a,b):
        #a,b = self._makeRat(a),self._makeRat(b)
        return Number.RationalNumber(a.num*b.num,a.denom*b.denom)
    def addInverse(self, a):
        #a = self._makeRat(a)
        return Number.RationalNumber(-a.num,a.denom)
    def mulInverse(self, a):
        #a = self._makeRat(a)
        return Number.RationalNumber(a.denom,a.num)
    
    def getElementType(self):
        return Number.RationalNumber
    
class PolynomialsOverField(Ring):
    def __init__(self, baseField):
        self.baseField = baseField
        
class FieldAdjunctionExtension(object):
    """
    Field Extensions of type F(x) where F is a field, x is a new symbol and either transcendental or algebraic over F
    """
    def __init__(self, baseField, variable, algebraic=False,minimalPolynomial=None):
        self.baseField = baseField
        self.variable = variable
        self.algebraic = algebraic
        self.minimalPolynomial = minimalPolynomial
        
        
import Number
import NTheoryUtils as NTU
F_Q = Rationals()
F_Q._setupZeroOne()

FIELDS = [F_Q]