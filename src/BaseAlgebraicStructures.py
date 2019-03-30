'''
Created on 19.01.2019

@author: Leonard
'''


class Ring(object):
    _zero=None
    _one=None
    __mulInverses=None
    @property
    def zero(self):
        if self._zero==None:
            raise NotImplementedError()
        return self._zero
    @property
    def one(self):
        if self._one==None:
            raise NotImplementedError()
        return self._one
    
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
        return a+self.addInverse(b)
    def div(self, a,b):
        return a*self.mulInverse(b)
    
    def getElementType(self):
        raise NotImplementedError()
    def __call__(self, value):
        return self.elementFromValue(value)
    def elementFromValue(self, value):
        if type(value)==self.getElementType():
            return value
        return self.getElementType().fromValue(value, self)
    def isField(self):
        return False
    def isIntegralDomain(self):
        return False
    def isEuclideanDomain(self):
        return False
    
    def __ne__(self,other):
        return not self==other
    
class IntegralDomain(Ring):
    def __init__(self):
        self.allowCalcInQuotientField = True
        self._quotientField = None
        
    def getQuotientField(self):
        raise NotImplementedError()
    
        
    def mulInverse(self, a):
        if self.allowCalcInQuotientField:
            return self.getQuotientField()(self.one, a)
        raise Exception("No unit")
    def isIntegralDomain(self):
        return True
        
class EuclideanDomain(IntegralDomain):
    def __init__(self):
        super(EuclideanDomain, self).__init__()
        
    def euclidFunction(self,a):
        """
        function f:R\{0}->N0 such there is a division with remainder and f(a)<=f(ab)
        """
        raise NotImplementedError()
    def divisionWithRemainder(self,a,b):
        """
        b!=0
        returns (s,r) with a=bs+r and f(r)<f(b)
        """
        raise NotImplementedError()
    def getAssociateRepresentant(self, a):
        """
        a~b :<=> a|b and b|a defines an equivalence relation on R. this function returns a 
        "simple" representant for the equivalence class [a]
        """
        raise NotImplementedError()
    def gcd(self, a, b):
        f = self.euclidFunction
        (A,B) = (a,b) if f(a)>=f(b) else (b,a)
        if B==self.zero:
            return self.getAssociateRepresentant(A)
        (_,r) = self.divisionWithRemainder(A, B)
        if r==self.zero:
            return self.getAssociateRepresentant(B)
        return self.gcd(B,r)
    
    def isEuclideanDomain(self):
        return True
    
class Field(EuclideanDomain):
    def isUnit(self, a):
        return a!=self.zero
    def euclidFunction(self, a):
        return 0 if a==self.zero else 1
    def divisionWithRemainder(self, a, b):
        return (a/b,0)
    def getAssociateRepresentant(self, a):
        return self.one
    def isField(self):
        return True
    
import DomainElement
class QuotientField(Field):
    def __init__(self, basering):
        self.basering = basering
        self._zero = self.elementFromValue((basering.zero,basering.one),simplify=False)
        self._one = self.elementFromValue((basering.one, basering.one),simplify=False)
        self.basering._quotientField = self
        
    def add(self, a, b):
        return self(a.a*b.b+b.a*a.b,a.b*b.b)
    def mul(self, a,b):
        return self(a.a*b.a,a.b*b.b)
    def addInverse(self, a):
        return self(-a.a,a.b)
    def mulInverse(self, a):
        return self(a.b,a.a)
    
    def equal(self, a,b):
        return a.a*b.b==a.b*b.a
    
    def _simplify(self,a,b):
        if a==self.basering.zero:
            return (a,self.basering.one)
        gcd = self.basering.gcd(a,b)
        na = a//gcd
        nb = b//gcd
        nnb = self.basering.getAssociateRepresentant(nb)
        r = nnb//nb
        nna = r*na
        return (nna,nnb)
    def simplify(self, a):
        if not isinstance(self.basering, EuclideanDomain):
            raise Exception("cant simplify in non-euclidean domains")
        nna,nnb = self._simplify(a.a, a.b)
        return self(nna,nnb)
        
        
    def getElementType(self):
        return DomainElement.QuotientFieldElement
    
    def __call__(self, *value):
        if type(value[0])==list or type(value[0])==tuple:
            return self.elementFromValue(value[0])
        if len(value)==1:
            return self.elementFromValue([value[0],self.basering.one])
        return self.elementFromValue(value)
    def elementFromValue(self, value,simplify=True):
        return self.getElementType().fromValue(value, self,simplify=simplify)
    
