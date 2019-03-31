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
        "simple" representant b for the equivalence class [a] and the factor x, that got it into
        this form, so b=ax
        """
        raise NotImplementedError()
    def gcd(self, a, b):
        f = self.euclidFunction
        (A,B) = (a,b) if f(a)>=f(b) else (b,a)
        if B==self.zero:
            return self.getAssociateRepresentant(A)[0]
        (_,r) = self.divisionWithRemainder(A, B)
        if r==self.zero:
            return self.getAssociateRepresentant(B)[0]
        return self.gcd(B,r)
    def euclid(self,a,b,d=None):
        (q,r) = self.divisionWithRemainder(a, b)
        if d==None:
            d = self.gcd(a,b)
        r,f = self.getAssociateRepresentant(r)
        if r==d:
            x = self.mulInverse(f)
            y = -q*x
            return (x,y)
        anew = b
        bnew = r
        xn,yn=self.euclid(anew,bnew,d=d)
        return (yn,xn-yn*q)
    def coprime(self,a,b):
        return self.gcd(a,b)==self.getAssociateRepresentant(self.one)[0]
    def divisibleBy(self, a,b):
        return self.divisionWithRemainder(a, b)[1]==self.zero
    
    
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
        return self.one,self.mulInverse(a)
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
        nnb,f = self.basering.getAssociateRepresentant(nb)
        nna = f*na
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
    
class Ideal(Ring):
    def __init__(self, ring, generator):
        self.ring = ring
        self.generator = generator
        
        self.add = ring.add
        self.mul = ring.mul
        self.isUnit = ring.isUnit
        self.getElementType = ring.getElementType
        self._zero = ring._zero
        
    def inIdeal(self, a):
        if self.ring.isEuclideanDomain():
            if self.isPrincipal():
                if (self.generator[0]==self.zero):
                    return a==self.zero
                return self.ring.divisibleBy(a,self.generator[0])
            raise NotImplementedError()
        raise NotImplementedError()
    
    def isPrincipal(self):
        return len(self.generator)==1
    def isMaximal(self):
        return False
    
class QuotientRing(Ring):
    """
    given a ring R and an ideal I in R, a~b :<=> a-b in I defines an equivalence relation on R
    the set of equivalence classes defines a new ring - the quotient ring
    if I is maximal the quotient ring is actually a field
    """
    def __init__(self,ring, ideal):
        self.ring = ring
        self.ideal = ideal
        self._zero = self(ring.zero)
        self._one = self(ring.one)
        if ideal.ring!=ring:
            raise Exception()
        
    def add(self, a,b):
        return self(a.val+b.val)
    def mul(self, a,b):
        return self(a.val*b.val)
    def addInverse(self,a):
        return self(-a.val)
    def mulInverse(self,a):
        if self.ring.isEuclideanDomain() and self.ideal.isPrincipal():
            return self(self.ring.euclid(a.val,self.ideal.generator[0])[0])
        raise NotImplementedError()
    
    def isUnit(self, a):
        if self.ring.isEuclideanDomain() and self.ideal.isPrincipal():
            return self.ring.coprime(a.val,self.ideal.generator[0])
        if self.ideal.isMaximal():
            return a!=self.zero
        raise NotImplementedError()
    def equal(self, a,b):
        return self.ideal.inIdeal(b.val-a.val)
    
    def isField(self):
        return self.ideal.isMaximal()
    
    def getSimpleRepresentant(self,val):
        raise NotImplementedError()
    
    def getElementType(self):
        return DomainElement.QuotientRingEquivalenceClass
        