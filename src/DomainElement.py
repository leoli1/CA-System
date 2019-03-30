'''
Created on 19.01.2019

@author: Leonard
'''
from __future__ import division
import AlgebraicStructures as AS
import BaseAlgebraicStructures as BAS
import polyTools

class DomainElement(object):

    def __init__(self, domain):
        """
        domain should be a ring or a field
        """
        self.domain = domain
    
    def __add__(self, other):
        self, other = DomainElement.__makeSameDomain(self, other)
        return self.domain.add(self,other)
    def __mul__(self, other):
        self, other = DomainElement.__makeSameDomain(self, other)
        return self.domain.mul(self, other)
    def __sub__(self, other):
        self, other = DomainElement.__makeSameDomain(self, other)
        return self.domain.sub(self,other)
    def __div__(self, other):
        
        return self/other
    def __truediv__(self, other):
        self, other = DomainElement.__makeSameDomain(self, other)
        return self.domain.div(self, other)
    def __floordiv__(self, other):
        self, other = DomainElement.__makeSameDomain(self, other)
        if isinstance(self.domain, BAS.EuclideanDomain):
            return self.domain.divisionWithRemainder(self,other)[0]
    def __neg__(self):
        return self.domain.addInverse(self)
    def __pow__(self, other):
        if type(other)!=int:
            raise Exception()
        if other==0:
            return self.domain.one
        if other<0:
            return self.inverse()*self.__pow__(other+1)
        else:
            return self*self.__pow__(other-1)
        
    def inverse(self):
        return self.domain.mulInverse(self)
    def isUnit(self):
        return self.domain.isUnit(self)
    
    def __eq__(self, other):
        raise NotImplementedError()
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def toQuotientFieldElement(self):
        if not self.domain.isIntegralDomain():
            raise Exception()
        return self.domain.getQuotientField()(self,self.domain.one)
    
    @staticmethod
    def __makeSameDomain(a,b):
        if type(a)==type(b):
            return a,b
        
        if not a.domain.isIntegralDomain() or not b.domain.isIntegralDomain():
            raise Exception()
        if a.domain.isField():
            #-> self is in quotientfield, other is in basering
            b = b.toQuotientFieldElement()
        else:
            #-> self is in basering, other in quotientfield
            a = a.toQuotientFieldElement()
        return a,b
    @staticmethod
    def fromValue(val, domain):
        raise NotImplementedError()
    
class QuotientFieldElement(DomainElement):
    def __init__(self, domain, basedomain, a,b,simplify=True):
        self.domain = domain
        self.basedomain = basedomain
        self.a, self.b = a,b
        if isinstance(self.basedomain, BAS.EuclideanDomain) and simplify:
            (self.a,self.b) = self.domain._simplify(a,b)
        if type(a)!=self.basedomain.getElementType():
            raise Exception()
        
    @classmethod
    def fromValue(cls,val, domain,simplify=True):
        if type(val)==domain.getElementType():
            return val
        if type(val)==list or type(val)==tuple:
            return cls(domain, domain.basering, *val,simplify=simplify)
        return cls(domain,domain.basering,val,domain.basering.one,simplify=simplify)
        
    def __eq__(self, other):
        if id(other)==id(None):
            return False
        return self.domain.equal(self,other)
    
    def __repr__(self):
        return "({})/({})".format(self.a,self.b)
    
    def asBaseRingElement(self):
        if self.b==self.basedomain.one:
            return self.a
        raise Exception()
class QuotientFieldElementRational(QuotientFieldElement):
    def __repr__(self):
        if self.b==1:
            return str(self.a)
        return "{}/{}".format(self.a,self.b)
        
    
class ModularInteger(DomainElement):
    def __init__(self, val, domain):
        super(ModularInteger, self).__init__(domain)
        self.val = val%domain.order
        
    
    def __eq__(self, other):
        if id(other)==id(None):
            return False
        return (other.val-self.val)%self.domain.order==0
    
    @staticmethod
    def fromValue(val, domain):
        return ModularInteger(val, domain)
    
    def __repr__(self):
        return str(self.val%self.domain.order)
    
        
class Polynomial(DomainElement):#TODO
    def __init__(self, domain, basedomain, coeffs=None):
        self.basedomain = basedomain
        self.domain = domain
        
        if coeffs==None:
            self.coeffs = [domain.zero]
        else:
            self.coeffs = map(basedomain.elementFromValue, coeffs)
            self.removeZeros()
        
    def __getitem__(self, index):
        return self.coeffs[index] if index<=self.degree else self.basedomain.zero
    @property
    def degree(self):
        return len(self.coeffs)-1
    @property
    def lcoeff(self):
        return self[self.degree]
    def removeZeros(self):
        deg=0
        for i in range(len(self.coeffs)-1,0,-1):
            if self[i]!=self.basedomain.zero:
                deg=i
                break
        self.coeffs = self.coeffs[0:deg+1]
        
    def __eq__(self, other):
        if id(other)==id(None):
            return False
        if self.domain!=other.domain:
            return False
        self.removeZeros()
        other.removeZeros()
        if self.degree!=other.degree:
            return False
        for i in range(self.degree+1):
            if self[i]!=other[i]:
                return False
        return True
    
    def evaluate(self, value):
        return sum((self[i]*(value**i) for i in range(self.degree+1)), self.basedomain.zero)
    def findRoots(self):
        if isinstance(self.basedomain, AS.ResidueClassField):
            return polyTools.rootsOverFiniteField(self)
        raise NotImplementedError()
    @staticmethod
    def fromValue(val, domain, basedomain):
        return Polynomial(domain, basedomain, val)
    
    def __repr__(self):
        out = ""
        for i in range(len(self.coeffs)-1, -1,-1):
            if self[i]!=self.basedomain.zero or self.degree==0:
                if self[i]!=self.basedomain.one or i==0:
                    out += str(self[i])
                if i!=0:
                    out += str(self.domain.symbol)
                    if i!=1:
                        out += "^{}".format(i)
                out += "+"
        out = out.strip("+")
        return out
                
        