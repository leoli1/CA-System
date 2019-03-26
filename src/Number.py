'''
Created on 19.01.2019

@author: Leonard
'''
from __future__ import division
import AlgebraicStructures as AS
import NTheoryUtils as NTU

class Number(object):

    def __init__(self, domain):
        """
        domain should be a ring or a field
        """
        self.domain = domain
    
    def __add__(self, other):
        return self.domain.add(self,other)
    def __mul__(self, other):
        return self.domain.mul(self, other)
    def __sub__(self, other):
        return self.domain.sub(self,other)
    def __div__(self, other):
        return self/other
    def __truediv__(self, other):
        return self.domain.div(self, other)
    def __neg__(self):
        return self.domain.addInverse(self)
    
    @staticmethod
    def fromValue(val):
        raise NotImplementedError()
    
class RationalNumber(Number):
    
    def __init__(self, num,denom):#, domain=AS.F_Q):
        super(RationalNumber, self).__init__(AS.F_Q)
        self.num = num
        self.denom = denom
        if self.denom==0:
            raise ValueError("Denominator can't be zero.")
        if self.denom<0:
            self.num *= -1
            self.denom *= -1
        
        if self.num==0:
            self.denom=1
        else:
            self.removeCommonFactors()    
        self.num = int(self.num)
        self.denom = int(self.denom)
        
    def removeCommonFactors(self):
        gcd = NTU.gcd(abs(self.num),abs(self.denom))
        if gcd!=1:
            self.num /= gcd
            self.denom /= gcd
    @staticmethod
    def fromValue(val):
        power = 0
        num = val
        while int(num)!=num:
            num*=10
            power+=1
        return RationalNumber(num,10**power)
    def __repr__(self):
        if self.denom==1:
            return str(self.num)
        return "{}/{}".format(self.num,self.denom)
    
class ModularInteger(Number):
    def __init__(self, val, domain):
        super(ModularInteger, self).__init__(domain)
        self.val = val%domain.order
        
    @staticmethod
    def fromValue(val, domain):
        return ModularInteger(val, domain)
    
    def __repr__(self):
        return str(self.val%self.domain.order)
    
        