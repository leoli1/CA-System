'''
Created on 27.03.2019

@author: Leonard
'''

class Symbol(object):
    AllSymbols = None
    def __new__(self,name):
        if Symbol.AllSymbols==None:
            Symbol.AllSymbols = []
        for sym in Symbol.AllSymbols:
            if sym.name==name:
                return sym
            
        obj = object.__new__(Symbol)
        obj.name = name
        Symbol.AllSymbols.append(obj)
        return obj
        
    def __repr__(self):
        return self.name
    def __eq__(self, other):
        if type(other)!=Symbol:
            return False
        return self.name==other.name
        