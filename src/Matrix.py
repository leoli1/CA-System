'''
Created on 19.01.2019

@author: Leonard
'''
from Utils import iterMap

class ElementaryRowOperation(object):
    pass
class RowOperation1(ElementaryRowOperation):
    def __init__(self, matrix, scalar, row):
        self.matrix = matrix
        self.scalar = scalar
        if scalar==self.matrix.domain.getZero():
            raise Exception()
        self.row = row
        
class RowOperation2(ElementaryRowOperation):
    def __init__(self, matrix, row1, row2):
        self.matrix = matrix
        self.row1 = row1
        self.row2 = row2
class RowOperation3(ElementaryRowOperation):
    def __init__(self, matrix, row1, row2, factor):
        self.matrix = matrix
        self.row1 = row1
        self.row2 = row2
        self.factor = factor
        if row1==row2:
            raise Exception()
        
class Matrix(object):

    def __init__(self, rows, columns, domain,data=None):

        if data==None:
            self.__data = [[domain.getZero() for j in range(columns)] for i in range(rows)]
        else:
            self.__data = iterMap(domain.elementFromValue, data)
        self.rows = rows
        self.columns = columns
        self.domain = domain
        
    def __setitem__(self, index, value):
        if type(index)==tuple:
            self.__data[index[0]][index[1]] = value
        else:
            self.__data[index] = map(self.domain.elementFromValue, value)
    def __getitem__(self, index):
        if type(index)==tuple:
            return self.__data[index[0]][index[1]]
        return self.__data[index]
    
    def __add__(self, other):
        if self.rows!=other.rows or self.columns!=other.columns:
            raise MatrixMultiplicationError("Wrong dimensions")
        if self.domain!=other.domain:
            raise MatrixMultiplicationError("Wrong domains")
        out = Matrix(self.rows,self.columns,self.domain)
        for i in range(self.rows):
            for j in range(self.columns):
                #out[i,j] = self.domain.add(self[i,j],other[i,j])
                out[i,j] = self[i,j]+other[i,j]
        return out
    
    def __mul__(self, other):
        if self.columns!=other.rows:
            raise MatrixMultiplicationError("Wrong dimensions")
        if self.domain!=other.domain:
            raise MatrixMultiplicationError("Wrong domains")
        out = Matrix(self.rows,other.columns, self.domain)
        for i in range(self.rows):
            for j in range(other.columns):
                s = self.domain.getZero()
                for k in range(self.columns):
                    #s = self.domain.add(s,self.domain.mul(self[i,k],other[k,j]))
                    s = s+self[i,k]*other[k,j]
                out[i,j] = s
        return out
    
    def __repr__(self):
        strings = reduce(list.__add__, [map(str, self[i]) for i in range(self.rows)])
        maxLength = len(max(strings, key=lambda s: len(s)))
        out = ""
        for i in range(self.rows):
            out+="["
            for j in range(self.columns):
                s = strings[j+self.columns*i]
                out += s+" "*(maxLength-len(s))
                if j!=self.columns-1:
                    out+="|"
            out+="]"
            if i!=self.rows-1:
                out += "\n"
        return out
    
class MatrixMultiplicationError(Exception):
    pass

if __name__=="__main__":
    import AlgebraicStructures
    Q = AlgebraicStructures.F_Q
    Z5 = AlgebraicStructures.ResidueClassField(5)
    m1 = Matrix(2,3,Q,data=[[1,2,3],[42,-1,0]])
    m2 = Matrix(3,4,Q,data=[[0.5,0,0,1],[0,1,1,-1],[100,20,0,-1]])
    m3 = Matrix(2,2,Z5,data=[[1,3],[0,2]])
    m4 = Matrix(2,2,Z5,data=[[4,0],[3,3]])
    print(m1)
    print(m2)
    print(m1*m2)
    print("---")
    print(m3)
    print(m4)
    print("-")
    print(m3+m4)
    print(".")
    print(m3*m4)
    print(".")
    print(m4*m3)
        