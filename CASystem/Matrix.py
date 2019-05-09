'''
Created on 19.01.2019

@author: Leonard
'''
from Utils import iterMap
import Symbol
import AlgebraicStructures as AS


class ElementaryRowOperation(object):
    def apply(self, matrix):
        raise NotImplementedError()
class RowOperation1(ElementaryRowOperation):
    def __init__(self, domain, row, factor):
        self.basedomain = domain
        self.factor = factor
        if factor==self.basedomain.zero:
            raise Exception()
        self.row = row
    def apply(self, matrix):
        m = matrix.copy()
        for j in range(m.columns):
            m[self.row,j]*=self.factor
        return m
    def __repr__(self):
        return "R{}->{}*R{}".format(self.row,self.factor,self.row)
        
class RowOperation2(ElementaryRowOperation):
    def __init__(self, domain, row1, row2, factor):
        self.basedomain = domain
        self.row1 = row1
        self.row2 = row2
        self.factor = factor
        if row1==row2:
            raise Exception()
        
    def apply(self, matrix):
        m = matrix.copy()
        for j in range(m.columns):
            m[self.row2,j] += m[self.row1,j]*self.factor
        return m
    def __repr__(self):
        return "R{}->R{}+{}*R{}".format(self.row2, self.row2, self.factor, self.row1)
            
class RowOperation3(ElementaryRowOperation):
    def __init__(self, domain, row1, row2):
        self.basedomain = domain
        self.row1 = row1
        self.row2 = row2
    def apply(self, matrix):
        m = matrix.copy()
        itemp = matrix[self.row1]
        jtemp = matrix[self.row2]
        m[self.row1] = jtemp
        m[self.row2] = itemp
        return m
    def __repr__(self):
        return "R{}<->R{}".format(self.row1,self.row2)
        
        
        
class Matrix(object):

    def __new__(cls, rows, columns, basedomain,data=None):
        
        obj = object.__new__(cls)
        if rows==columns:# and _makeSquareOpt:
            obj = object.__new__(SquareMatrix)
            obj._init(rows,basedomain, data)#object.__new__(SquareMatrix, rows,basedomain,data) 
            return obj
        obj._init(rows,columns,basedomain,data)
        return obj
    #def __init__(self,data=None,*args):
    #    pass
    def _init(self, rows, columns, basedomain, data=None):
        if data==None:
            self.__data = [[basedomain.zero for j in range(columns)] for i in range(rows)]
        else:
            self.__data = iterMap(basedomain.elementFromValue, data)
            for i in range(rows):
                if len(self.__data[i])<columns:
                    self.__data[i] += [basedomain.zero]*(columns-len(self.__data[i]))
        self.rows = rows
        self.columns = columns
        self.basedomain = basedomain
        self.__rrefTransform = None
        
    def __setitem__(self, index, value):
        if type(index)==tuple:
            self.__data[index[0]][index[1]] = self.basedomain.elementFromValue(value)
        else:
            self.__data[index] = map(self.basedomain.elementFromValue, value)
    def __getitem__(self, index):
        if type(index)==tuple:
            return self.__data[index[0]][index[1]]
        return self.__data[index]
    def getRow(self, index):
        return Matrix(1,self.columns,self.basedomain, data=[self[index]])
    
    def __add__(self, other):
        if self.rows!=other.rows or self.columns!=other.columns:
            raise MatrixMultiplicationError("Wrong dimensions")
        if self.basedomain!=other.basedomain:
            raise MatrixMultiplicationError("Wrong domains")
        out = Matrix(self.rows,self.columns,self.basedomain)
        for i in range(self.rows):
            for j in range(self.columns):
                #out[i,j] = self.basedomain.add(self[i,j],other[i,j])
                out[i,j] = self[i,j]+other[i,j]
        return out
    def __sub__(self, other):
        return self+other.scalarMultiplication(-self.basedomain.one)
    
    def __mul__(self, other):
        if self.columns!=other.rows:
            raise MatrixMultiplicationError("Wrong dimensions")
        if self.basedomain!=other.basedomain:
            raise MatrixMultiplicationError("Wrong domains")
        out = Matrix(self.rows,other.columns, self.basedomain)
        for i in range(self.rows):
            for j in range(other.columns):
                s = self.basedomain.zero
                for k in range(self.columns):
                    #s = self.basedomain.add(s,self.basedomain.mul(self[i,k],other[k,j]))
                    s = s+self[i,k]*other[k,j]
                out[i,j] = s
        return out
    
    def scalarMultiplication(self, scalar):
        m = self.copy()
        for i in range(self.rows):
            for j in range(self.columns):
                m[i,j]*=scalar
        return m
    
    def __neg__(self):
        return self.scalarMultiplication(self.basedomain.one.addInverse())
    
    def isZero(self):
        for i in range(self.rows):
            for j in range(self.columns):
                if self[i,j] != self.basedomain.zero:
                    return False
        return True
    
    def rank(self):
        return self.getRREF()[1]
    def nullity(self):
        return self.rows-self.getRREF()[1]
    def hasFullRank(self):
        (rref, rank, rowOps) = self.getRREF()
        return rank == min(self.rows,self.columns)
    
    def determinant(self, alg="gauss"):
        """
        returns the determinant of self
        """
        if self.rows!=self.columns:
            raise Exception("no square matrix")
        if alg=="gauss":
            (_, rank, rowOps,_) = self.getRREF()
            if rank!=self.rows:
                return self.basedomain.zero
            det = self.basedomain.one#sum(self[i,i] for i in range(self.rows), self.basedomain.zero
            for ops in rowOps:
                if type(ops)==RowOperation1:
                    det /= ops.factor
                elif type(ops)==RowOperation3:
                    det *= -self.basedomain.one
            return det
        
    def invert(self):
        """
        returns the inverse of self, if self is invertible
        """
        if self.columns!=self.rows:
            raise MatrixInversionError("no square matrix")
        (_, rank, rowOps,_) = self.getRREF()
        if rank!=self.columns:
            raise MatrixInversionError("not full rank -> not invertible")
        m = getUnitMatrix(rank, self.basedomain)
        for ops in rowOps:
            m = ops.apply(m)
        return m
    def charPoly(self):
        """
        returns the characteristic polynomial p of self=A, p = det(tI-A)
        """
        if self.columns!=self.rows:
            raise MatrixInversionError("no square matrix")
        polyDomain = AS.PolynomialDomain(self.basedomain, Symbol.Symbol("t"))
        rationalsDomain = AS.RationalFunctionsDomain(polyDomain)
        charMatrix = Matrix(self.rows, self.columns, rationalsDomain)
        t = polyDomain([self.basedomain.zero, self.basedomain.one])#rationalsDomain(polyDomain([self.basedomain.zero, self.basedomain.one]),polyDomain.one)
        for i in range(self.rows):
            for j in range(self.columns):
                charMatrix[i,j] = polyDomain([-self[i,j]])#rationalsDomain(polyDomain([-self[i,j]]),polyDomain.one)
                if i==j:
                    charMatrix[i,j] += t
                    
        return charMatrix.determinant().asBaseRingElement()
    
    def getNullspaceBasis(self):
        """
        returns a basis for the nullspace of self
        """
        (rref, _,_,pivots) = self.getRREF()
        basis = []
        for i in range(self.columns):
            if i not in pivots:
                vector = Matrix(self.columns,1,self.basedomain)
                vector[i,0]=self.basedomain.one
                for j in range(i):
                    if j in pivots:
                        vector[j,0] = -rref[pivots.index(j),i]
                basis.append(vector)
                
        return basis
    
    def getLeftNullSpaceMatrix(self):
        trans = self.transpose()
        nb = trans.getNullspaceBasis()
        return map(Matrix.transpose,nb)
    
    def getEigenValues(self):
        """
        returns a list [[l1,a1,g1],[l2,a2,g2],...] where the li are the different eigenvalues and the ai/gi 
        are their algebraic/geometric multiplicities        
        """
        cpoly = self.charPoly()
        roots = cpoly.findRoots()
        evs = []
        for (root,algMult) in roots:
            evs.append([root,algMult,(getUnitMatrix(self.rows, self.basedomain).scalarMultiplication(root)-self).nullity()])
        evs.reverse()
        return evs
    def isDiagonalizable(self):
        evs = self.getEigenValues()
        return sum(x[2] for x in evs)==self.rows
    def isTrigonalizable(self):
        evs = self.getEigenValues()
        return sum(x[1] for x in evs)==self.rows
    def transpose(self):
        m = Matrix(self.columns,self.rows,self.basedomain)
        for i in range(self.rows):
            for j in range(self.columns):
                m[j,i] = self[i,j]
        return m
    def getDiagonalizerMatrix(self):
        """
        returns (P,D), such that D is diagonal and self = PDP^(-1)
        """
        if not self.isDiagonalizable():
            raise Exception("Not diagonalizable")
        
        eigenVectors = []
        eigenValues = []
        for (ev,m,_) in self.getEigenValues():
            eigenVectors += (getUnitMatrix(self.rows, self.basedomain).scalarMultiplication(ev)-self).getNullspaceBasis()
            eigenValues += [ev]*m
        return (Matrix.getMatrixFromColumnVectors(eigenVectors), Matrix.diag(self.basedomain, eigenValues))
        
    

    def getRREF(self):
        """
        returns (rref,rank,rowOps,pivots) where
            - rref is the reduced row echelon form of self
            - rank is the rank of self
            - rowOps are the rowOperations that took it to transfrom self into rref
            - pivots are the column-indices of the leading 1s in rref
        """
        if self.__rrefTransform!=None:
            return self.__rrefTransform
        r=0
        l=0
        pivots = []
        rowOperations = []
        if self.isZero():
            return (self, r,rowOperations,pivots)
        m = self.copy()
        (m, rowOperations) = self._RREF_Helper_moveZeroRows(m)
        
        while l<self.rows:
            jl = -1
            k = -1
            for j in range(self.columns):
                (nonzero, k) =  self._RREF_Helper_hasNonZeroEntriesAfterRow(m, j, l)
                if nonzero:
                    jl=j
                    break
            if jl==-1:
                raise Exception()
            pivots.append(jl)
            if k!=l:
                rowOps = RowOperation3(m.basedomain, k,l)
                m = rowOps.apply(m)
                rowOperations.append(rowOps)
            if m[l,jl]!=m.basedomain.one:
                rowOps = RowOperation1(m.basedomain, l, m[l,jl].inverse())
                m = rowOps.apply(m)
                rowOperations.append(rowOps)
            for k in range(self.rows):
                if k==l:
                    continue
                if m[k,jl]==m.basedomain.zero:
                    continue
                rowOps = RowOperation2(m.basedomain, l, k, -m[k,jl])
                m = rowOps.apply(m)
                rowOperations.append(rowOps)
            (m, nrowOps) = self._RREF_Helper_moveZeroRows(m)
            rowOperations += nrowOps
            for i in range(l+1, m.rows):
                if not m.getRow(i).isZero():
                    l+=1
                    break
            else:
                r=l+1
                break
            
        self.__rrefTransform = (m,r,rowOperations,pivots)
        return self.__rrefTransform
            
        
    def _RREF_Helper_moveZeroRows(self,m):
        swaps = [] # type 3
        for i1 in range(m.rows):
            for i in range(i1,m.rows):
                if m.getRow(i).isZero():
                    if i<(m.rows-1) and not m.getRow(i+1).isZero():
                        rowops = RowOperation3(m.basedomain, i,i+1)
                        m = rowops.apply(m)
                        swaps.append(rowops)
        return (m, swaps)
    
    def _RREF_Helper_hasNonZeroEntriesAfterRow(self, m, column, row):
        for i in range(row,self.rows):
            if not m[i,column]==m.basedomain.zero:
                return (True, i)
            
        return (False, -1)
        
    
    @staticmethod
    def diag(basedomain, diagElements):
        """
        returns matrix with diagElements on the diagonal
        """
        n = len(diagElements)
        m = Matrix(n,n,basedomain)
        for i in range(n):
            m[i,i] = diagElements[i]
        return m
    
    @staticmethod
    def getMatrixFromColumnVectors(vecs):
        """
        given vectors vecs=[v1,v2,...,vn]
        returns Matrix m=[v1 v2 ... vn]
        """
        m = Matrix(vecs[0].rows, len(vecs), vecs[0].basedomain)
        for i in range(m.rows):
            for j in range(m.columns):
                m[i,j] = vecs[j][i,0]
        return m
    
    def copy(self):
        m = Matrix(self.rows,self.columns,self.basedomain)
        for i in range(self.rows):
            for j in range(self.columns):
                m[i,j]=self[i,j]
        return m
    
    def __eq__(self, other):
        if id(self)==id(None):
            return False
        if self.rows!=other.rows or self.columns!=other.columns:
            return False
        for i in range(self.rows):
            for j in range(self.columns):
                if self[i,j]!=other[i,j]:
                    return False
        return True
    
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
    
import DomainElement
class SquareMatrix(Matrix, DomainElement.DomainElement):
    def __init__(self,a=None,b=None,c=None,data=None,*args):
    #def __init__(self,*args):
        pass
    def _init(self, size, basedomain, data=None, domain=None):
        super(SquareMatrix, self)._init(size, size, basedomain, data)
        if domain!=None:
            self.domain = domain
        self.size = size

def getUnitMatrix(size, domain):
    one = Matrix(size, size, domain)
    for i in range(size):
        one[i,i]=domain.one
    return one

class MatrixMultiplicationError(Exception):
    pass
class MatrixInversionError(Exception):
    pass

if __name__=="__main__":
    Q = AS.F_Q
    Z5 = AS.IntegerResidueClassField(5)
    m1 = Matrix(2,3,Q,data=[[1,2,3],[42,-1,0]])
    m2 = Matrix(3,4,Q,data=[[0,0,0,1],[0,1,1,-1],[100,20,0,-1]])
    m3 = Matrix(2,2,Z5,data=[[1,3],[0,2]])
    m4 = Matrix(2,2,Z5,data=[[4,0],[3,3]])
    """print(m1)
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
    print(m4*m3)"""
    m = Matrix(3,3,Q,data=[[2,2,4],[0,1,2],[0,0,1]])
    print(m)
    print(m.invert())
        