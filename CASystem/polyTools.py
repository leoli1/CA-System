'''
Created on 30.03.2019

@author: Leonard
'''
import AlgebraicStructures as AS
import Matrix


def rootsOverFiniteField(poly):
    field = poly.basedomain
    if poly.degree==0:
        return []
    if poly.degree==1:
        return [[-poly[0]/poly[1],1]]
        
    for i in range(field.order):
        n = field(i)
        if poly.evaluate(n)==field.zero:
            lin = poly.domain([-n,field.one])
            q = poly//lin
            roots = rootsOverFiniteField(q)
            for r in roots:
                if r[0]==n:
                    r[1]+=1
                    break
            else:
                roots.append([n,1])
            return roots
        
    return []

def isSquareFree(poly):
    return poly.domain.coprime(poly.derivative(),poly)

def SFF(f): # squarefree factorization
    factors = []
    c = f.domain.gcd(f.derivative(),f)
    w = f//c
    
    i=1
    while w!=f.domain.one:
        y = f.domain.gcd(w,c)
        fac = w//y
        factors.append((fac,i))
        w = y
        c = c//y
        i+=1
        
    if c!=f.domain.one: # for finite fields
        p = f.domain.basedomain.order
        newCoeffs = []
        for i in range(0,f.degree+1,p):
            newCoeffs.append(c[i])
        g = f.domain(newCoeffs)
        newFactors = SFF(g)
        for i in range(len(newFactors)):
            newFactors[i] = (newFactors[i][0],newFactors[i][1]*p)
        factors += newFactors
        factors.sort(key=lambda x: x[1])
        #raise NotImplementedError()
        
    ffactors = []
    for fac in factors:
        if fac[0]!=f.domain.one:
            ffactors.append(fac)
    return ffactors
    """if isSquareFree(poly):
        return [(poly,1)]
    c = poly.domain.gcd(poly.derivative(),poly)
    w = poly//c
    #if w==poly.domain.one:
     #   return []
    y = poly.domain.gcd(c,w)
    a1 = w//y
    rest = SFF(c)
    an = [(a1,1)]+rest
    for i in range(1,len(an)):
        an[i] = (an[i][0],an[i][1]+1)
    return an"""
def modPow(base, exp, mod):
    out = base.domain.one
    if base==out:
        return base
    while exp>0:
        if exp%2==1:
            exp -= 1
            out =(out*base)%mod
        else:
            base = (base*base)%mod
            exp/=2
    return out
def DDF(f):
    if not isSquareFree(f) or not f.isMonic():
        raise Exception()
    if not isinstance(f.domain.basedomain,AS.IntegerResidueClassRing):
        raise Exception()
    q = f.domain.basedomain.order
    i=1
    S=[]
    fs = f
    x = f.domain.getVar()
    while fs.degree>=2*i:
        g = f.domain.gcd(fs,modPow(x,q**i,fs)-(x%fs))
        if not g==f.domain.getAssociateRepresentant(f.domain.one)[0]:
            S.append((g,i))
            fs = fs//g
        i+=1
    
    if fs!=f.domain.one:
        S.append((fs,fs.degree))
    if S==[]:
        return (f,1)
    else:
        return S
    
def getBerlekampMatrix(f,p):
    return Matrix.Matrix(f.degree,f.degree,f.basedomain,[(modPow(f.domain.getVar(), p*i, f)).coeffs for i in range(f.degree)])
    
def BerlekampFactorization(f):
    if not isSquareFree(f) or not f.isMonic():
        raise Exception()
    if not isinstance(f.domain.basedomain,AS.IntegerResidueClassRing):
        raise Exception()
    p = f.domain.basedomain.order
    B = getBerlekampMatrix(f, p)
    nb = (B-Matrix.getUnitMatrix(f.degree, f.basedomain)).getLeftNullSpaceMatrix()
    r = len(nb)
    vs = []
    for i in range(r):
        vs.append(f.domain(nb[i][0]))
        
    V = []
    for x in vs:
        if x!=x.domain.one:
            V.append(x)
    F = [f]
    FI = []
    while len(F)<r and len(V)>0:
        v = V[0]
        V.remove(v)
        Fnew = []
        for fi in F:
            for j in range(0,p):
                gij = f.domain.gcd(fi,v-j)
                if gij!=gij.domain.one:
                    Fnew.append(gij)
        F = []
        for g in Fnew:
            if g.degree==1:
                FI.append(g)
            else:
                F.append(g)
        
    
    return F+FI
        
    
    