'''
Created on 30.03.2019

@author: Leonard
'''
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
        
