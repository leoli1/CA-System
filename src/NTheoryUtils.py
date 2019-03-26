'''
Created on 26.03.2019

@author: Leonard
'''
def euclid(a,b):
    r = a % b
    s = (a-r)/b
    if r==1:
        return (1,-s)
    (x,y) = euclid(b,r)
    return (y,x-y*s)

def gcd(a,b):
    (A,B) = (a,b) if a>=b else (b,a)
    if B==0:
        return A
    r = A%B  # A=s*B+r
    if r==0:
        return B
    else:
        return gcd(B,r)
    
def isPrime(a):
    for i in range(2,a/2+1):
        if a%i==0:
            return False
    return True