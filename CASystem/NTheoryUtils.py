'''
Created on 26.03.2019

@author: Leonard
'''

def isPrime(a):
    for i in range(2,a/2+1):
        if a%i==0:
            return False
    return True