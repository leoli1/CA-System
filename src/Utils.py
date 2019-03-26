'''
Created on 26.03.2019

@author: Leonard
'''

def iterMap(func, lst):
    return map(lambda l: map(func, l), lst)