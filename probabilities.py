import operator as op
from functools import reduce

## n: Number of nodes
def factorial(n):
    return reduce(op.mul, range(1, n+1), 1)

## Number of r-combinations of n elements
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = factorial(r)
    return numer // denom  

## Function P_adjacent
def prob_adj(ki,kj,m):
    return sum([((-1)**(l+1))*ncr(ki,l)*ncr(kj,l)*factorial(l)/(reduce(op.mul, [2*m+1-2*s for s in range(1,l+1)], 1)) for l in range(1,min(ki,kj)+1)])

## Function P_common_neighbor
def prob_CN(ki,kj,ks,m):
    return sum([((-1)**(l+1))*prob_adj(kj,ks-l,m-l)*ncr(ki,l)*ncr(ks,l)*factorial(l)/(reduce(op.mul, [2*m+1-2*s for s in range(1,l+1)], 1)) for l in range(1,min(ki,ks)+1)])

## Function P_triangle
def prob_triangle(ki,kj,ks,m):
    return sum([((-1)**(l+1))*prob_CN(ki-l,kj-l,ks,m-l)*ncr(ki,l)*ncr(kj,l)*factorial(l)/(reduce(op.mul, [2*m+1-2*s for s in range(1,l+1)], 1)) for l in range(1,min(ki,kj)+1)])