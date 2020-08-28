#!/usr/bin/env python3
# coding=utf-8

from gurobipy import *

'''
According to Gurobi manual, 
PoolSearchMode = 1 is to find additional solutions, but
with no guarantees about the quality of those solutions.
PoolSearchMode = 2 is to find the n best solutions, where n is determined by the value
of the PoolSolutions parameter.

It seems we should use PoolSearchMode = 1, 
however, the presolve process of Gurobi may throw some solutions. 
We cannot get all the solutions if we only set PoolSearchMode = 1.
We have to set Presolve = 0 (close it) simultaneously. 

Another method is to set PoolSearchMode = 2, 
we confirm its correctness from Dr. Zonghao Gu, one of the co-founder of Gurobi. 
'''
# PoolSearchMode = 2
m = Model()
m.setParam( 'OutputFlag', 0 )
m.setParam( 'PoolSearchMode',  2 )
m.setParam( 'PoolSolutions', 2000000000 )

print()
print( "*" * 50 )
print()
print( "PoolSearchMode = {}".format( 2 ) )

a = m.addVar(vtype = GRB.BINARY) 
b = m.addVar(vtype = GRB.BINARY) 
c = m.addVar(vtype = GRB.BINARY) 
d = m.addVar(vtype = GRB.BINARY)
e = m.addVar(vtype = GRB.BINARY)
f = m.addVar(vtype = GRB.BINARY)

m.addConstr( a + b + c + d + e + f == 1 ) 

m.optimize()

counter = m.getAttr( 'SolCount' ) 
print( "{} solutions...".format( counter ) )

for i in range( counter ):
    m.setParam( 'SolutionNumber', i )
    #print( "a: {}  b: {}  c: {}  d: {}  e: {}  f: {}".format( a.Xn, b.Xn, c.Xn, d.Xn, e.Xn, f.Xn ) )
    #print( "Round number..." )
    print( "a: {}  b: {}  c: {}  d: {}  e: {}  f: {}".format( 
        round( a.Xn ), round( b.Xn ), round( c.Xn ), round( d.Xn ), round( e.Xn ), round ( f.Xn ) ) )
    

print()
print( "*" * 50 )
print()

# PoolSearchMode = 1, Presolve is defautly 1
m = Model()
m.setParam( 'OutputFlag', 0 )
m.setParam( 'PoolSearchMode',  1 )
m.setParam( 'PoolSolutions', 2000000000 )
print( "PoolSearchMode = {}".format( 1 ) )

a = m.addVar(vtype = GRB.BINARY) 
b = m.addVar(vtype = GRB.BINARY) 
c = m.addVar(vtype = GRB.BINARY) 
d = m.addVar(vtype = GRB.BINARY)
e = m.addVar(vtype = GRB.BINARY)
f = m.addVar(vtype = GRB.BINARY)

m.addConstr( a + b + c + d + e + f == 1 ) 

m.optimize()

counter = m.getAttr( 'SolCount' ) 
print( "{} solutions...".format( counter ) )

for i in range( counter ):
    m.setParam( 'SolutionNumber', i )
    #print( "a: {}  b: {}  c: {}  d: {}  e: {}  f: {}".format( a.Xn, b.Xn, c.Xn, d.Xn, e.Xn, f.Xn ) )
    #print( "Round number..." )
    print( "a: {}  b: {}  c: {}  d: {}  e: {}  f: {}".format( 
        round( a.Xn ), round( b.Xn ), round( c.Xn ), round( d.Xn ), round( e.Xn ), round ( f.Xn ) ) )

print()
print( "*" * 50 )
print()

# PoolSearchMode = 1, Presolve = 0
m = Model()
m.setParam( 'OutputFlag', 0 )
m.setParam( 'Presolve', 0 )
m.setParam( 'PoolSearchMode',  1 )
m.setParam( 'PoolSolutions', 2000000000 )
print( "PoolSearchMode = {}, Presolve = {}".format( 1, 0 ) )

a = m.addVar(vtype = GRB.BINARY) 
b = m.addVar(vtype = GRB.BINARY) 
c = m.addVar(vtype = GRB.BINARY) 
d = m.addVar(vtype = GRB.BINARY)
e = m.addVar(vtype = GRB.BINARY)
f = m.addVar(vtype = GRB.BINARY)

m.addConstr( a + b + c + d + e + f == 1 ) 

m.optimize()

counter = m.getAttr( 'SolCount' ) 
print( "{} solutions...".format( counter ) )

for i in range( counter ):
    m.setParam( 'SolutionNumber', i )
    #print( "a: {}  b: {}  c: {}  d: {}  e: {}  f: {}".format( a.Xn, b.Xn, c.Xn, d.Xn, e.Xn, f.Xn ) )
    #print( "Round number..." )
    print( "a: {}  b: {}  c: {}  d: {}  e: {}  f: {}".format( 
        round( a.Xn ), round( b.Xn ), round( c.Xn ), round( d.Xn ), round( e.Xn ), round ( f.Xn ) ) )

print()
print( "*" * 50 )
