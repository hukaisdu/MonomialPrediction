# Monomial Prediction
This is the repository for the paper "An Algebraic Formulation of the Division Property: Revisiting Degree Evaluations, Cube Attacks, and Key-Independent Sums" accepted by the Asiacrypt 2020.

## Contents

1.  [Code for recovering the superpoly for 840-, 841- and 842-round Trivium](https://github.com/hukaisdu/MonomialPrediction/tree/master/Cube)
2. [Code for compute the exact degree of Trivium up to 834 rounds](https://github.com/hukaisdu/MonomialPrediction/tree/master/Degree)
3. [A demo for PoolSearchMode by a + b + c + d + e + f = 1 example](https://github.com/hukaisdu/MonomialPrediction/blob/master/demo_poolsearchMode.py)
4. [Superpolies that contain two many monomials (PI6, PI7, PI9, PI10) used in Section 5](https://github.com/hukaisdu/MonomialPrediction/blob/master/superpoly.pdf)

 ## Usage of the Codes in Cube and Degree

### Dependencies

To run our code, you should first install the [Gurobi solver](https://www.gurobi.com) and set the proper license. 

### Compile 

After you install the solver, then you need to edit Makefile to modify --lgurobixx to your own version. Then type 

`make`

to compilen our code.

### Run the Code

#### Compute the Exact Degree

After you compile the code, just type 

`./trivium`

to compute the degree of Trivium up to 834 rounds.

#### Recover the Superpoly

After you compile the code, type 

`./trivium [ROUND] [INDEX]`  

The possible combinations of (ROUND, INDEX) are listed as follows, 
1. ROUND = 840, INDEX = 1:
    Recover the superpoly for [0,1,...,79]/{70, 72, 74, 76, 78} of 840-round Trivium

2. ROUND = 840, INDEX = 2:
    Recover the superpoly for [0,1,...,79]/{72, 74, 76, 78} of 840-round Trivium

3. ROUND = 840, INDEX = 3:
    Recover the superpoly for [0,1,...,79]/{70, 74, 76, 78} of 840-round Trivium

4. ROUND = 841, INDEX = 1:
    Recover the superpoly for [0,1,...,79]/{70, 72, 76, 78} of 841-round Trivium

5. ROUND = 841, INDEX = 2:
    Recover the superpoly for [0,1,...,79]/{72, 76, 78} of 841-round Trivium

6. ROUND = 842, INDEX = 1:
    Recover the superpoly for [0,1,...,79]/{72, 74, 76, 78} of 842-round Trivium

7. ROUND = 842, INDEX = 2:
    Recover the superpoly for [0,1,...,79]/{74, 76, 78} of 842-round Trivium

