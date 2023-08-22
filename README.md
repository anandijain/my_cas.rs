# my_cas.rs

a basic cas that is used to generate ODE functions used by cvode-wrap/sundials-sys (eventually diffeq.rs)

the goal is to add a spice parser to do real-time circuit simulation of audio circuits


### notes on pantelides paper 
#### SDR walkthrough
A = {a, b, c}

A1 = {a, b} 
A2 = {a, c}
A3 = {b, c}

m = 3 

I = [1]
A1 = a,b, |Ai| = 2 |I| = 1
I = [2]
A2 = a,c, |Ai| = 2 |I| = 1
I = [3]
A3 = b,c, |Ai| = 2 |I| = 1
 
I = [1,2]
A1 U A2 = {abc} |Ai| = 3 |I| = 2 
I = [1,3]
A1 U A3 = {abc} |Ai| = 3 |I| = 2

I = [123]
A1 U A2 U A3 = {abc} |Ai| = 3 |I| = 3

=> {A1,A2,A3} has an SDR 

V = 1:6 with X = 1:3 Y = 4:6
X is the subsets (equations)
Y is the variables (states)

then the bipartite graph is 
1 -> 4
1 -> 5
2 -> 4
2 -> 6
3 -> 5
3 -> 6

trivial vertex cover is always just V 
so i think 1:3 is a vertex cover of g 
similarly for 4:6 

a (maximum) matching of g would be edges {(1, 4) (2, 6) (3, 5)}
so technically a matching is {} (no edges?)


order lowering walkthrou

sys:
1) Tx - der(x, 2)
2) Ty - g - der(y, 2)
3) x^2 + y^2 - L^2 

-> 
 
1) der(x, 1) - w 
2) der(y, 1) - z 
3) der(w, 1) - Tx
4) der(z, 1) - Ty - g
5) x^2 + y^2 - L^2

X = [x, y, w, z, der(x), der(y), der(w), der(z), T]
A = [5, 6, 7, 8, 0,      0,      0,      0,      0]
