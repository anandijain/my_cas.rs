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

ode_order_lowering(ex: Ex) -> Vec<Ex> ?
ie find all higher order 

ex1 = der(der(x)) - Tx
ode_order_owering(ex1) ->

vec![
der(x) - d_x, 
der(d_x) - Tx,
]
