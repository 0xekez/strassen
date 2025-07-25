import z3
from itertools import product

multiplications = 7

ones = [
    (0,0,0),
    (1,2,0),

    (0,1,1),
    (1,3,1),

    (2,0,2),
    (3,2,2),

    (2,1,3),
    (3,3,3),
]
w = [[z3.Int(f'w_{l}^{k}') for l in range(4)] for k in range(multiplications)]
u = [[z3.Int(f'u_{i}^{k}') for i in range(4)] for k in range(multiplications)]
v = [[z3.Int(f'v_{j}^{k}') for j in range(4)] for k in range(multiplications)]

equations = [
    (1 if (i,j,l) in ones else 0) == z3.Sum([w[k][l]*u[k][i]*v[k][j] for k in range(multiplications)])
    for i,j,l in product(range(4),range(4),range(4))
]

s = z3.Solver()
for k in range(multiplications):
    for t in (u[k]+v[k]+w[k]):
        s.add(z3.Or(t==-1, t==0, t==1))
for e in equations:
    s.add(e)

# this is a 'cheat' to assert only two matrix elements appear in each
# summation. we know this doesn't make the problem unsat because the
# strassen solution has this property. without this the SMT solver
# takes too long.
for k in range(multiplications):
    s.add(z3.Sum([z3.If(u[k][i]!=0,1,0) for i in range(4)]) <= 2)
    s.add(z3.Sum([z3.If(v[k][j]!=0,1,0) for j in range(4)]) <= 2)

assert s.check()==z3.sat
m = s.model()

for k in range(multiplications):
    a = [f'a_{i}' if m[u[k][i]].as_long()==1 else f'-a_{i}' for i in range(4) if m[u[k][i]].as_long()!=0]
    b = [f'b_{i}' if m[v[k][i]].as_long()==1 else f'-b_{i}' for i in range(4) if m[v[k][i]].as_long()!=0]
    print(f'P_{k}=({" + ".join(a)})({" + ".join(b)})')
