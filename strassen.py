import z3
from itertools import product

def lex_leq(a, b):
    if not a:
        return True
    return z3.Or(a[0] <  b[0],
                z3.And(a[0] == b[0],
                       lex_leq(a[1:], b[1:])))

multiplications = 7
a_dim = [2,2]
b_dim = [2,2]
assert a_dim[1]==b_dim[0]
c_dim = [a_dim[0],b_dim[1]]

a_size = a_dim[0]*a_dim[1]
b_size = b_dim[0]*b_dim[1]
c_size = c_dim[0]*c_dim[1]

gamma = [[z3.Int(f'gamma_{l}^{k}') for l in range(c_size)] for k in range(multiplications)]
alpha = [[z3.Int(f'alpha_{i}^{k}') for i in range(a_size)] for k in range(multiplications)]
beta = [[z3.Int(f'beta_{j}^{k}') for j in range(b_size)] for k in range(multiplications)]

# with row-major indexing: T_ijl is 1 iff a_ib_j appears in the
# standard formula for c_l. the code below is derived by converting
# indexing in the standard formula for c_ij into row-major.
ones = []
for i,j in product(range(c_dim[0]),range(c_dim[1])):
    ones += [(i*a_dim[1]+k,k*b_dim[0]+j,i*c_dim[1]+j) for k in range(a_dim[1])]

terms = [
    (1 if (i,j,l) in ones else 0) == z3.Sum([gamma[k][l]*alpha[k][i]*beta[k][j] for k in range(multiplications)])
    for i,j,l in product(range(a_size),range(b_size),range(c_size))
]

for k in range(multiplications):
    for t in (alpha[k]+beta[k]+gamma[k]):
        terms += [z3.Or(t==-1, t==0, t==1)]

# one may flip the sign of P_k and the sign of all gamma^k_l to create a
# symmetry. P_k is (Σ)(Σ) is (-1)(Σ)(-1)(Σ), another symmetry. only
# one of (-1)(Σ) and (Σ) will have a positive first term. by prefering
# one, we prefer a sign for P_k, breaking the symmetries.
for k in range(multiplications):
    terms += [
        z3.Or(*[
            z3.And(alpha[k][i]==1,
                   *[alpha[k][j]==0 for j in range(i)])
            for i in range(a_size)
    ])]
    terms += [
        z3.Or(*[
            z3.And(beta[k][i]==1,
                   *[beta[k][j]==0 for j in range(i)])
            for i in range(b_size)
    ])]

# P_k should be lexographically ordered, otherwise permuting the order
# of products creates a symmetry.
for k in range(multiplications-1):
    terms += [lex_leq(alpha[k]+beta[k], alpha[k+1]+beta[k+1])]

## Optionally, search for only "Strassen-like" solutions.
# for k in range(multiplications):
#     s.add(z3.Sum([z3.If(u[k][i]!=0,1,0) for i in range(4)]) <= 2)
#     s.add(z3.Sum([z3.If(v[k][j]!=0,1,0) for j in range(4)]) <= 2)

s = z3.Solver()
for e in terms:
    s.add(e)

while s.check() == z3.sat:
    m = s.model()
    s.add(z3.Or([v() != m[v] for v in m]))
    for k in range(multiplications):
        a = [f'a_{i}' if m[alpha[k][i]].as_long()==1 else f'-a_{i}' for i in range(a_size) if m[alpha[k][i]].as_long()!=0]
        b = [f'b_{i}' if m[beta[k][i]].as_long()==1 else f'-b_{i}' for i in range(b_size) if m[beta[k][i]].as_long()!=0]
        print(f'P_{k}=({" + ".join(a)})({" + ".join(b)})')
    for l in range(c_size):
        p = [f'P_{k}' if m[gamma[k][l]].as_long()==1 else f'-P_{k}' for k in range(multiplications) if m[gamma[k][l]].as_long()!=0]
        print(f'c_{l}=({" + ".join(p)})')
    print()
