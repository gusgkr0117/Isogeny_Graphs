import json
# p + 1 = 20 = 2^2 * 5
p = 659
Fp.<i> = GF(p^2, modulus = x^2 + 1)
ell = 5

graph_nodes = set([])
graph_links = set([])

total_visit = [Fp(-1)]
edges_storage = []
unique_node = {}
BFS_depth = 10
# Curve, Odd Prime Degree(n != 2), Direction([0,n])
def Isog(A, n, x, eval_point = 0):
    A = Fp(A)
    E = EllipticCurve(Fp, [0, A, 0, 1, 0])
    P, Q = E.gens()
    P, Q = ((p+1)//n) * P, ((p+1)//n) * Q
    if x==n: R = Q
    else: R = P + x * Q

    if eval_point ==0 : eval_point = P[0], Q[0]

    hs_1, hs_2 = Fp(1), Fp(1)
    hs_P1, hs_P2 = Fp(1), Fp(1)
    hs_Q1, hs_Q2 = Fp(1), Fp(1)
    T = R
    for _ in range((n-1)//2):
        hs_1 *= (1-T[0])
        hs_2 *= (-1-T[0])
        hs_P1 *= (eval_point[0] - T[0])
        hs_P2 *= (eval_point[0]^(-1) - T[0])
        hs_Q1 *= (eval_point[1] - T[0])
        hs_Q2 *= (eval_point[1]^(-1) - T[0])

        T += 2*R

    d = (((A-2) / (A+2))^n) * (hs_1 / hs_2)^8
    A_co = 2*((1+d)/(1-d))
    A_co_red = ReducedMontCoeff(A_co)

    E_co = EllipticCurve(Fp, [0, A_co, 0, 1, 0])
    E_co_red = EllipticCurve(Fp, [0, A_co_red, 0, 1, 0])

    eval_P = E.lift_x(eval_point[0])
    eval_Q = E.lift_x(eval_point[1])

    phi_P, phi_Q = E_co(0), E_co(0)
    if eval_P.weil_pairing(R, n) != 1 : phi_P = E_co.lift_x(eval_point[0]^n * (hs_P2 / hs_P1)^2)
    if eval_Q.weil_pairing(R, n) != 1 : phi_Q = E_co.lift_x(eval_point[1]^n * (hs_Q2 / hs_Q1)^2)

    isom = E_co.isomorphism_to(E_co_red)

    isom_phi_P = isom(phi_P)
    isom_phi_Q = isom(phi_Q)

    return A_co_red, isom_phi_P[0], isom_phi_Q[0]

def ReducedMontCoeff(A):
    J = EllipticCurve(Fp, [0, Fp(A), 0, 1, 0]).j_invariant()
    if J in unique_node : return unique_node[J]
    else : unique_node[J] = Fp(A)
    return Fp(A)

def SubgroupToDirect(A, xP, xQ, n):
    E = EllipticCurve(Fp, [0, Fp(A), 0, 1, 0])
    P, Q = E.gens()
    P, Q = ((p+1)//n) * P, ((p+1)//n) * Q

    K = E(0)
    if xP != 0: K = E.lift_x(xP)
    else : K = E.lift_x(xQ)
    for x in range(n):
        R = P + x*Q
        if R.weil_pairing(K, n) == 1: return x

    R = Q
    if R.weil_pairing(K, n) == 1: return n

    # Unreachable
    raise Exception("Unreachable : No direction for given subgroup")

def MontInv(A):
    return EllipticCurve(Fp, [0, A, 0, 1, 0]).j_invariant()

f = open(str(p) + '_' + str(ell) + '_graph.dot', 'w')

f.write('digraph{\n')
cand_coeff = [Fp(0)]
for depth in range(BFS_depth):
    new_cand_coeff = []
    #print(cand_coeff)
    for elem in cand_coeff:
        src_A = elem
        E_src = EllipticCurve(Fp, [0, src_A, 0, 1, 0]).j_invariant()
        graph_nodes.add(str(E_src)) # For drawing
        for direct in range(0, ell + 1):
            if (E_src, direct) in edges_storage : continue
            dst_A, eval_P, eval_Q = Isog(src_A, ell, direct)            
            E_dst = EllipticCurve(Fp, [0, dst_A, 0, 1, 0]).j_invariant()
            udirect = SubgroupToDirect(dst_A, eval_P, eval_Q, ell)
            edges_storage.append((E_src, direct))
            if (E_dst, udirect) not in edges_storage:
                edges_storage.append((E_dst, udirect))
                f.write('"' + str(E_src) + '"->"' + str(E_dst) + '"' + '[label="'+ str(direct) + str(udirect) + '"]\n')
                graph_links.add((str(E_src), str(E_dst), str(direct) + str(udirect)))
            if E_dst not in total_visit:
                new_cand_coeff.append(dst_A)
                total_visit.append(E_dst)
    cand_coeff = new_cand_coeff

f.write('}')
f.close()

dataset = {}
dataset["nodes"] = []
dataset["links"] = []
for elem in list(graph_nodes):
    dataset["nodes"].append({"id":elem})

for elem in list(graph_links):
    dataset["links"].append({"source" : elem[0], "target" : elem[1], "value" : elem[2]})

with open(str(p) + '_' + str(ell) + '_graph.json','w') as f:
    f.write(json.dumps(dataset))
