from sage.libs.gap.element import GapElement
import matplotlib.pyplot as plt
import numpy as np
from time import time

def gen_combinations(g,N,l):
    if l==0:
        yield g
        return
    for gm in gen_combinations(g,N,l-1):
        for n in N:
            old = gm[l-1]
            gm[l-1] = old * n
            yield gm
            gm[l-1] = old
def Quot(G,N):
    phi = G.NaturalHomomorphismByNormalSubgroup(N)
    GbyN = phi.ImagesSource()
    return GbyN,phi
def QuotElemNSmallGensReps(G,N):
    GbyN,phi = Quot(G,N)
    elems = list(GbyN.AsList())
    gens = list(libgap.SmallGeneratingSet(GbyN))
    gens = Phinv(gens,phi)
    elems = Phinv(elems,phi)
    return GbyN,elems,gens,phi
def Phinv(elem,phi):
    if isinstance(elem,list):
        return [phi.PreImagesRepresentative(x) for x in elem]
    return phi.PreImagesRepresentative(elem)
def Phi(g,phi):
    if isinstance(g,list):
        return [phi.ImagesRepresentative(x) for x in g]
    return phi.ImagesRepresentative(g)
def is_group_by_gens(group, gens):
    G = libgap.GroupByGenerators(gens)
    return group == G
def is_Qgroup_by_reps(phi,GbyN,reps): return (GbyN == libgap.GroupByGenerators([phi.ImagesRepresentative(x) for x in reps]))
def lift(G_by_Gim1_mingen_reps,
         Gim1_by_Gi,
         G_by_Gi,
         phi_G_by_Gi,
         phi_Gim1_by_Gi,
         debug=False
         ):
    """
        G_by_Gim1 = G / G_{i-1}
        Gim1_by_Gi = G_{i-1} / G_{i}
        phi_G_by_Gi is the homomorphism defining the cosets of Gi in G. We use this to find g_j G_i
        We want to find mingen of G / G_{i} by lifting GbyGim1 with Gim1byGi
    """

    def gen_combinations(g,N,l):
        if l==0:
            yield g
            return
        for gm in gen_combinations(g,N,l-1):
            for n in N:
                old = gm[l-1]
                gm[l-1] = old * n
                yield gm
                gm[l-1] = old

    g = G_by_Gim1_mingen_reps
    if debug : print("g :",g)

    l = len(g)
    if debug : print("l :",l)

    old_G_phi = phi_G_by_Gi
    if debug : print("old_G_phi :",old_G_phi)

    old_G = G_by_Gi
    if debug : print("old_G :",old_G)

    Gim1_by_Gi_L = list(Gim1_by_Gi.AsList())
    if debug : print("Gim1_by_Gi_L :",Gim1_by_Gi_L)

    Gim1_by_Gi_elem_reps = [phi_Gim1_by_Gi.PreImagesRepresentative(x) for x in Gim1_by_Gi_L]
    if debug : print("Gim1_by_Gi_elem_reps",Gim1_by_Gi_elem_reps)

    Gim1_by_Gi_gen = list(libgap.SmallGeneratingSet(Gim1_by_Gi))
    if debug : print("Gim1_by_Gi_gen",Gim1_by_Gi_gen)

    Gim1_by_Gi_gen_reps = [phi_Gim1_by_Gi.PreImagesRepresentative(x) for x in Gim1_by_Gi_gen]
    if debug : print("Gim1_by_Gi_gen_reps",Gim1_by_Gi_gen_reps)

    N = Gim1_by_Gi
    if debug : print("N",N)

    N_list = Gim1_by_Gi_elem_reps
    if debug : print("N_list",N_list)

    n = Gim1_by_Gi_gen_reps
    if debug : print("n",n)

    if N.IsAbelian().sage():
        if (old_G == libgap.GroupByGenerators([old_G_phi.ImagesRepresentative(x) for x in g])):
            if debug : print("unmodified g works")
            return g

        for i in range(l):
            for j in range(len(n)):
                temp = g[i]
                g[i] = g[i]*n[j]
                if (old_G == libgap.GroupByGenerators([old_G_phi.ImagesRepresentative(x) for x in g])):
                    if debug : print("{g1,g2...gi*nj...gl} works")
                    return g
                g[i] = temp
        if debug: print("returning g U { n_0 }")
        return g + [n[0]]
    
    for raw_gens in gen_combinations(g, N_list, l):
        if (old_G == libgap.GroupByGenerators([old_G_phi.ImagesRepresentative(x) for x in raw_gens])):
            if debug : print("returning {g1n1,g2n2...glnl}")
            return raw_gens

    for raw_gens in gen_combinations(g+[N_list[0]], N_list, l+1):
        if (old_G == libgap.GroupByGenerators([old_G_phi.ImagesRepresentative(x) for x in raw_gens])):
            if debug : print("returning {g1n1,g2n2...glnl,n_{l+1}}")
            return raw_gens
    if debug : assert False ,("This stage shouldn't ever be reached.")
def LIFT(cs,j,i,mingenset_j_reps=None,debug=False):
    """
    'cs' is the chief series [G0,G1,G2...]

    Given the representatives of minimum generating set for G/Gj ,
    this function finds the representatives of the minimum generating set of G/Gi
    where i > j
    """
    assert len(cs)>i and i>j and j>0
    G = cs[0]
    Gj = cs[j]
    Gi = cs[i]
    if debug:
        print("G :",G )
        print("Gi:",Gi)
        print("Gj:",Gj)
    phi_GbyGj = G.NaturalHomomorphismByNormalSubgroup(Gj)
    GbyGj = phi_GbyGj.ImagesSource()
    if debug:
        print("GbyGj    ",GbyGj    )
        print("phi_GbyGj",phi_GbyGj)
    phi_GbyGi = G.NaturalHomomorphismByNormalSubgroup(Gi)
    GbyGi = phi_GbyGi.ImagesSource()
    if debug:
        print("GbyGi    ",GbyGi    )
        print("phi_GbyGi",phi_GbyGi)
    if mingenset_j_reps is None : 
        mingenset_j = list(libgap.SmallGeneratingSet(GbyGj))
        mingenset_j_reps = [phi_GbyGj.PreImagesRepresentative(x) for x in mingenset_j]
    mingenset_k_reps = mingenset_j_reps
    for k in range(j+1,i+1): 
        mingenset_km1_reps = mingenset_k_reps
        Gk = cs[k]
        Gkm1 = cs[k-1]
        phi_GbyGk = G.NaturalHomomorphismByNormalSubgroup(Gk)
        GbyGk = phi_GbyGk.ImagesSource()
        phi_Gkm1byGk = Gkm1.NaturalHomomorphismByNormalSubgroup(Gk)
        Gkm1byGk = phi_Gkm1byGk.ImagesSource()
        mingenset_k_reps = lift(
            mingenset_km1_reps,
            Gkm1byGk,
            GbyGk,
            phi_GbyGk,
            phi_Gkm1byGk
        )
    assert (GbyGi == libgap.GroupByGenerators([phi_GbyGi.ImagesRepresentative(x) for x in mingenset_k_reps]))
    return mingenset_k_reps
def LIFT_TEST():
    D = DihedralGroup(3).gap()
    D2 = D.DirectProduct(D)
    G = D2.DirectProduct(D)
    print("G :",G)
    cs = G.ChiefSeries()
    l = len(cs)
    print("l :",l)
    gens = LIFT(cs,1,l-1)
    print("gens :",gens)
def lift_test():
    G = DihedralGroup(4).gap()
    print("G :",G)
    cs = G.ChiefSeries()
    print("Chief series length :",len(cs))
    G1byG2,G1byG2_elems,G1byG2_gens,G1byG2_phi = QuotElemNSmallGensReps(cs[1],cs[2])
    print("G1byG2       : ",G1byG2      )
    print("G1byG2_elems : ",G1byG2_elems)
    print("G1byG2_gens  : ",G1byG2_gens )
    print("G1byG2_phi   : ",G1byG2_phi  )
    
    G0byG1,G0byG1_elems,G0byG1_gens,G0byG1_phi = QuotElemNSmallGensReps(cs[0],cs[1])
    print("G0byG1       : ",G0byG1      )
    print("G0byG1_elems : ",G0byG1_elems)
    print("G0byG1_gens  : ",G0byG1_gens ) 
    print("G0byG1_phi   : ",G0byG1_phi  ) 
    G0byG2,G0byG2_phi = Quot(cs[0],cs[2])
    print("G0byG2     : ",G0byG2    )
    print("G0byG2_phi : ",G0byG2_phi)
    gens = lift(
        G0byG1_gens,
        G1byG2,
        G0byG2,
        G0byG2_phi,
        G1byG2_phi,
    )
    print("gens :",gens)
    print("Works :",is_Qgroup_by_reps(G0byG2_phi,G0byG2,gens))
def test_mingen():
    A5 = AlternatingGroup(5).gap()
    A5_2 = A5.DirectProduct(A5)
    Z5 = PermutationGroup([(1,2,3,4,5)]).gap()
    P = PermutationGroup([(1,2,3), (2,3), (4,5)]).gap()
    L = [
        (A5,2),
        (A5_2,2),
        (Z5,1),
        (P,2),
        ]
    for (G,l) in L:
        gens = minimum_generating_set(G)
        print("Group : ",G)
        print("mingen :",gens)
        if not is_group_by_gens(G,gens):
            print("fail\n")
            break
        if not len(gens) == l:
            print("fail\n")
            break
        print("Pass\n")

def Z_p_S_3(p):
    S = SymmetricGroup(3).gap()
    Z = libgap.CyclicGroup(p)
    SZ = S.DirectProduct(Z)
    return SZ
def Z_2_to_n(n):
    Z_2 =  PermutationGroup([(1,2)]).gap()
    G = Z_2
    for i in range(1,n):
        G = G.DirectProduct(Z_2)
    return G
def Z_n(n):
    Z = PermutationGroup([tuple((i+1 for i in range(n)))]).gap()
    #Z =  libgap.CyclicGroup(n)
    return Z
def S_n(n):
    S = SymmetricGroup(n).gap()
    return S
def D_n(n):
    D = DihedralGroup(n).gap()
    return D
def H_n_2(n):
    H = groups.matrix.Heisenberg(n,2)
    H = H.gap()
    return H
def A_5_to_n(n):
    A5 = AlternatingGroup(5).gap()
    G = A5
    for i in range(n-1):
        G = G.DirectProduct(A5)
    return G

def minimum_generating_set(G,debug=False,report_time=False)->list:
    assert isinstance(G, GapElement)
    if not G.IsFinite().sage(): raise NotImplementedError("only implemented for finite groups")
    try:return list(libgap.MinimalGeneratingSet(G))
    except:pass
    cs = G.ChiefSeries()
    l = len(cs)-1
    if report_time : t0 = time()
    #gens = LIFT(cs,1,l)
    i = l
    j = 1
    """
    LIFT
    'cs' is the chief series [G0,G1,G2...]

    Given the representatives of minimum generating set for G/Gj ,
    this function finds the representatives of the minimum generating set of G/Gi
    where i > j
    """

    def lift(G_by_Gim1_mingen_reps,
             Gim1_by_Gi,
             G_by_Gi,
             phi_G_by_Gi,
             phi_Gim1_by_Gi,
             debug=False
             ):
        """
            G_by_Gim1 = G / G_{i-1}
            Gim1_by_Gi = G_{i-1} / G_{i}
            phi_G_by_Gi is the homomorphism defining the cosets of Gi in G. We use this to find g_j G_i
            We want to find mingen of G / G_{i} by lifting GbyGim1 with Gim1byGi
        """

        def gen_combinations(g,N,l):
            if l==0:
                yield g
                return
            for gm in gen_combinations(g,N,l-1):
                for n in N:
                    old = gm[l-1]
                    gm[l-1] = old * n
                    yield gm
                    gm[l-1] = old

        g = G_by_Gim1_mingen_reps
        if debug : print("g :",g)

        l = len(g)
        if debug : print("l :",l)

        old_G_phi = phi_G_by_Gi
        if debug : print("old_G_phi :",old_G_phi)

        old_G = G_by_Gi
        if debug : print("old_G :",old_G)

        Gim1_by_Gi_L = list(Gim1_by_Gi.AsList())
        if debug : print("Gim1_by_Gi_L :",Gim1_by_Gi_L)

        Gim1_by_Gi_elem_reps = [phi_Gim1_by_Gi.PreImagesRepresentative(x) for x in Gim1_by_Gi_L]
        if debug : print("Gim1_by_Gi_elem_reps",Gim1_by_Gi_elem_reps)

        Gim1_by_Gi_gen = list(libgap.SmallGeneratingSet(Gim1_by_Gi))
        if debug : print("Gim1_by_Gi_gen",Gim1_by_Gi_gen)

        Gim1_by_Gi_gen_reps = [phi_Gim1_by_Gi.PreImagesRepresentative(x) for x in Gim1_by_Gi_gen]
        if debug : print("Gim1_by_Gi_gen_reps",Gim1_by_Gi_gen_reps)

        N = Gim1_by_Gi
        if debug : print("N",N)

        N_list = Gim1_by_Gi_elem_reps
        if debug : print("N_list",N_list)

        n = Gim1_by_Gi_gen_reps
        if debug : print("n",n)

        if N.IsAbelian().sage():
            if (old_G == libgap.GroupByGenerators([old_G_phi.ImagesRepresentative(x) for x in g])):
                if debug : print("unmodified g works")
                return g

            for i in range(l):
                for j in range(len(n)):
                    temp = g[i]
                    g[i] = g[i]*n[j]
                    if (old_G == libgap.GroupByGenerators([old_G_phi.ImagesRepresentative(x) for x in g])):
                        if debug : print("{g1,g2...gi*nj...gl} works")
                        return g
                    g[i] = temp
            if debug: print("returning g U { n_0 }")
            return g + [n[0]]

        for raw_gens in gen_combinations(g, N_list, l):
            if (old_G == libgap.GroupByGenerators([old_G_phi.ImagesRepresentative(x) for x in raw_gens])):
                if debug : print("returning {g1n1,g2n2...glnl}")
                return raw_gens

        for raw_gens in gen_combinations(g+[N_list[0]], N_list, l+1):
            if (old_G == libgap.GroupByGenerators([old_G_phi.ImagesRepresentative(x) for x in raw_gens])):
                if debug : print("returning {g1n1,g2n2...glnl,n_{l+1}}")
                return raw_gens
        if debug : assert False ,("This stage shouldn't ever be reached.")

    Gj = cs[j]
    Gi = cs[i]
    if debug:
        print("G :",G )
        print("Gi:",Gi)
        print("Gj:",Gj)
    phi_GbyGj = G.NaturalHomomorphismByNormalSubgroup(Gj)
    GbyGj = phi_GbyGj.ImagesSource()
    if report_time :
        t = time()
        dt = t- t0
        t0 = t
        print("GbyG1:",dt)
    if debug:
        print("GbyGj    ",GbyGj    )
        print("phi_GbyGj",phi_GbyGj)
    
    phi_GbyGi = G.NaturalHomomorphismByNormalSubgroup(Gi)
    GbyGi = phi_GbyGi.ImagesSource()
    if report_time :
        t = time()
        dt = t- t0
        t0 = t
        print("GbyGl:",dt)
    if debug:
        print("GbyGi    ",GbyGi    )
        print("phi_GbyGi",phi_GbyGi)
    mingenset_j = list(libgap.SmallGeneratingSet(GbyGj))
    mingenset_j_reps = [phi_GbyGj.PreImagesRepresentative(x) for x in mingenset_j]
    if report_time :
        t = time()
        dt = t- t0
        t0 = t
        print("mingenset_j_reps:",dt)
    mingenset_k_reps = mingenset_j_reps
    for k in range(j+1,i+1): 
        mingenset_km1_reps = mingenset_k_reps
        Gk = cs[k]
        Gkm1 = cs[k-1]
        phi_GbyGk = G.NaturalHomomorphismByNormalSubgroup(Gk)
        GbyGk = phi_GbyGk.ImagesSource()
        if report_time :
            t = time()
            dt = t- t0
            t0 = t
            print(f"GbyG{k}:",dt)
        phi_Gkm1byGk = Gkm1.NaturalHomomorphismByNormalSubgroup(Gk)
        Gkm1byGk = phi_Gkm1byGk.ImagesSource()
        if report_time :
            t = time()
            dt = t- t0
            t0 = t
            print(f"G{k-1}byG{k}:",dt)
        mingenset_k_reps = lift(
            mingenset_km1_reps,
            Gkm1byGk,
            GbyGk,
            phi_GbyGk,
            phi_Gkm1byGk
        )
        if report_time :
            t = time()
            dt = t- t0
            t0 = t
            print(f"lift to get mingenset_{k}_reps:",dt)
    assert (GbyGi == libgap.GroupByGenerators([phi_GbyGi.ImagesRepresentative(x) for x in mingenset_k_reps]))

    gens = mingenset_k_reps

    assert (G == libgap.GroupByGenerators(gens))
    return gens 
MGS = minimum_generating_set

def TAP():
    ToWrite = r"\begin{matrix}"
    D = {
        H_n_2:(9,2,1,1,r"H(n,2)"),
        Z_p_S_3:(60,1,3,3,r"Z_n\times S_3"),
        Z_2_to_n:(20,2,2,5,"Z_2^n"),
        S_n:(10,2,1,2,'S_n'),
        Z_n:(500,2,10,5,'Z_n'),
        D_n:(100,1,1,5,'D_n'),
        A_5_to_n:(10,1,1,2,'A_5^n'),
    }
    ToWrite += " & ".join(["Group Type",r"\text{len}(G)",r"\text{len}(\text{mingen}(G))",r"\text{mingen}(g)"]) + r'\\'+'\n'
    I = 0
    for Gfunc in D:
        N,N0,d,iterations,name = D[Gfunc]
        I += 1
        plt.figure()

        y = []
        x = []
        for n in range(N0,N,d):
            G = Gfunc(n)
            assert G is not None , Gfunc
            to = time()
            for _ in range(iterations):
                g = minimum_generating_set(G)
            gs = ",".join([str(x) for x in g])
            if len(gs) > 100: gs = r"\text{too long to write}"
            ToWrite += ' & '.join([name,str(G.Size().sage()),str(len(g)),gs]) + r'\\'+'\n'
            y.append((time()-to)/iterations)
            x.append(int(G.Size().sage()))
        y = np.log2(np.array(y))
        x = np.array(x)
        plt.plot(x,y,'-o')

        # Curve fitting
        Ln = np.log2(x)
        X = np.array([[ln,1] for ln in Ln])
        XT = X.T
        XTX = X.T @ X
        XTXi = np.linalg.inv(XTX)
        pseudo_inverse = XTXi @ XT
        theta = np.dot(pseudo_inverse,y)
        a,b = theta
        x = np.linspace(x[0],x[-1],1000)
        Ln = np.log2(x)
        L_pred = Ln*a + b
        plt.plot(x,L_pred)

        plt.xlabel("$ |G| $")
        plt.ylabel(r"$\log_2(t)$ ($ t $ in seconds)")
        plt.title(f"Time ($ t $) to find minimum generating set for $ G = {name} $")
        plt.legend(["Actual",f"{round(a,2)}log_2|G| + ({round(b,2)})"])
        plt.savefig(f"/home/hp/Figure {I}")
        plt.close()
        print(f"{name} done")
    ToWrite += r"\end{matrix}"
    file = open("/home/hp/OutPut.txt",'w')
    file.write(ToWrite)
def mingtest():
    A5 = AlternatingGroup(5).gap()
    G = A5
    n = 9
    for i in range(n):
        print(f"A_5^{i+1} 's mingen :")
        g = MGS(G)
        print(g)
        H = libgap.GroupByGenerators(g)
        assert H == G
        print(len(g))
        if i < n-1:
            G = G.DirectProduct(A5)
