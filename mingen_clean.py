from sage.libs.gap.element import GapElement
def minimum_generating_set(G)->list:
    r"""
    INPUT:

        - ``G`` -- The group whose minimum generating set we want. This must be converted to a 'gap based' group first if it's not via ``G=G.gap()`` .
    
    OUTPUT:

        - minimum generating set of ``G``, which is the the set `g` of elements of `G` of smallest cardinality such that `G` is generated by `g`, i.e. `\braket{g}=G`
    
    ALGORITHM:

    First we cover the cases when the Chief series is of length 1, that is if `G` is simple. This case is handled by ``libgap.MinimalGeneratingSet`` , so we assume that ``libgap.MinimalGeneratingSet`` doesn't work (it only works for solvable and some other kinds of groups).
    So, we are guaranteed to find a chief series of length at least 2, since `G` is not simple. Then we proceed as follows..

    `S := ChiefSeries(G)`

    `l := len(S)-1` (this is index of last normal subgroup, namely `G_l=\{e\}`)

    Let `g` be the set of representatives of the minimum generating set of `G/S[1]` . (This can be found using ``libgap.MinimalGeneratingSet`` since `G/S[1]` is simple group)

    for k = 2 to 'l':
        
        Compute ``GbyGk`` := `G/S[k]`

        Compute ``GbyGkm1`` := `G/S[k-1]`

        Compute ``Gkm1byGk`` := `S[k-1]/S[k]`

        `g := lift(g,Gkm1byGk,Gkm1byGk,GbyGk)` . The lift function is discussed in the code.

    return `g`

    TESTS:

    1. `A_5^7`
    ```
    sage: def A_5_to_n(n):
    ....:     A5 = AlternatingGroup(5).gap()
    ....:     G = A5
    ....:     for i in range(n-1):
    ....:         G = G.DirectProduct(A5)
    ....:     return G
    ....: 
    sage: G = A_5_to_n(7)
    sage: g = minimum_generating_set(G)
    sage: g
    [(1,5,4,3,2)(8,9,10)(12,14,13)(17,19,18)(22,24,23)(27,29,28)(32,34,33),
     (2,4,3)(6,8,10,7,9)(11,14,15,13,12)(16,19)(18,20)(21,25,24,22,23)(26,30,29,28,27)(31,35,34)]
    sage: %timeit g = minimum_generating_set(G)
    900 ms ± 90.5 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
    ```
    2. `Z_3^10`
    ```
    sage: def Z_2_to_n(n):
    ....:     Z_2 =  PermutationGroup([(1,2)]).gap()
    ....:     G = Z_2
    ....:     for i in range(1,n):
    ....:         G = G.DirectProduct(Z_2)
    ....:     return G
    ....: 
    sage: G = Z_2_to_n(10)
    sage: g = minimum_generating_set(G)
    sage: len(g)
    10
    sage: g
    [(19,20),
     (17,18),
     (15,16),
     (13,14),
     (11,12),
     (9,10),
     (7,8),
     (5,6),
     (3,4),
     (1,2)]
    ```
    """
    assert isinstance(G, GapElement)
    if not G.IsFinite().sage(): raise NotImplementedError("only implemented for finite groups")
    try:return list(libgap.MinimalGeneratingSet(G))
    except:pass
    def lift(G_by_Gim1_mingen_reps , Gim1_by_Gi , G_by_Gi , phi_G_by_Gi , phi_Gim1_by_Gi):
        r"""
        Note that G_k is the same as S[k] in the ``minimum_generating_set`` algorithm and cs[k] in its code. 

        INPUT:

            - ``G_by_Gim1_mingen_reps`` -- representative elements of the minimum generating set of `G/G_{i-1}`

            - ``G_by_Gim1`` --  `G / G_{i-1}`  Quotient Group

            - ``Gim1_by_Gi`` --  `G_{i-1} / G_{i}` Quotinet Group

            - ``phi_G_by_Gi`` -- the homomorphism defining the cosets of `G_i` in `G`.

            - ``phi_Gim1_by_Gi`` -- the homomorphism defining the cosets of `G_{i}` in `G_{i-1}`.

        OUTPUT:
        
            representative elements of the minimum generating set of `G / G_{i}`.

        ALGORITHM:

            The inputs:
            
                `g = \{g_1,g_2,\dots g_l\}` is the set of representatives of the (supposed) minimum generating st of `G/G_{i-1}`, what we are calling ``G_by_Gim1_mingen_reps`` in the co

                `\bold{n} =\{n_1,n_2\dots n_k\}` where `\{n_1 G_i,n_2G_i \dots n_kG_{i}\}` is any generating set of `G_{i-1}/G_i
            
                `\bold{N} = \{N_1,N_2\dots N_m\}` where `G_{i-1}/G_i = \{N_1G_i,N_2G_2\dots N_m G_m\}` . 

            We wish to find the representatives of minimum generating set of `G/G_i` , which can be done as follows:

            if `G_{i-1}/G_i` is abelian :
            
                if `\braket{gG_i}= G/G_i` :

                    return `g`

            for `1 \le p \le l`  and `n_j \in \bold{n}` :
            
                `g^* = \{g_1,g_2\dots g_{p-1} ,g_p n_j,g_{p_1},\dots\}`

                if `\braket{g^* G_i} = G/G_i` :
            
                    return `g^*` 

            else:

                for any (not necessarily distinct) elements `N_{i_1},N_{i_2}\dots N_{i_t} \in \bold{N}` : 
                
                    `g^* = \{g_1N_{i_1},g_{i_2}N_{i_3}\dots g_{i_t}N_t,g_{t+1}\dots g_l\}`

                    if `\braket{g^*G_i}\; = G/G_i`:
                    
                        return `{g^*}`

                for any (not necessarily distinct) elements `N_{i_1},N_{i_2}\dots N_{i_t} N_{i_{t+1}} \in \bold{N}` : 

                    `g^* = \{g_1N_{i_1},g_{i_2}N_{i_3}\dots g_{i_t}N_t,g_{t+1}\dots g_l\}`

                    if `\braket{g^*G_i}\; = G/G_i`:
                
                        return `{g^*}`

            By now, we must have exhausted our search.

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
        l = len(G_by_Gim1_mingen_reps)
        Gim1_by_Gi_L = list(Gim1_by_Gi.AsList())
        Gim1_by_Gi_elem_reps = [phi_Gim1_by_Gi.PreImagesRepresentative(x) for x in Gim1_by_Gi_L]
        Gim1_by_Gi_gen = list(libgap.SmallGeneratingSet(Gim1_by_Gi))
        Gim1_by_Gi_gen_reps = [phi_Gim1_by_Gi.PreImagesRepresentative(x) for x in Gim1_by_Gi_gen]
        if Gim1_by_Gi.IsAbelian().sage():
            if (G_by_Gi == libgap.GroupByGenerators([phi_G_by_Gi.ImagesRepresentative(x) for x in G_by_Gim1_mingen_reps])): return G_by_Gim1_mingen_reps
            for i in range(l):
                for j in range(len(Gim1_by_Gi_gen_reps)):
                    temp = G_by_Gim1_mingen_reps[i]
                    G_by_Gim1_mingen_reps[i] = G_by_Gim1_mingen_reps[i]*Gim1_by_Gi_gen_reps[j]
                    if (G_by_Gi == libgap.GroupByGenerators([phi_G_by_Gi.ImagesRepresentative(x) for x in G_by_Gim1_mingen_reps])): return G_by_Gim1_mingen_reps
                    G_by_Gim1_mingen_reps[i] = temp
            return G_by_Gim1_mingen_reps + [Gim1_by_Gi_gen_reps[0]]
        for raw_gens in gen_combinations(G_by_Gim1_mingen_reps, Gim1_by_Gi_elem_reps, l):
            if (G_by_Gi == libgap.GroupByGenerators([phi_G_by_Gi.ImagesRepresentative(x) for x in raw_gens])): return raw_gens
        for raw_gens in gen_combinations(G_by_Gim1_mingen_reps+[Gim1_by_Gi_elem_reps[0]], Gim1_by_Gi_elem_reps, l+1):
            if (G_by_Gi == libgap.GroupByGenerators([phi_G_by_Gi.ImagesRepresentative(x) for x in raw_gens])): return raw_gens
    cs = G.ChiefSeries()
    l = len(cs)-1
    phi_GbyG1 = G.NaturalHomomorphismByNormalSubgroup(cs[1])
    GbyG1 = phi_GbyG1.ImagesSource()
    mingenset_k_reps = [phi_GbyG1.PreImagesRepresentative(x) for x in list(libgap.SmallGeneratingSet(GbyG1))] # k=1 initially
    for k in range(2,l+1): 
        mingenset_km1_reps = mingenset_k_reps
        Gk,Gkm1 = cs[k], cs[k-1]
        phi_GbyGk = G.NaturalHomomorphismByNormalSubgroup(Gk)
        GbyGk = phi_GbyGk.ImagesSource()
        phi_Gkm1byGk = Gkm1.NaturalHomomorphismByNormalSubgroup(Gk)
        Gkm1byGk = phi_Gkm1byGk.ImagesSource()
        mingenset_k_reps = lift(mingenset_km1_reps,Gkm1byGk,GbyGk,phi_GbyGk,phi_Gkm1byGk)
    assert (G == libgap.GroupByGenerators(mingenset_k_reps))
    return mingenset_k_reps 
