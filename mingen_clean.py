from sage.libs.gap.element import GapElement
def minimum_generating_set(G)->list:
    assert isinstance(G, GapElement)
    if not G.IsFinite().sage(): raise NotImplementedError("only implemented for finite groups")
    try:return list(libgap.MinimalGeneratingSet(G)) 
    except:pass
    def lift(G_by_Gim1_mingen_reps , Gim1_by_Gi , G_by_Gi , phi_G_by_Gi , phi_Gim1_by_Gi):
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