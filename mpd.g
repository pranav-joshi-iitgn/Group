SporadicGroupMinimalFaithfulPermutationDegrees :=  rec(
    M11 := 11,
    M12 := 12,
    M22 := 22,
    M23 := 23,
    M24 := 24,
    Co1 := 98280,
    Co2 := 2300,
    Co3 := 276,
    McL := 275,
    HS := 100,
    Suz := 1782,
    Fi22 := 3510,
    Fi23 := 31671,
    Fi24 := 306936,
    M := 97239461142009186000,
    B := 13571955000,
    Th := 143127000,
    HN := 1140000,
    He := 2058,
    J1 := 266,
    J2 := 100,
    J3 := 6156,
    J4 := 173067389,
    ON := 122760,
    Ly := 8835156,
    Ru := 4060,
);

MinimalFaithfulPermutationDegreeOfSimpleGroupWithIsomorphismType := function (info)
    # This function is derived from table 4 of this paper :
    # https://www.ams.org/journals/tran/2015-367-11/S0002-9947-2015-06293-X/S0002-9947-2015-06293-X.pdf

    # `info` is the a record containing information about the type of the simple group,
    # like one obtained from `IsomorphismTypeInfoFiniteSimpleGroup`. You need to give the series,
    # parameter, and shortname values in the record.
    local
        series,       # series of simple groups
        d,            # mostly the dimension of vector space for classical groups
        q,            # elements in the field over which group is defined
        m,            # first parameter
        b;            # q = p^b for some prime p

    series := info.series;
    if series = "Spor" then
        return SporadicGroupMinimalFaithfulPermutationDegrees.(info.shortname);
    else
        if IsList(info.parameter) then
            q := info.parameter[2];
            m := info.parameter[1];
        else
          q := info.parameter;
        fi;
    fi;
    if series = "Z" then #Cyclic group of prime order
        return q;
    elif series = "A" then # Alternating group
        return q;
    elif series = "L" then # PSL(m,q)
        d := m;
        if (d = 4 and q = 2) then return 8; fi;
        if d = 2 then
            if q = 9 then return 6; fi;
            if q in [5,7,11] then return q; fi;
        fi;
        return (q^d-1)/(q-1);
    elif series = "2A" then # PSU(m+1,q)
        d := m + 1;
        if d = 3 and q = 5 then return 50; fi;
        if d = 3 then return q^3 + 1; fi;
        if d = 4 then return (q+1)*(q^3 + 1); fi;
        if d mod 2 = 0 and q = 2 then return (2^(d-1)*(2^d - 1))/3; fi;
        return ((q^d - (-1)^d)*(q^(d-1) + (-1)^d))/(q^2-1);
    elif series = "B" then # P\Omega(2*m+1,q) or O
        if q = 3 and m > 2 then return (3^m)*(3^m - 1)/2;
        elif q > 4 and m > 2 then return (q^(2*m)-1)/(q-1);
        elif q = 3 and m = 2 then return 27; #Special case : B(2,3) ~ 2A(3,2) = PSU(4,2)
        elif q = 2 and m = 2 then # B(2,2) is not a simple group.
          Error("B(2,2) is not a simple group. This shouldn't be happening.\n");
        elif m = 2 then return (q^(2*m) -1)/(q-1); #Special case : B(2,q) ~ C(2,q) = PSp(4,q)
        elif q = 2 then return (2^(m-1))*(2^m -1); #Special case : B(m,2) ~ C(m,2) = PSp(2*m,2)
        else Error("series B and m,q not of proper form\n"); fi;
    elif series = "2B" then # Sz or _2 B^2
        b := Log2Int(q);
        if 2^b = q and b mod 2 = 1 then return q^2 + 1;
        else Error("2B series without q of proper form\n"); fi;
    elif series = "C" then # PSp(2*m,q)
        d := 2*m;
        if d=4 and q=2 then return 6;fi;
        if d=4 and q=3 then return 27;fi;
        if m>2 and q=2 then return (2^(m-1))*(2^m -1);fi;
        if m>1 and q>2 then return (q^d -1)/(q-1);fi;
        Error("series C and m,q are not of proper form\n");
    elif series = "D" then # POmega(+1,2*m,q) or O+
        if m > 3 then
            if q < 4 then return q^(m-1)*(q^m -1)/(q-1);
            else return (q^(m-1) + 1)*(q^m-1)/(q-1); fi;
        else Error("series D and m < 4\n"); fi;
    elif series = "2D" then # POmega(-1,2*m,q) or O-
        if m > 3 then return (q^m+1)*(q^(m-1)-1)/(q-1);
        else Error("series 2D and m < 4\n"); fi;
    elif series = "3D" then # ^3 D_4
        return (q^8 + q^4 + 1)*(q+1);
    elif series = "E" then #E_n(q)
        if m = 6 then return (q^9 - 1)*(q^8 + q^4 + 1)/(q-1);
        elif m = 7 then return (q^14 - 1)*(q^9 + 1)*(q^5 -1)/(q-1);
        elif m = 8 then return (q^30 - 1)*(q^12 + 1)*(q^10 + 1)*(q^6 + 1)/(q-1);
        else Error("series E and m is not 6,7 or 8\n");
        fi;
    elif series = "2E" then #2E(6,q)
        return (q^12 -1)*(q^6 - q^3 + 1)^(q^4 +1)/(q-1);
    elif series = "F" then #F(4,q)
        return (q^12 -1)*(q^4 + 1)/(q-1);
    elif series = "2F" then #2F(4,q)
        if q = 2 then return 1600; #special case : 2F4(2) ~ Tits
        else return (q^6 + 1)*(q^3 + 1)*(q+1); fi;
    elif series = "G" then #G(2,q)
        if q = 3 then return 351;
        elif q = 4 then return 416;
        else return (q^6 -1)/(q-1); fi;
    elif series = "2G" then #2G(2,q)
        b := PValuation(q,3);
        if q = 3^b and b mod 2 = 1 then return q^3 + 1;
        else Error("series 2G and q not of proper form\n"); fi;
    fi;
    Error("series `",series,"` is not valid\n");
end;

MinimalFaithfulPermutationDegreeOfSimpleGroup := function(G)
    local info;
    info := IsomorphismTypeInfoFiniteSimpleGroup(G); #This requires computing Size(G)
    return MinimalFaithfulPermutationDegreeOfSimpleGroupWithIsomorphismType(info);
end;

#########################################################################
##
#F  CanLiftPermutationDegreeToNormaliserOfSimpleMatrixGroup(<S>,<G>,<d>,<q>,<type>)
##
##  Checks whether mu(A) = mu(S) when S is isomorphic to the
##  simple matrix group specified by "d","q" and "type",
##  where A is an automorphism subgroup of S,
##  made of all conjugators ^n for all n in normaliser of S in G,
##  and mu(X) is the Minimal Faithful Permutation Degree of X.
##  This is implemented, as described in this research paper:
##  https://dl.acm.org/doi/10.1145/3618260.3649641
CanLiftPermutationDegreeToNormaliserOfSimpleMatrixGroup := function(S,G,d,q,type)
    # S is a simple subgroup of G, which is isomorphic to PSL(d,q) , or PSp(d,q) , or POmega(1,2d,q)
    # type is the type of simple group we are dealing with ("PSL","PSp","POmega+")
    # For the sake of sanity, the varaibles are chosen to represent the process for type="PSL" only,
    # but work functionally, for every type of simple matrix group specified.
    local
        Sldq,           # SL(d,q) with generators Slgen
        Slgen,          # Very specific generating set of SL(d,q), or correspoding group based on "type"
        p,              # q = p^e
        s,g,            # Any element in S
        m2,             # A homomorphism between SL(d,q) and PSL(d,q)
        i,j,x,r,        # Dummy variables
        pS,             # PSL(d,q)
        AutpSgens,      # Generators of AutpS
        lambda,         # Any automorphism from Sldq to Sldq
        C,              # Centraliser of S in G,
                        # and later, the Center of Sldq,
        m2lambm2inv,    # A function that returns the composition m2inv * lambda * m2 ,
                        # given lambda, where m2inv maps to the elements in coset with order p
        AutSlgens,      # Generators of AutSl
        alpha,          # The automorphism corresponding to m2lambm2inv
        alphaSlgen,     # Images of Slgen under alpha
        StopS,          # An isomorphism from S to pS
        pStoS,          # An isomorphism from pS to S
        A,              # Group of conjugation automorphisms of S with elements from G
        U,              # Any matrix in Slgen
        N,              # Normaliser of S in G
        NbyC,           # Natural Homomorphism from N, by normal subgroup C
        Ngens,          # Generators of N
        AutS,           # The group made of conjugator automorphisms ^n for all n in N,
                        # also called A.
        NgensNeeded,    # A minimal subset of Ngens that still creates AutS
        good,           # A boolean
        n,              # any element of N
        AutSgens,       # Generators of group AutS (also called A)
        pSgens,         # Generators of pS
        StoS,           # Any automorphism of S, can also be called lambda
        alphapSgens,    # Images of Hgens under alpha
        h,              # any element of pS
        pStopS,         # Any automorphism of pS, can also be called alpha
        lambda_g,       # Image of g under lambda (StoS)
        alpha_h,        # Image of h under alpha (pStoHpS)
        AutpS,          # The Automorphism subgroup of pS isomorphic to A
        g_list,         # preimages of pSgens under the isomorphism StopS
        lambda_g_list,  # Images of g_list under any isomorphism of S, StoS
        L,              # Any list
        conj,           # true if everything in A is a conjugator,
                        # upto a field automorphism
        M,              # Any matrix
        f,              # Any field automorphism of GF(q)
        alphafinv,      # alpha * f^-1
        scalars,        # all required elements from GF(q)
        flisinit,       # initial list of field automorphisms
        finvlisinit,    # f^-1 for f in flisinit
        findlis,        # list of indices of flisinit
        finvlis,finv,
        FirstOfimagesOfscalarMatrices,
        scalarMatrices,imagesOfscalars,
        MakeFieldAutomorphismOfSldq,finvmatrixlist,O,phi,PO,POgens;

    # Computing Normaliser and Centraliser
    N := Normaliser(G,S);
    C := Centraliser(G,S);
    Ngens := SmallGeneratingSet(N);

    # Computing NgensNeeded
    NbyC := NaturalHomomorphismByNormalSubgroupNC(N,C);
    NgensNeeded := ImagesSet(NbyC,Ngens);
    NgensNeeded := List(NgensNeeded,x->PreImagesRepresentative(NbyC,x));

    # Finding generators of AutS (A)
    AutSgens := List(NgensNeeded,n -> ConjugatorAutomorphismNC(S,n));

    # Calculating Sldq
    p := PrimeDivisors(q)[1];
    Slgen := [];
    if type = "PSL" then
        for i in [1..d] do
            for j in [1..d] do
                if i=j then continue; fi;
                for r in [0,1..p-1] do
                    U := IdentityMatrix(GF(q),d);
                    U[i][j] := Z(q)^r;
                    Add(Slgen,U);
                od;
            od;
        od;
    elif type="POmega+" then
        for j in [2..d] do
            for i in [1..j-1] do
                for r in [0,1..q-2] do
                    # I + z(e_{i,j} - e_{-j,-i})
                    U := IdentityMatrix(GF(q),2*d);
                    U[i][j] := Z(q)^r;
                    U[d+j][d+i] := -Z(q)^r;
                    Add(Slgen,U);
                    # I + z(e_{j,i} - e_{-i,-j})
                    U := IdentityMatrix(GF(q),2*d);
                    U[j][i] := Z(q)^r;
                    U[d+i][d+j] := -Z(q)^r;
                    Add(Slgen,U);
                    # I + z(e_{i,-i} - e_{j,-i})
                    U := IdentityMatrix(GF(q),2*d);
                    U[i][d+j] := Z(q)^r;
                    U[j][d+i] := -Z(q)^r;
                    Add(Slgen,U);
                    # I + z(e_{-i,j} - e_{-j,i})
                    U := IdentityMatrix(GF(q),2*d);
                    U[d+i][j] := Z(q)^r;
                    U[d+j][i] := -Z(q)^r;
                    Add(Slgen,U);
                od;
            od;
        od;
        d := 2*d;
    elif type="PSp" then
        for j in [2..d] do
            for i in [1..j-1] do
                for r in [0,1..q-2] do
                    # I + z(e_{i,j} - e_{-j,-i})
                    U := IdentityMatrix(GF(q),2*d);
                    U[i][j] := Z(q)^r;
                    U[d+j][d+i] := -Z(q)^r;
                    Add(Slgen,U);
                    # I - z(e_{-i,-j} - e_{j,i})
                    U := IdentityMatrix(GF(q),2*d);
                    U[j][i] := Z(q)^r;
                    U[d+i][d+j] := -Z(q)^r;
                    Add(Slgen,U);
                    # I + z(e_{i,-i} + e_{j,-i})
                    U := IdentityMatrix(GF(q),2*d);
                    U[i][d+j] := Z(q)^r;
                    U[j][d+i] := Z(q)^r;
                    Add(Slgen,U);
                    # I - z( - e_{-i,j} - e_{-j,i})
                    U := IdentityMatrix(GF(q),2*d);
                    U[d+i][j] := Z(q)^r;
                    U[d+j][i] := Z(q)^r;
                    Add(Slgen,U);
                od;
            od;
        od;
        for i in [1..d] do
            # I + z e_{i,-i}
            U := IdentityMatrix(GF(q),2*d);
            U[i][d+i] := Z(q)^r;
            Add(Slgen,U);
            # I + z e_{-i,i}
            U := IdentityMatrix(GF(q),2*d);
            U[d+i][i] := Z(q)^r;
            Add(Slgen,U);
        od;
        d := 2*d;
    fi;
    Sldq := GroupWithGenerators(Slgen);

    # Calculating m2 and pS
    if type="PSL" or type="PSp" then
        C := Center(Sldq);
    elif type="POmega+" then
        C := Center(GO(1,d,q));
        C := Filtered(C,x->x in Sldq); #There are only O(q) elements in C
        C := AsGroup(C);
    fi;
    m2 := NaturalHomomorphismByNormalSubgroupNC(Sldq,C);
    pS := Image(m2); #PSL(d,q),PSp(d,q),or POmega+(2d,q)

    # Calculating StopS
    StopS := IsomorphismSimpleGroups(S,pS); # Isomorphism between S and pS
    if StopS = fail then StopS := IsomorphismGroups(S,pS);fi;

    # Calculating AutpSgens
    AutpSgens := List(AutSgens,x -> InducedAutomorphism(StopS,x));

    # Function to get alpha for any lambda
    m2lambm2inv := function(U,lambda,m2,C,p)
        local UZ,VZ,VZrep,V;
        UZ := Image(m2,U);
        VZ := Image(lambda,UZ);
        VZrep := PreImagesRepresentative(m2,VZ);
        if Size(C)=1 then return VZrep;fi;
        VZ := RightCoset(C,VZrep);
        for V in VZ do
            if Order(V) = p then return V; fi;
        od;
        Error("No elemnt found");
    end;

    # Calculating AutSl
    AutSlgens := [];
    for lambda in AutpSgens do
        alphaSlgen := List(Slgen,U -> m2lambm2inv(U,lambda,m2,C,p));
        alpha := GroupHomomorphismByImagesNC(Sldq,Sldq,Slgen,alphaSlgen);
        Add(AutSlgens,alpha);
    od;

    # Checking if AutSlgens contains an automorphism that is not a
    # conjugator automorphism times a field automorphism.
    scalars := Filtered(GF(q), x -> x^d = One(GF(q)));
    scalars := SmallGeneratingSet(AsGroup(scalars));
    scalarMatrices := List(scalars,x -> x * IdentityMatrix(GF(q),d));
    f := FrobeniusAutomorphism(GF(q));
    flisinit := List([1..Order(f)],x->f^x);
    finvlisinit := List(flisinit,f -> f^-1);
    imagesOfscalars := List(flisinit,f -> List(scalars,f));
    MakeFieldAutomorphismOfSldq := function(f)
      local row,M,applyfonrow,applyf,fgens,x;
      applyfonrow := row -> List(row,x -> Image(f,x));
      applyf := M -> List(M,applyfonrow);
      fgens := List(Slgen,applyf);
      return GroupHomomorphismByImagesNC(Sldq,Sldq,Slgen,fgens);
    end;
    finvmatrixlist := List(finvlisinit,MakeFieldAutomorphismOfSldq);
    if type="PSL" or type="PSp" then
        for alpha in AutSlgens do
          conj := false;
          FirstOfimagesOfscalarMatrices := List(scalarMatrices,M -> alpha(M)[1][1]);
          findlis := Filtered([1..Length(flisinit)], i -> FirstOfimagesOfscalarMatrices = imagesOfscalars[i]);
          finvlis := List(findlis,i -> finvmatrixlist[i]);
          for finv in finvlis do
            alphafinv := alpha * finv;
            if IsConjugatorAutomorphism(alphafinv) then
              conj := true;
              break;
            fi;
          od;
          if not conj then return false; fi;
        od;
    elif type="POmega+" then
        O := GO(1,d,q);
        C := Center(GO(1,d,q));
        phi := NaturalHomomorphismByNormalSubgroupNC(O,C);
        PO := Image(phi);
        POgens := GeneratorsOfGroup(PO);
        POgens := List(POgens,x -> PreImagesRepresentative(phi,x));
        POgens := List(POgens,x-> ConjugatorAutomorphismNC(Sldq,x));
        PO := GroupWithGenerators(POgens);
        if d > 8 then
            for alpha in AutSlgens do
                if not alpha in PO then return false; fi;
            od;
        elif d = 8 then

        fi;
    fi;
    return true;
end;


NonClassicalLieGroupThings := function(d,q,type)
K :=GF(q);
L := SimpleLieAlgebra(type,d,K);
R := RootSystem(L);
PhiP :=PositiveRoots(R);
PhiN := NegativeRoots(R);
VecPhiP :=PositiveRootVectors(R);
VecPhiN :=NegativeRootVectors(R);
VecPi :=SimpleSystem(R);
LK := TensorProductOfAlgebraModules(K,L);
o := One(K);
PhiP := List(PhiP,x->TensorProduct(o,k));
end;

##################################################################################
##
#F  MinimalFaithfulPermutationDegreeOfAlmostSimpleGroupWithSimpleSubgroup(<A,S>)
##
##  Returns mu(A) where S < A < Aut(S) for simple group S
##  if it can compute it easily, else returns -1.
##  For more details, take a look at this paper :
##  https://www.sciencedirect.com/science/article/pii/S0747717118300993
MinimalFaithfulPermutationDegreeOfAlmostSimpleGroupWithSimpleSubgroup := function(A,S)
    local info,d,q,m,series,name,parameter,mu,Aut,b,sizeS,sizeA;

    # Just for convenience
    Aut := AutomorphismGroup;
    mu := MinimalFaithfulPermutationDegreeOfSimpleGroup;

    if IsInt(A) then sizeA := A;
    else sizeA := Size(A); fi;
    sizeS := Size(S);
    info := IsomorphismTypeInfoFiniteSimpleGroup(S);
    series := info.series;
    name := info.shortname;

    if series = "Spor" then 
        if name = "M12" and
            sizeA = Size(Aut(S))  #A = Aut(S)
            then return 2 * mu(S);
        elif name = "ON" and
            sizeA = Size(Aut(S))  #A = Aut(S)
            then return 2 * mu(S);
        fi;
    else # Set q and m
        parameter := info.parameter;
        if IsList(parameter) then
            q := parameter[2];
            m := parameter[1];
        else q := parameter; fi;
    fi;

    if series = "A" and q = 6 then
        #if not A <~ SymmetricGroup(6) then return 10; fi;
        if 720 mod sizeA <> 0 then return 10; 
        else return -1; fi;

    elif series = "L" then
        d := m;
        if d = 2 and q = 7 then
            # if A ~ PGL(2,7) then return 8; fi;
            # PGL(2,7) is isomorphic to Aut(PSL(2,7)) which is Aut(S) which has a subgroup isomorphic to A
            # Thus A <~ PGL(2,7); So, the only requirement for isomorphism is that Size(A) = Size(PGL(2,7)) = 336
            if sizeA = 336 then return 8; fi;
        elif d > 2 and not (q = 2 and (d = 3 or d = 4)) then
            #if not A <~ GammaL(d,q) then return 2 * mu(S);
            if Size(GammaL(d,q)) mod sizeA <> 0 then return 2 * mu(S);
            else return -1; fi;
        fi;

    elif series = "2A" # PSU(3,5)
        and m+1 = 3
        and q = 5 then
        #if not A <~ PSigmaU(3,5) then return 126; fi;
        return -1;

    elif series = "D" then #P\Omega^+
        d := 2*m;
        if d = 8 and q = 2
            and sizeA/sizeS mod 3 = 0
            then return 3 * mu(S);
        elif d = 8 and q = 3
            and sizeA/sizeS mod 3 = 0
            and sizeA/sizeS mod 12 <> 0
            then return 3 * mu(S);
        elif d = 8 and q = 3
            and sizeA/sizeS mod 12 = 0
            then return 3360;
        elif d = 8
            and sizeA/sizeS mod 3 = 0
            then return 3 * mu(S);
        elif q = 3 and d > 7 and
            sizeA/sizeS mod 3 <> 0 then
            if Size(PGO(1,d,3)) mod sizeA <>0 # (*)
            then return (3^(m-1) + 1) * (3^m -1) / 2 ;
            else return -1; fi;
        fi;

    elif series = "G" then
        #if q = 3 and Iso(A,Aut(S)) then return 2 * mu(S);fi;
        if q = 3 and sizeA = Size(Aut(S)) then return 2 * mu(S); fi;
        #This is because A is isomorphic to a subgroup in Aut(S) that contains Inn(S) \cong S at any point.
        b := LogInt(q,3);
        if 3^b = q and b >1 then
            #if not A <~ GammaG(q) then return 2 * mu(S); fi;
            if Size(GF(q))*Size(GaloisGroup(GF(q))) mod sizeA <> 0 then return 2 * mu(S);
            else return -1; fi;
        fi;
    elif series = "C" #PSp
        and m = 2
        and 2^Log2Int(q) = q
        and q > 3 then
        #if not A <~ PGammaSp(d,q) then return 2*mu(S); fi;
        return -1;
    elif series = "F" and 2^Log2Int(q) = q then
        #if not A <~ GammaF(d,q) then return 2*mu(S); fi;
        if sizeS * Size(GaloisGroup(GF(q))) mod sizeA <> 0 then return 2 * mu(S);
        else return -1; fi;
    elif series = "E" and m = 6 then
        #if not A <~ GammaE(d,q) then return 2*mu(S); fi;
        if sizeS * Size(GaloisGroup(GF(q))) mod sizeA <> 0 then return 2 * mu(S);
        else return -1; fi;
    fi;
    return mu(S);
end;

##################################################################################
##
#F  MinimalFaithfulPermutationDegreeOfSemiSimpleGroup(<G>)
##
##  Returns the minimal faithful permutation degree for a semi-simple group G,
##  which is a group that doesn't contain any abelian normal subgroups.
##  The function is based on this research paper :
##  https://dl.acm.org/doi/10.1145/3618260.3649641
MinimalFaithfulPermutationDegreeOfSemiSimpleGroup := function(G)
    local
        mu,     # Objective value to compute
        Slis,   # Decomposition of N as direct product of simple groups
        N,      # Any minimal normal subgroup of G
        NGS,    # Normaliser of S in G
        CGS,    # Centraliser of S in G
        phi,    # Homomorphism from NGS to NGS/CGS
        A,      # NGS/CGS =~ {^n for n in N_G(S)}
        muA,    # Minimal Permutation Degree of A
        sizeA,  # Size(A)
        S,      # Slis[1] =~ Slis[i]
        info,   # Information about isomorphism type of S
        MNS,m,q;
    mu := 0;
    MNS := MinimalNormalSubgroups(G);
    if Length(MNS) = 1 and MNS[1] = G then return MinimalFaithfulPermutationDegreeOfSimpleGroup(G); fi;
    for N in MNS do
        Slis := DirectFactorsOfGroup(N);
        S := Slis[1];
        NGS := Normalizer(G,S);
        CGS := Centralizer(G,S);
        sizeA := Size(NGS)/Size(CGS);
        muA := MinimalFaithfulPermutationDegreeOfAlmostSimpleGroupWithSimpleSubgroup(sizeA,S);
        if muA = -1 then
            info := IsomorphismTypeInfoFiniteSimpleGroup(S);
            if not IsList(info.parameter) then q := info.parameter;
            else m := info.parameter[1]; q := info.parameter[2]; fi;
            if info.series = "L" #PSL
                and m > 2
                and not (q = 2 and (m in [3,4])) then
                #Print("\t PSL type simple subgroup\n");
                if CanLiftPermutationDegreeToNormaliserOfSimpleMatrixGroup(S,G,m,q,"PSL") then
                    muA := MinimalFaithfulPermutationDegreeOfSimpleGroupWithIsomorphismType(info);
                else
                    muA := 2*MinimalFaithfulPermutationDegreeOfSimpleGroupWithIsomorphismType(info);
                fi;
            elif info.series ="C" #PSp
                and m = 2
                and 2^Log2Int(info.parameter[2]) = info.parameter[2] then
                if CanLiftPermutationDegreeToNormaliserOfSimpleMatrixGroup(S,G,m,q,"PSp") then
                    muA := MinimalFaithfulPermutationDegreeOfSimpleGroupWithIsomorphismType(info);
                else
                    muA := 2*MinimalFaithfulPermutationDegreeOfSimpleGroupWithIsomorphismType(info);
                fi;
            elif info.series ="D" #POmega+
                and q = 3
                and m >= 4 then
                if CanLiftPermutationDegreeToNormaliserOfSimpleMatrixGroup(S,G,m,q,"POmega+") then
                    muA := MinimalFaithfulPermutationDegreeOfSimpleGroupWithIsomorphismType(info);
                else
                    muA := (3^(m-1) + 1) * (3^m -1) / 2 ;
                fi;
            else
                #Print("\t Leaving to Lattice Algo\n");
                phi := NaturalHomomorphismByNormalSubgroupNC(NGS,CGS);
                A := Image(phi);
                muA := MinimalFaithfulPermutationDegree(A);
            fi;
        fi;
        mu := mu + (Length(Slis) * muA);
    od;
    return mu;
end;

# tests

CheckInList := function(L,Mus,compare)
    local i,stop,interval,limit,size,GLis,t0,t,mu,mu2,G,G1,G2,sizenindex,g;
    stop := false;
    for i in [1..Length(L)] do
        G := L[i];
        if IsList(Mus) then mu := Mus[i];
        else mu := -1; fi;
        if stop then break; fi;

        if Length(String(G)) < 7 then Print(G," : ");
        else Print("Group ",i," : "); fi;
        if IsSimple(G) then Print(" (Simple)");fi;
        Print("\n");

        FlushCaches();
        t0 := Runtime();
        mu2 := MinimalFaithfulPermutationDegreeOfSemiSimpleGroup(G);
        t := Runtime();
        Print("\t dhara algo time ",t-t0,"\n");

        if compare or mu = -1 then 
            g := GeneratorsOfGroup(G);
            if not IsList(g) then g := [g];fi;
            G2 := GroupByGenerators(g);
            FlushCaches();
            t0 := Runtime();
            mu := DoMinimalFaithfulPermutationDegree(G2,false);
            t := Runtime();
            Print("\t lattice algo time ",t-t0,"\n");
        fi;

        Print("\t Result : ");
        if mu2 = mu then
            Print("pass\n\n");
        else
            Print(" F.A.I.L \n");
            Print("\t",G,"\n");
            stop := true;
            break;
        fi;
    od;
end;

#TestCases := Filtered(List(SimpleGroupsIterator(2,50000)),G->IsPSL(G) or (Size(G) < 1000 and not IsAbelian(G)));
TestCases := List(SimpleGroupsIterator(2,300));
Degs := List(TestCases,MinimalFaithfulPermutationDegree);
l := Length(TestCases);
for i in [1..l] do
    for j in [1..i] do
        G := DirectProduct(TestCases[i],TestCases[j]);
        if Gcd(Size(TestCases[i]),Size(TestCases[j])) = 1 or true then
            Add(TestCases,G);
            Add(Degs,Degs[i] + Degs[j]);
        fi;
    od;
od;
SemiSimpleGroupsTillSize1000 := [
    SmallGroup(60,5),
    SmallGroup(120,34),
    SmallGroup(168,42),
    SmallGroup(336,208),
    SmallGroup(360,118),
    SmallGroup(504,156),
    SmallGroup(660,13),
    SmallGroup(720,763),
    SmallGroup(720,764),
    SmallGroup(720,765),
    #sizes 512 and 768 skipped
];
BadApples := [
    [PSL(2,7), "PSL(2,7)", 168, 7 ],
    [AlternatingGroup(6), "A6", 360, 6 ],
    [PSL(3,3), "PSL(3,3)", 5616, 13 ],
    [SimpleGroup("M12"), "M12", 95040, 12 ],
    [PSU(3,5), "PSU(3,5)", 126000, 50 ],
    [PSp(4,4), "PSp(4,4)", 979200, 85 ],
    [SimpleGroup("G",2,3), "G(2,3)", 4245696, 351 ],
    [PSL(4,3), "PSL(4,3)", 6065280, 40 ],
    [PSL(5,2), "PSL(5,2)", 9999360, 31 ],
    [POmega(1,8,2), "POmega+(8,2)", 174182400, 120 ],
    [PSp(4,8), "PSp(4,8)", 1056706560, 585 ],
    #[SimpleGroup("ON"), "ON", 460815505920, 122760 ], # skipping
    [POmega(1,8,3), "POmega+(8,3)", 4952179814400, 1080 ],
    [SimpleGroup("G",2,9), "G(2,9)", 22594320403200, 66430 ],
    [SimpleGroup("F",4,2), "F(4,2)", 3311126603366400, 69615 ],
    [POmega(1,8,4), "POmega+(8,4)", 67010895544320000, 5525 ],
    [POmega(1,8,5), "POmega+(8,5)", 8911539000000000000, 19656 ],
    [POmega(1,10,3), "POmega+(10,3)", 1289512799941305139200, 9801 ]
  ];
CheckBadApples := function()
    local B,i,j,mu,G,mu2;
    for i in [1..Length(BadApples)] do
        B := BadApples[i][1];
        Print("\n",BadApples[i][2]," : \n\n");
        for j in [1..Length(TestCases)] do
            Print(String(j)," : ");
            mu := Degs[j] + BadApples[i][4];
            G := DirectProduct(B,TestCases[j]);
            mu2 := MinimalFaithfulPermutationDegreeOfSemiSimpleGroup(G);
            if not mu=mu2 then Print(" FAIL \n"); break;
            else Print("pass \n"); fi;
        od;
    od;
end;

Print("\n\nSemi Simple Groups below size 1000 : \n\n");
CheckInList(SemiSimpleGroupsTillSize1000,-1,true);

Print("\n\nProducts of Simple Groups : \n\n");
CheckInList(TestCases,Degs,false);

Print("\n\n Corner cases : \n\n");
CheckBadApples();