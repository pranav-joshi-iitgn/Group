from numpy import *

class Relation:
    def __init__(self,T:ndarray,Elements=None,calculate=True,level=1) -> None:
        self.level = level
        if Elements is None:
            Elements = [f"x_{i}" for i in range(T.shape[0])]
        if T.dtype != int32:
            try:
                T = T.astype(int32)
            except:
                raise TypeError("Not correct datatype")
        if (T.shape[0] != T.shape[1]):
            raise TypeError("Matrix must be square")
        if T.shape[0] != len(Elements):
            raise TypeError("Provide names for all elements please")
        
        self.T = T
        self.Elements = Elements
        if calculate:
            self.calculate_properties()
            self.calculated = True
        else:
            self.calculated = False
    def inverses(self,ElId):
        if not self.calculated:
            self.calculate_properties()
        inv = [None]*len(ElId)
        for i in range(len(ElId)):
            x = ElId[i]
            inv[i] = self.inv_ind[x]
        return inv
    def __getitem__(self,i):
        return self.Elements[i]
    def __repr__(self) -> str:
        if sum([len(str(name)) for name in self.Elements]) > 100:
            return "Relation("+",".join([str(name) for name in self.Elements])+")"
        s = " "*max([len(str(name)) for name in self.Elements])+"\t   "
        for name in self.Elements:
            s += str(name) + '\t'
        s += "\n"+" "*max([len(str(name)) for name in self.Elements])+"\t   "
        for name in self.Elements:
            s += '-'*len(str(name)) + '\t'
        i = 0
        for name in self.Elements:
            s += "\n"+str(name)+"\t|  "
            for j in range(self.T.shape[0]):
                s += str(self.Elements[self.T[i][j]]) + '\t'
            i += 1
        return s
    def Element(self,i):
        return Element(self,i)
    def calculate_properties(self):
        self.identity_ind = get_identity(self)
        if self.identity_ind is None:
            print("Relation doesn't have an identity")
            return
        self.inv_ind = get_inverses(self,self.identity_ind)
    def generate(self,generators,check=True,ret_vis=False,visited=None):
        T=self.T
        n = T.shape[0]
        if visited is None:
            visited = [False]*n
        visited[self.identity_ind] = True
        while True:
            changes = 0
            for x in range(n):
                if not visited[x]:
                    continue
                for g in generators:
                    try:
                        gx = T[g][x]
                    except:
                        raise TypeError(f"g={g},x={x}")
                    if not visited[gx]:
                        visited[gx] = True 
                        changes += 1
            if not changes:
                break
        if ret_vis:
            return visited
        L = [i for i in range(n) if visited[i]]
        try:
            G = Group(self,L,check)
        except:
            raise TypeError("Can't make a group from"+str(L)+"generated from"+str(generators))
        return G
    def __mul__(R1,R2):
        T1 = R1.T
        T2 = R2.T
        E1 = R1.Elements
        E2 = R2.Elements
        n1 = len(E1)
        n2 = len(E2)
        E = [None]*(n1*n2)
        for i in range(n1):
            for j in range(n2):
                x = E1[i]
                y = E2[j]
                if isinstance(x,tuple):
                    if isinstance(y,tuple):
                        z = x + y
                    else:
                        z = x + (y,)
                else:
                    if isinstance(y,tuple):
                        z = (x,) + y
                    else:
                        z = (x,y)
                E[n2*i+j] = z
        n = len(E)
        T = [[None]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                ai,bi = i//n2,i%n2
                aj,bj = j//n2,j%n2
                ak = T1[ai][aj]
                bk = T2[bi][bj]
                k = ak*n2 + bk
                T[i][j] = k
        T = array(T)
        R = Relation(T,E)
        return R

class Element:
    def __init__(self,R:Relation,i:int) -> None:
        self.R = R
        self.i = i
        self.element = R.Elements[i]
    def __repr__(self) -> str:
        return str(self.element)
    def __mul__(self,other):
        if other.R is not self.R:
            raise TypeError("Not having same parent group.")
        ind = self.R.T[self.i][other.i]
        return Element(self.R,ind)

def get_identity(R:Relation,ElInd=None)->bool:
    T = R.T
    n = T.shape[0]
    if ElInd is None:
        ElInd = range(n)
    for x in ElInd:
        is_x_e = True
        for y in ElInd:
            xy = T[x][y]
            if xy != y:
                is_x_e = False
                break
        if is_x_e:
            return x
    return None    
def is_associative(R:Relation,ElInd=None)->bool:
    T = R.T
    n = T.shape[0]
    if ElInd is None:
        ElInd = range(n)
    for x in ElInd:
        for y in ElInd:
            for z in ElInd:
                yz = T[y][z]
                xy = T[x][y]
                x_yz = T[x][yz]
                xy_z = T[xy][z]
                if x_yz != xy_z:
                    return False
    return True
def get_inverses(R:Relation,identity_ind:int,ElInd=None)->list:
    T= R.T
    n = T.shape[0]
    if ElInd is None:
        ElInd = range(n)
    n=len(ElInd)
    inverses = [None]*n
    for i in range(n):
        x = ElInd[i]
        for y in ElInd:
            if T[x][y] == identity_ind:
                inverses[i] = y
                break
    return inverses
def contains_duplicates(L)->bool:
    n = len(L)
    for i in range(n):
        for j in range(i):
            if L[i] == L[j]:
                return True
    return False
def is_commutative(R:Relation,ElInd=None)->bool:
    T = R.T
    n = T.shape[0]
    if ElInd is None:
        ElInd = range(n)
    Ts = T[ElInd][:,ElInd]
    return array_equal(Ts.T,Ts)

class Group:
    def __init__(self,R:Relation,ElInd=None,check=True) -> None:
        self.level = R.level
        T = R.T
        n = T.shape[0]
        if ElInd is None:
            ElInd = list(range(n))
        self.ElInd = ElInd
        self.R = R
        if not R.calculated:
            R.calculate_properties()
        id = R.identity_ind
        if id is None:
            raise TypeError("The relation doesn't have an identity element")
        inv = R.inverses(ElInd)
        assert len(inv) == len(ElInd)
        if None in inv:
                raise TypeError("An element doesn't have it's inverse.")
        if check:
            if not is_associative(R,ElInd):
                raise TypeError("The relation provided isn't associative")
            #inv = get_inverses(R,id,ElInd)
            if contains_duplicates(inv):
                raise TypeError("Two elements have the same inverses")
        self.identity_ind = id
        self.inv_ind = inv
    def perform_checks(self):
        ElInd = self.ElInd
        R = self.R
        id = get_identity(R,ElInd)
        if id is None:
            print("The relation doesn't have an identity element")
            return False
        if not is_associative(R,ElInd):
            print("The relation provided isn't associative")
            return False
        inv = get_inverses(R,id,ElInd)
        if None in inv:
            print("An element doesn't have it's inverse.")
            return False
        if contains_duplicates(inv):
            print("Two elements have the same inverses")
            return False
        self.identity_ind = id
        self.inv_ind = inv
        return True
    def __getitem__(self,i):
        return self.R.Elements[self.ElInd[i]]
    def __repr__(self) -> str:
        strEl = [str(self.R.Elements[i]) for i in self.ElInd]
        l = sum([len(s) for s in strEl])
        if l>100:
            return "Group(\n"+",\n".join(strEl)+")"
        return "Group("+",".join(strEl)+")"
    def Element(self,i):
        return Element(self.R,i)
    def has_subgroup(self,H,check=False,debug=True):
        if not isinstance(H,Group):
            if debug:
                print("Not a group")
            return False
        if self.R is not H.R:
            if debug:
                print("Working on different relations :",self.R.Elements,H.R.Elements,sep='\n')
            return False
        for h in H.ElInd:
            if h not in self.ElInd:
                if debug:
                    print("An element in H is not in G")
                return False
        if check and not H.perform_checks():
            return False
        return True
    def has_normal_subgroup(self,N,check=False):
        assert isinstance(N,Group)
        if not self.has_subgroup(N,check):
            return False
        R = self.R
        for g in self.ElInd:
            if g in N.ElInd:
                continue
            giNg = N.conjugate(g)
            if not giNg == N:
                return False
        return True
    def conjugate(self,g_ind):
        ElInd = list(self.ElInd)
        n = len(ElInd)
        newElInd = [None]*n
        R = self.R
        T = self.R.T
        g = g_ind
        g_inv = R.inv_ind[g]
        for i in range(n):
            x = ElInd[i]
            xg = T[x][g]
            g_inc_xg = T[g_inv][xg]
            newElInd[i] = g_inc_xg
        assert None not in newElInd
        giGg = Group(self.R,newElInd,False)
        return giGg
    def __eq__(self,other):
        if not isinstance(other,Group):
            return False
        if self.R is not other.R: 
            return False
        n1 = len(self.ElInd)
        n2 = len(other.ElInd)
        if n1 != n2:
            return False
        ell1 = sorted(self.ElInd)
        ell2 = sorted(other.ElInd)
        for i in range(n1):
            if ell1[i]!=ell2[i]:
                return False
        return True
    def Cosets(self,H,return_division=False)->list:
        assert self.has_subgroup(H)
        R = self.R
        T = R.T
        n = T.shape[0]
        visited = [0]*n
        G = self.ElInd
        cosets = []
        i = 1
        for g in G:
            if visited[g]:
                continue
            gH = Coset(H,g) #this will by default be in standard form
            for gh in gH.expand():
                visited[gh] = i
            cosets.append(gH)
            i += 1
        if return_division:
            return cosets,visited
        return cosets
    def __truediv__(self,N):
        G = self
        assert self.has_subgroup(N)
        assert self.has_normal_subgroup(N)
        cosets,division = self.Cosets(N,True)
        T = self.R.T
        M = []
        n= len(cosets)
        for i in range(n):
            g1N = cosets[i]
            g1 = g1N.g_ind
            row = [None]*n
            for j in range(n):
                g2N = cosets[j]
                g2 = g2N.g_ind
                g1g2 = T[g1][g2]
                row[j] = division[g1g2]-1
            assert None not in row, str((row,g1))
            M.append(row)
        M = array(M)
        Rnew = Relation(M,cosets,level=self.level+1)
        GbyN = Group(Rnew)
        return GbyN
    def normal_closure(self,xL):
        if not iterable(xL):
            xL = [xL]
        R = self.R
        T = R.T
        n = T.shape[0]
        visited = [False]*n
        for x in xL:
            found = False
            for i in range(len(self.ElInd)):
                g = self.ElInd[i]
                g_inv = self.inv_ind[i]
                if g==x:
                    found = True
                xg = T[x][g]
                gixg = T[g_inv][xg]
                visited[gixg] = True
            assert found, str((x,g))
        L = [i for i in range(n) if visited[i]]
        xtoG = R.generate(L,False)
        return xtoG
    def MinimumNormalSubGroup(self):
        R = self.R
        T = R.T
        #n = T.shape[0]
        min_l = len(self.ElInd)
        min_N = self
        for x in self.ElInd:
            if x == self.identity_ind:
                continue
            N = self.normal_closure(x)
            l = len(N.ElInd)
            if l < min_l:
                min_l = l
                min_N = N
        return min_N
    def MinimalNormalSubGroup(self):
        R = self.R
        T = R.T
        n = T.shape[0]
        visited = [False]*n
        visited[self.identity_ind] = True
        ToDo = self.ElInd
        N = self
        while True:
            for x in ToDo:
                if visited[x]:
                    continue
                N = self.normal_closure([x])
                ToDoNew = N.ElInd
                visited[x] = True
                sam = (len(ToDo) == len(ToDoNew))
                #sam = True
                #for i in range(n):
                #    same[i] = False
                #for el in ToDoNew:
                #    same[el] = True
                #for el in ToDo:
                #    if not same[el]:
                #        sam = False
                #        break
                if not sam:
                    ToDo = ToDoNew
                    break
            return N      
    def GeneratingSet(self)->list:
        if len(self.ElInd)==1:
            return self.ElInd
        R = self.R
        T=R.T
        n = T.shape[0]
        visited = [False]*n
        visited[self.identity_ind] = True
        generators = []
        for g in self.ElInd:
            if visited[g]:
                continue
            generators.append(g)
            visited = R.generate(generators,False,True,visited)
            #gps = []
            #x = g
            #while not visited[x]:
            #    gps.append(x)
            #    x = T[x][g]
            #for x in self.ElInd:
            #    if not visited[x]:
            #        continue
            #    for y in self.ElInd:
            #        if not visited[y]:
            #            continue
            #        for gp in gps:
            #            gx = T[gp][x]
            #            ygx = T[y][gx]
            #            if not visited[ygx]:
            #                visited[ygx] = True 
        return generators    
    def MinimalGeneratingSet(self,gen=None):
        if gen is None:
            gen = self.GeneratingSet()
        i = 0
        while i < len(gen):
            genn = gen[:i]+gen[i+1:]
            if self.has_generating_set(genn):
                gen = genn
            else:
                i += 1
        return gen
    def order(self,i:int):
        G = self
        I = i
        for o in range(1,len(G)+1):
            if G.identity_ind == I:
                return o
            I = G.R.T[I][i]
        print("Not good")
    def cyclic_generator(G):
        n = G.R.T.shape[0]
        visited = [False]*n
        visited[G.identity_ind] = True
        for i in G.ElInd:
            if visited[i]:
                continue
            I = i
            broken = False
            for o in range(1,len(G)):
                visited[I] = True
                if G.identity_ind == I:
                    broken = True
                    break
                I = G.R.T[I][i]
            if not broken:
                return i
        return None
    def is_cyclic(self):
        return (self.cyclic_generator is not None)
    def is_simple(self):
        if self.MinimumNormalSubGroup() == self:
            return True
        return False  
    def is_abelian(self):
        return is_commutative(self.R,self.ElInd)
    def has_generating_set(self,generating_set):
        G = self.R.generate(generating_set)
        return self==G 
    def __len__(self):
        return len(self.ElInd)
    def explode(self,gL,t)->list:
        if t==0:
            return [gL]
        old = self.explode(gL[1:],t-1)
        gN = Coset(self,gL[0])
        gNL = gN.expand()
        new = old.copy()
        for c in old:
            new.extend([([gn]+c) for gn in gNL if ((gn not in c) and (gn != self.identity_ind))])
        return new
    def __pow__(self,n:int):
        """
            Finding G^n via fast exponentiation
        """
        L = [None]*(int(log2(n))+1)
        B = [None]*(int(log2(n))+1)
        L[0]=self
        for i in range(len(L)-1):
            L[i+1] = L[i].Cross(L[i])
        G = None
        for i in range(len(B)):
            B[i] = n%2
            n = n//2
            if B[i]:
                if G is None:
                    G = L[i]
                else:
                    G = G.Cross(L[i])
        return G    
    def Cross(G1,G2):
        R1 = G1.R
        R2 = G2.R
        R = R1*R2
        E1 = G1.ElInd
        E2 = G2.ElInd
        n1 = R1.T.shape[0]
        n2 = R2.T.shape[0]
        l1 = len(E1)
        l2 = len(E2)
        E = [None]*(l1*l2)
        for i in range(l1):
            for j in range(l2):
                x = E1[i]
                y = E2[j]
                z = x*n2 + y
                E[i*l2+j] = z
        G = Group(R,E)
        return G
    def MinimumGeneratingSet(self,debug=False):
        G = self
        if len(G) == 1:
            if debug:
                print("G has only identity")
            return G.ElInd
        if debug:
            print("Finding minimum generating set for :\n",self)
        g = G.cyclic_generator()
        if g is not None:
            if debug:
                print("G is cyclic")
            return [g]
        T = G.R.T
        N = G.MinimalNormalSubGroup()
        if debug:
            print("N :\n",N)
        if G==N:
            for x in G.ElInd:
                if x==G.identity_ind:
                    continue
                if G.has_generating_set([x]):
                    return [x]
            for x in G.ElInd:
                if x==G.identity_ind:
                    continue
                for y in G.ElInd:
                    if y==G.identity_ind:
                        continue
                    if G.has_generating_set([x,y]):
                        return [x,y]
            print("This shouldn't be happening")
        n = N.MinimalGeneratingSet()
        if debug:
            print("n :\n",n)
        m = len(n)
        GbyN = G/N
        if debug:
            print("G/N :\n",GbyN)
        mingenGbyN = GbyN.MinimumGeneratingSet(debug)
        g = [(GbyN[i]).g_ind for i in mingenGbyN]
        if debug:
            print("g :\n",g)
        l = len(g)
        if N.is_abelian():
            if debug:
                print("N is abelian.")
            #proceed = not G.has_generating_set(g)
            for i in range(l):
                for j in range(m):
                    modified_g = g[:i]+[T[g[i]][n[j]]]+g[i+1:]
                    for el in modified_g:
                        assert (isinstance(el,int) or isinstance(el,int32)),modified_g
                    if G.has_generating_set(modified_g):
                        return modified_g
            for i in N.ElInd:
                if N.identity_ind != i:
                    return g + [i]
            else:
                print("idk wat to do now.")
        else:
            t = 13/5 + log(len(G))/log(len(N))
            if debug:
                print("N is not abelian")
                print("t = ",t,'l = ',l)
            if t<= l:
                for candidate_gen in N.explode(g,t):
                    if G.has_generating_set(candidate_gen):
                        return candidate_gen
            else:
                for candidate_gen0 in N.explode(g,l):
                    for nl in N.ElInd:
                        candidate_gen = candidate_gen0 + [nl]
                        if G.has_generating_set(candidate_gen):
                            return candidate_gen

class Coset:
    """
    Class for left cosets
    """
    def __init__(self,H:Group,g_ind:int,G=None,check=False) -> None:
        self.H = H
        if check :
            assert G.has_subgroup(H) , "H is not a subgroup of G"
        self.g_ind = g_ind
    def expand(self)->list:
        g = self.g_ind
        H = self.H
        R = H.R
        T = R.T
        H_el_ind = H.ElInd
        n = len(H_el_ind)
        elements = [None]*n
        for i in range(n):
            h = H_el_ind[i]
            elements[i] = T[g][h]
        return elements
    def __repr__(self) -> str:
        H = self.H
        g = self.H.R[self.g_ind]
        sH = str(H)
        if len(sH) > 20:
            sH = f"Group[level={H.level}]"
        return "("+str(g)+" * "+sH+")"
    def __mul__(self,other):
        assert isinstance(other,Coset),"Cosets can only be multiplied with cosets"
        assert other.H is self.H,"The cosets must have the same subgroup object (H) generating the"
        g1 = self.g_ind
        g2 = other.g_ind
        H = self.H
        R = H.R
        T = R.T
        g3 = T[g1][g2]
        return Coset(H,g3)
    def standardise(self):
        self.g_ind = min(self.expand())

def AdditiveGroupOnIntegersModulo(n:int)->Group:
    T = array([
        [(i+j)%n for j in range(n)]
        for i in range(n)
    ])
    names = [i for i in range(n)]
    R = Relation(T,names)
    return Group(R)
def DihegralGroup(n:int)->Group:
    names = ['e']+[f'r^{i}' for i in range(1,n)] + ['s']+[f'sr^{i}' for i in range(1,n)]
    T = [[None]*(2*n) for _ in range(2*n)]
    for i in range(2*n):
        for j in range(2*n):
            if i<n and j<n:
                T[i][j] = (i+j)%n
            elif i>=n and j<n:
                T[i][j] = n+(i+j)%n
            elif i<n and j>=n:
                T[i][j] = n+(n-i+j)%n
            elif i>=n and j>=n:
                T[i][j] = (n-i+j)%n
    for row in T:
        assert None not in row, str(T)
    T = array(T,dtype=int32)
    R = Relation(T,names,True)
    return Group(R)
def generatePermuations(n:int)->list:
    """
    Generates permutations of numbers 0,1,...n-1 in increasing lexicographical order
    """
    if n==0:
        return []
    if n==1:
        return [[0]]
    Lo = generatePermuations(n-1)
    lo = len(Lo)
    l = lo*n
    L = [None]*l
    for i in range(n):
        for j in range(lo):
            po = Lo[j]
            p = [None]*n
            p[0] = i
            for r in range(1,n):
                if po[r-1]>=i:
                    p[r]=1+po[r-1]
                else:
                    p[r]=po[r-1]
            assert None not in p
            L[i*lo+j] = p

    assert None not in L
    return L
def perutationIndex(p:int,facts:list):
    """
    facts is a list containing factorials of numbers 0,1,...n-1 where n = len(p)
    """
    p = p.copy()
    n = len(p)
    assert len(facts)==n
    I = 0
    for i in range(n):
        x = p[i]
        I += x*facts[n-1-i]
        for j in range(i+1,n):
            if p[j] > p[i]:
                p[j] -= 1
            elif p[j] == p[i]:
                raise ValueError(str((p[j],p[i],j,i)))
    return I
def multiplyPermutations(x,y):
    """
    returns xy where x and y are permutations of same size
    """
    assert len(x) == len(y)
    z = [y[i] for i in x]
    return z
def PermutationGroup(n:int)->Group:
    """
    generates S(n)
    """
    facts = [None]*n
    facts[0] = 1
    for i in range(1,n):
        facts[i] = i*facts[i-1]
    P = generatePermuations(n)
    P = [tuple(x) for x in P]
    l = len(P)
    #D = {P[i]:i for i in range(l)}
    T = [[None]*l for _ in range(l)]
    for i in range(l):
        for j in range(l):
            x = P[i]
            y = P[j]
            z = multiplyPermutations(x,y)
            zind = perutationIndex(z,facts)
            #zind = D[z]
            T[i][j] = zind
    T = array(T)
    R = Relation(T,["["+"".join([str(x) for x in p])+']' for p in P])
    G = Group(R)
    return G
def AlternatingGroup(n:int)->Group:
    P = PermutationGroup(n)
    A = P.MinimalNormalSubGroup()
    assert len(P) == 2*len(A)
    return A
def RelationFromElements(L:list,group_op=lambda x,y:x*y)->Relation:
    """
    L:list of group elements
    group_op: Group operation over these elements. 
    """
    n = len(L)
    T = [[None]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            x = L[i]
            y = L[j]
            xy = group_op(x,y)
            xyi = L.index(xy)
            T[i][j] = xyi
    T = array(T)
    R = Relation(T,L)
    return R
def TableGroupFromElements(L:list,group_op=lambda x,y:x*y)->Group:
    R = RelationFromElements(L,group_op)
    G = Group(R)
    return G

if __name__=='__main__':
    from matplotlib.pyplot import *
    from time import time
    G0 = AdditiveGroupOnIntegersModulo(2)
    N = 6
    iterations = 3
    G = G0
    L = []
    for n in range(1,N+1):
        t0 = time()
        for i in range(iterations):
            g = G.MinimumGeneratingSet()
        t = time()
        dt = (t-t0)/iterations
        L.append(dt)
        print([G.R[x] for x in g])
        if n!=N:
            G = G.Cross(G0)
    L = array(L)
    L = log(L)
    plot([2**n for n in range(1,N+1)],L)
    show()
    exit()
    G = DihegralGroup(3)
    print("Relation :")
    print(G.R,'\n')
    g = G.MinimumGeneratingSet(True)
    print("Output :",g)
    if g is not None:
        print("Pretty Output :",",".join([str(G.R[i]) for i in g]))

    print('\n\n')

    G = PermutationGroup(4)
    print("Relation :")
    print(G.R)
    g = G.MinimumGeneratingSet(True)
    print("Output :",g)
    if g is not None:
        print("Pretty Output :",",".join([str(G.R[i]) for i in g]))

