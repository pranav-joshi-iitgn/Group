from numpy import *

class Relation:
    def __init__(self,T:ndarray,Elements=None,calculate=True) -> None:
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

    def __getitem__(self,i):
        return self.Elements[i]
    def __repr__(self) -> str:
        s = "\t  "
        for name in self.Elements:
            s += str(name) + '\t'
        s += "\n\t  "
        for name in self.Elements:
            s += '-'*len(str(name)) + '\t'
        i = 0
        for name in self.Elements:
            s += "\n"+str(name)+'\t'+"| "
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
    def generate(self,generators):
        T=self.T
        n = T.shape[0]
        visited = [False]*n
        visited[self.identity_ind] = True
        while True:
            changes = 0
            for x in range(n):
                if not visited[x]:
                    continue
                for g in generators:
                    gx = T[g][x]
                    if not visited[gx]:
                        visited[gx] = True 
                        changes += 1
            if not changes:
                break
        L = [i for i in range(n) if visited[i]]
        try:
            G = Group(self,L)
        except:
            raise TypeError("Can't make a group from"+str(L)+"generated from"+str(generators))
        return G

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

def get_identity(R:Relation,ElInd=None):
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
    
def is_associative(R:Relation,ElInd=None):
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

def get_inverses(R:Relation,identity_ind,ElInd=None)->list:
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

def contains_duplicates(L):
    n = len(L)
    for i in range(n):
        for j in range(i):
            if L[i] == L[j]:
                return True
    return False

def is_commutative(R:Relation,ElInd=None):
    T = R.T
    n = T.shape[0]
    if ElInd is None:
        ElInd = range(n)
    Ts = T[ElInd][:,ElInd]
    return array_equal(Ts.T,Ts)

class Group:
    def __init__(self,R:Relation,ElInd=None,check=True) -> None:
        T = R.T
        n = T.shape[0]
        if ElInd is None:
            ElInd = list(range(n))
        self.ElInd = ElInd
        self.R = R
        if not check:
            return
        id = get_identity(R,ElInd)
        if id is None:
            raise TypeError("The relation doesn't have an identity element")
        if not is_associative(R,ElInd):
            raise TypeError("The relation provided isn't associative")
        inv = get_inverses(R,id,ElInd)
        assert len(inv) == len(ElInd)
        if None in inv:
            raise TypeError("An element doesn't have it's inverse.")
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
        return "Group("+",".join([str(self.R.Elements[i]) for i in self.ElInd])+")"
    def Element(self,i):
        return Element(self.R,i)
    def has_subgroup(self,H,check=False):
        if not isinstance(H,Group):
            return False
        if self.R is not H.R:
            return False
        for h in H.ElInd:
            if h not in self.ElInd:
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
        Rnew = Relation(M,cosets)
        return Group(Rnew)
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
        xtoG = R.generate(L)
        return xtoG
    def MinimumNormalSubGroup(self):
        R = self.R
        T = R.T
        n = T.shape[0]
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
        
    def MinimalGeneratingSet(self)->list:
        R = self.R
        T=R.T
        n = T.shape[0]
        visited = [False]*n
        visited[self.identity_ind] = True
        generators = []
        for g in self.ElInd:
            if visited[g]:
                continue
            Changes = 0
            while True:
                changes = 0
                for x in self.ElInd:
                    if not visited[x]:
                        continue
                    gx = T[g][x]
                    if not visited[gx]:
                        visited[gx] = True 
                        changes += 1
                Changes += changes
                if not changes:
                    break
            if Changes:
                generators.append(g)
            done = True
            for x in self.ElInd:
                if not visited[x]:
                    done = False
                    break
            if done:
                break
        return generators
        
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

    def MinimumGeneratingSet(self,debug=False):
        if debug:
            print("Finding minimum generating set for",self)
        T = self.R.T
        N = self.MinimumNormalSubGroup()
        if debug:
            print("N :",N)
        if self==N:
            return self.MinimalGeneratingSet()
        n = N.MinimalGeneratingSet()
        if debug:
            print("n :",n)
        m = len(n)
        GbyN = G/N
        if debug:
            print("G/N :",GbyN)
        mingenGbyN = GbyN.MinimumGeneratingSet(debug)
        g = [(GbyN[i]).g_ind for i in mingenGbyN]
        if debug:
            print("g :",g)
        l = len(g)
        if N.is_abelian():
            if debug:
                print("N is abelian.")
            proceed = not G.has_generating_set(g)
            for i in range(l):
                for j in range(m):
                    modified_g = g[:i]+[T[g[i]][n[j]]]+g[i+1:]
                    proceed = proceed and not G.has_generating_set(modified_g)
            if proceed:
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
                    for nl in N:
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
        return str(g)+" * "+str(H)
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

def AdditiveGroupOnIntegersModulo(n):
    T = array([
        [(i+j)%n for j in range(n)]
        for i in range(n)
    ])
    names = [str(i) for i in range(n)]
    R = Relation(T,names)
    return Group(R)

def DihegralGroup(n):
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

if __name__=='__main__':
    G = DihegralGroup(3)
    print("Relation :")
    print(G.R,'\n')
    g = G.MinimumGeneratingSet(True)
    print("Output :",g)
    print("Pretty Output :",[G.R[i] for i in g])
