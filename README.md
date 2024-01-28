This is a coustom made group theory library to implement the first method to find the minumum generating set of a group in polynomial time. 

Output upon running the test case:
```
Relation :
           e    r^1     r^2     s       sr^1    sr^2
           -    ---     ---     -       ----    ----
e       |  e    r^1     r^2     s       sr^1    sr^2
r^1     |  r^1  r^2     e       sr^2    s       sr^1
r^2     |  r^2  e       r^1     sr^1    sr^2    s
s       |  s    sr^1    sr^2    e       r^1     r^2
sr^1    |  sr^1 sr^2    s       r^2     e       r^1
sr^2    |  sr^2 s       sr^1    r^1     r^2     e        

Finding minimum generating set for :
 Group(e,r^1,r^2,s,sr^1,sr^2)
N :
 Group(e,r^1,r^2)
n :
 [1]
G/N :
 Group((e * Group(e,r^1,r^2)),(s * Group(e,r^1,r^2)))
Finding minimum generating set for :
 Group((e * Group(e,r^1,r^2)),(s * Group(e,r^1,r^2)))
N :
 Group((e * Group(e,r^1,r^2)),(s * Group(e,r^1,r^2)))
g :
 [3]
N is abelian.
Output : [3, 1]
Pretty Output : s,r^1



Relation :
Relation([0123],[0132],[0213],[0231],[0312],[0321],[1023],[1032],[1203],[1230],[1302],[1320],[2013],[2031],[2103],[2130],[2301],[2310],[3012],[3021],[3102],[3120],[3201],[3210])
Finding minimum generating set for :
 Group(
[0123],
[0132],
[0213],
[0231],
[0312],
[0321],
[1023],
[1032],
[1203],
[1230],
[1302],
[1320],
[2013],
[2031],
[2103],
[2130],
[2301],
[2310],
[3012],
[3021],
[3102],
[3120],
[3201],
[3210])
N :
 Group([0123],[1032],[2301],[3210])
n :
 [7, 16]
G/N :
 Group(
([0123] * Group[level=1]),
([0132] * Group[level=1]),
([0213] * Group[level=1]),
([0231] * Group[level=1]),
([0312] * Group[level=1]),
([0321] * Group[level=1]))
Finding minimum generating set for :
 Group(
([0123] * Group[level=1]),
([0132] * Group[level=1]),
([0213] * Group[level=1]),
([0231] * Group[level=1]),
([0312] * Group[level=1]),
([0321] * Group[level=1]))
N :
 Group(([0123] * Group[level=1]),([0231] * Group[level=1]),([0312] * Group[level=1]))
n :
 [3]
G/N :
 Group((([0123] * Group[level=1]) * Group[level=2]),(([0132] * Group[level=1]) * Group[level=2]))
Finding minimum generating set for :
 Group((([0123] * Group[level=1]) * Group[level=2]),(([0132] * Group[level=1]) * Group[level=2]))
N :
 Group((([0123] * Group[level=1]) * Group[level=2]),(([0132] * Group[level=1]) * Group[level=2]))
g :
 [1]
N is abelian.
g :
 [1, 3]
N is abelian.
Output : [6, 3]
Pretty Output : [1023],[0231]
```
