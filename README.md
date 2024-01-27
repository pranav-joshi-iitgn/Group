This is a coustom made group theory library to implement the first method to find the minumum generating set of a group in polynomial time. 

Output upon running the test case:
```
Relation :
          e     r^1     r^2     s       sr^1   sr^2
          -     ---     ---     -       ----   ----
e       | e     r^1     r^2     s       sr^1   sr^2
r^1     | r^1   r^2     e       sr^2    s   sr^1
r^2     | r^2   e       r^1     sr^1    sr^2   s
s       | s     sr^1    sr^2    e       r^1   r^2
sr^1    | sr^1  sr^2    s       r^2     e   r^1
sr^2    | sr^2  s       sr^1    r^1     r^2   e        

Finding minimum generating set for Group(e,r^1,r^2,s,sr^1,sr^2)
N : Group(e,r^1,r^2)
n : [1]
G/N : Group(e * Group(e,r^1,r^2),s * Group(e,r^1,r^2))
Finding minimum generating set for Group(e * Group(e,r^1,r^2),s * Group(e,r^1,r^2))
N : Group(e * Group(e,r^1,r^2),s * Group(e,r^1,r^2))
g : [3]
N is abelian.
Output : [3, 1]
Pretty Output : s,r^1
```
