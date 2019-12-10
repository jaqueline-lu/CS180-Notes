## Divide and Conquer (chap.5)
### Counting Inversions

### Closest Pair of Points
1. Sort every point in increasing x-coordinate and increasing y coordinate, put it as Px and Py
2. Let L represents the first n/2 positions in Px and R represents the second n/2 positions in Px
3. Pass through Px and Py, we create: Lx including the points in L sorted in increasing x coordinate, Ly including points in L sorted in increasing y coordinate; and so do it for R
4. Recursively find the closest pair of points in L and R with accessing to Lx and Ly, Rx and Ry respectively. Then we have the closest pair of points from L as l1 and l2, and r1 and r2 respectively
5. Get the minimum from d(l1,l2) and d(r1,r2) as d
6. Find the x* as the x coordinate of the rightmost point in L, and we let V represents the vertical line by x=x* which separate L from R
7. Let ${S \subseteq P}$ this set, and let Sy denote the list including sorted increasing y coordinate list points
8. Single pass through Py, we get Sy
9. single pass through Sy, and for each s in Sy we compute its distance to each of the next 15 points
10. Compare the distance from the closest pair in the middle and on the side

Proof:
- if there exists ${l \in L}$ and ${r \in R}$ and ${d(l,r) < \sigma}$ then each of l and r lies withn a distance ${\sigma}$ of V
- if ${s,s' \in S}$ have the proprety that d(s,s') < ${\sigma}$ then s and s' are within 15 positions of each other in the sortd list Sy

Runtime:
O(nlogn)
## Dynamic Programming (chap.6)
1. There are only a polynomial number of subproblems
2. the solution to the original problem can be computed from the solutions to the subproblems
3. a natural ordering on subproblem from "smallest" to "largest"

### Weighted Interval Scheduling
Recurrence: $Opt(j)={max(v_j+Opt(p(j)),Opt(j-1)))}$
base case: when j=0, opt is 0.

Runtime:
O(n): each time the procedure invodes the recurrence, giving 2 recursive calls to Opt(), it fills in a new entry. Since M has n+1 entries, there can be at most O(n) calls to Opt(n).
Takes O(logn)time to sort, so total **O(nlogn)**.

if wants to look for a solution other than just value:
algorithm:
```python
Find Solution(j):
if j=0 then output nothing
elif vj+M[p(j)]>= M[j-1] :
  output j with result of find-solution(p(j))
else output result of find-solution(j-1)
```

Runtime: O(nlogn)

### Knapsacks, Subset Sums
Suppose the weight of the knapsacks it can holds W and the number of item is n
if ${n \notin \alpha}$ then ${Opt(n,W)=Opt(n-1,W)}$
if ${n \in \alpha}$ then ${Opt(n,W)=w_n+Opt(n-1,W-w_n)}$ 
recurrence: ${Opt(i,w)=max(Opt(i-1,w),w_i+Opt(i-1,w-w_i))}$
```python
Given an array M[0...n,0...W]
Initialize M[0,w] =0 for each w=0,1...w #base case
For i=1,2..n:
  For w =0,...W:
    Use recurrence for M[i,w]
Return M[n,W]
```
runtime:
O(nW) time. pseudo-polynoial, can be efficient when numbers {wi} are small, but become less practical when w is large.

### RNA Secondary Structure
- no sharp turns
- {A,U} or {C,G}
- no base appears in more than 1 pair
- no crossing

Maximum number of matching
Recurrence: ${Opt(i,j)=max(Opt(i,j-1),max(1+Opt(i,t-1)+Opt(t+1,j-1)))}$ assuming t and j are matching
```python
Initialize Opt(i,j)=0 whenever i>=j-4
for k=5,6..n-1:
  for i =1,2,...n-k:
    set j=i+k
    Use recurrence for Opt(i,j)
Return Opt(1,n)
```
Runtime:
${O(n^2)}$ subproblems to solve and the recurrence takes O(n) for each. Total running time is ${O(n^3)}$
### Sequence Alignment
optimal alignment M:
- ${(m,n)\in M}$ or
- the mth position of X is not matched
- the nth position of Y is not matched
Recurrence:
${LCS(i,j)=max(LCS(i-1,j-1)+1,LCS(i-1,j),LCS(i,j-1))}$

Proof by contradiction:
if the last a in X is matching with the a before in Y:
  - if we change the crossing, then we can change it to point at the last one and get the maximum

Runtime:
${O(|x||y|)}$ time, it's polynomial since |x| is given and |Y| is given

### Shortest Paths in a Graph
if G has no negative cycles, then there is a shortest path from s to t that is simple (no repeat nodes) and hence has at most n-1 edges
  - since every cycle has nonnegative cost, repeating will have larger cost
Recurrence:
${SP_l(s,x)=min_{1\le j\le d}(SP_{l-1}(s,y_j)+w_j)}$
Go through all the x and all the l length of edges included in the path, there will be n of such neighbors like y_j to x. 
```python
Define M[0,t]=0 and M[0,v]=infinity for all other v in V
For i =1...n-1
  for v in V that w has an edge to:
    compute M[i,v] using the recurrence
return M[n-1,s]
```
runtime:
If we check whether there is an edge between two nodes, and there are only m edges in total and n nodes, the runtime will be O(MN).
The worst case is ${O(n^3)}$
## Network Flow (chap.7)
### Maximum-Flow Problem & Ford-Fulkerson Algorithm
- There is at least one edge incident to each node
- All capacities are integers

```python
Ford-Fulkerson Algorithm:
Initialize flow for each e in G as 0
While there is an s-t path in the residual graph G: 
  Let P be a simple s-t path in G 
  f'= augment (f,P)
  Update f to be f'
  Update the residual graph G_f to be G_f'
Return f
```
Runtime:
C iterations, m>=n/2 edges,using BFS/DFS, finding an s-t path is O(m+n)=O(m)time. The augment(f,p)takes O(n)time, since there are n-1 edges in path P. 

### Maximum Flows and Min-Cuts
- Let f be any s-t flow, and (A,B) any s-t cut, then flow of some vertex v is smaller than cut of (A,B). ${v(f)\le(A,B)}$. The value of every flow is upper-bounded by the capacity of every cut.

#### Proof
1. f is a max flow N
2. The residual network ${N_f}$ has no augmenting path (which is the stopping condition for Ford-Fulkerson)
3. |f| = C(S,T) for some cut (s,t) of N

1->2: 
- Given a graph with a max flow f. If there is still an augmented path, then we can increase the flow by at least 1. 
- But this means that f is not maximum flow, which is contradict to the assumption. Done.

2->3:
- Assume that the residual network has no augmenting path, need to show that |f| is equal to C(s,t), or show the flow is capacity of some cut.
- Find the saturated edge, whose load is equal to its capacity. Remove the saturated edges. There will be no augmented path from s to t after removing.
- Since now s is separated from t, we put everything that has a path to s to s side, and do the same for t.
- Saturated edge's capacity is the flow. And the cut will have the number as the sum of the capacity of the saturated edges, which is the same as the sum of the flow.

3->1:
- Assume |f|=C(s,t), then f is a max flow.
- For any flow and cut, because all the flow must go through this cut(it makes the whole graph disconnected), so |f| must be no greater than cut.
- Then if |f| is already equal to cut, then f is already the largest possible, so it's the maximum flow.

### Choosing Good Augmenting Paths
Main idea: prioritize taking edges with larger capacites to avoid ending up with a path that has a small bottleneck.

Algorithm:
- Let U be the largest edge capacity in initial flow graph
- Let a be the largest power of 2 less than or equal to U
- Repeatedly find augmenting paths with remaining capacity ${\ge}$ a until no more paths satisfy
- Decrease the value of a by dividing it by 2 (a/=2) and repeat while a > 0

```python
Scaling Max-Flow:
Initialize f(e)=0 for all edge in G
Initialize p to be the largest power of 2 that is no larger than the maximum capacity out of s
while (p is no smaller than 1):
  while there is an s-t path in the graph G_f(p):
    Let A be a simple s-t path in G_f(p)
    f'=augment(f,P)
    Update f to be f' and update G_f(p)
  p = p/2
return f
```

Runtime:
using DFS, ${O(m^2log_2(U))}$
- there are m edges, and we can only use each of the edges out of s only for at most one augmentation in that phase. Since p'= 2p when we change from one phase to another, the maximum flow f' is at most ${v(f')\le v(f_p)+mp'=v(f_p)+2mp}$. 
- And find a maximum flow in at most ${2m(1+log_2C)}$ augmentations.
### Bipartite Matching
Begin with a graph G.
Direct all edges in G from X to Y. Then add a node s, and an edge (s,x) from s to each node in X. Then add a node t, and an edge (y,t) from each node in Y to t. Then give each edge in G a capacity of 1.

Runtime:
Using Ford-Fulkerson algorithm, the maximum matching in a bipartite graph is O(mn)time.

### Disjoint Paths in Directed and Undirected Graphs
Question: Given a directed graph G=(V,E) with two distinguished nodes s,t ${\in}$ V, the directed edge-disjoint paths problem is to find the maximum number of edge-disjoint s-t paths in G. The Undirected Edge-disjoint paths is to find the maximum number of edge-disjoint s-t path in an undirected graph G.
Algorithm:
- Suppose there are k edge-disjoint s-t paths. Make each of these paths carry 1 unit of flow.
- We set the flow to be f(e)=1 for each edge e on any of the paths, and f(e')=0 on all other edges.
- If threre are k edge-disjoint paths in a directed graph G from s to t, then the value of the maximum s-t flow in G is at least k.

Proof:
Show if there is a flow of value k, then there exist k edge-joint s-t paths.
- if f is a 0-1 valued flow of value v, then the set of edges with flow value f(e)=1 contains a set of v edge-disjoint paths.
  - base case: if v=0,trivial.
  - since (s,u) carries a unit of flow, then there is some edge (u,v) carries one unit of flow and so forth. And by doing this, we either reach t, or we reach a node v for a second time.
  - if we reach t, then we find a path from s to t.

> There are k edge-disjoint paths in a directed graph G from s to t iff the maximum value of an s-t flow in G is at least k.

Runtime:
Ford-Fulkerson gives maximum set of edge-disjoint s-t path in a directed graph G in O(mn) time.

### Survey Design
### Airline Scheduling
- for each flight i, the graph G will have two nodes ui and vi. Set up a source node s and sink node t
- 
## NP-Completeness (chap.8)
### Polynomial-Time Reductions
- if problem Y can be solved using a polynomial number of standard computational steps, plus a polynomial number of calls to some black box that solves problem X, then Problem Y is polynomial-time reducible to X, which written as: ${Y\le_p X}$
- Suppose ${Y\le_p X}$, if X can be solved in poly. time, then Y can be solved in poly. time
- Suppose ${Y\le_p X}$, if Y cannot be solved in polynomial time, then X cannot be solved in polynomial time (**Suppose ${Y\le_p X}$, if Y is an NP-Complete problem, then X is NP-Complete**)

#### Maximum Independent Set and Vertex over
- a set of nodes ${S\subseteq V}$ is **independent** if no two nodes in S are joined by an edge.
- a set of nodes ${S\subseteq V}$ is **vertex cover** if every edge e ${\in}$ E has at least one end in S.
- Let G=(V,E) be a graph, then S is an independent set iff ifs complement V-S is a vertex cover
- **Independent Set ${\le_p}$ Vertex Cover and Vertex Cover ${\le_p}$ Independent Set**

#### 3-SAT
- Given a set of clauses C1,...Ck, each of length 3, over a set of variables X={x1,...xn}, does threre exist a satisfying truth assignment?
- **3-SAT ${\le_p}$ Independent Set**
  - choose one term from each clauses, then find a truth assignment that causes all these terms to evaluate to 1, so as to satisfying all clauses. 
  - Succeed if you can select a term from each clause so that no two selected terms in the way that "one is equal to a variable xi and the other is equal to its complement ${\bar{x}}$"
  - Construct G=(V,E) consisting of 3k nodes grouped into k triangles. For i=1,2,...k we construct 3 verices ${v_{i1}v_{i2},v_{i3}}$ edge one another. So ${v_{ij}}$ is term labeled with jth term from the clause ${C_i}$ of the 3-SAT
  - choosing 1 vertex from each triangle = choosing a term in each clause that evaluate to 1
  - for each pair of vertices whose labels refer to terms that conflict, add an edge between them.
- Proof
  - if 3-SAT is satisfiable, then each triangel in the graph contains at least one node that label evaluates to 1. We set S to be the set consisting of 1 such node from each triangle. Then S is independent because if there were an edge between 2 nodes in S then their labels would be conflict. But since they evaulate to 1, it's not possible
  - Conversely, 
#### Transformation
Proof Strategy:
1. Prove ${X\in NP}$
2. Choose a problem Y that is known to be NP-Complete
3. Think of an instance of problem Y and show how construct an instance of problem X in polynomial time so that:(instance X and Y have the same answer)
   - if instance of Y is yes, then instance of X is also yes
   - if instance of X is yes, then instance of Y is also yes

- if ${Z\le_p Y}$, and ${Y\le_p X}$, then ${Z\le_p X}$
- NP-Completen Problems: Travelingsalesman problem, Hamiltonian Cycle, 3D Matching, 
- 3-SAT ${\le_p}$ Hamiltonian Cycle
- Hamiltonian Cycle ${\le_p}$ Hamiltonian Path
- Hamiltonian Cycle ${\le_p}$ Traveling Salesman
- 3-SAT ${\le_p}$ 3D-Matching
- 3-SAT ${\le_p}$ 3-coloring
