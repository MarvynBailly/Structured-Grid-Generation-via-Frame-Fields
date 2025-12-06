# Mixed-Integer Quadrangulation (MIQ)

## Objective
Describe the main objective and approach of this idea.

## Papers Referenced
- List relevant papers in the `papers/` directory

## Code Organization
- Describe the main scripts and their purposes

## Results
- Document key findings
- Reference images in the `images/` directory

## Notes

### Mixed-Integer Quadrangulation
- Process consists of two steps
  - cross field generation
    - can take sparse/dense orientation constraints
    - automatically places singularities which enable smoothest cross field
  - global parametrization
    - generate quad meshes satisfying various constraints
      - orientation
      - alignment 
      - integer singularity locations
    - An optional anisotropic stretch metric allows trade of squareness of the quads for improved feature alignment.
  - both steps can be reduced to Mixed-Integer Problem
    - A minimization problem where some of the variables are strictly integers
- Paper proposes a Greedy Mixed-Integer Solver (MIS)
  - The authors found traditional MIS, known as *direct rounding*, would get stuck at local minima
  - They propose a *greedy rounding* algorithm which is like Local Gauss-Seidel
    - The idea is basically to minimize and round the continuous variable that causes the smallest residual change
    - this variable is now assumed to be constant for the rest of the minimization
  - while this local method does not always work, they found it to be suitable for most cases
- The paper talks about how the search for the cross field takes into account the curvature of the triangulation
  - I'm a little worried that the results will be strange for strictly flat field but it should be okay
- Measuring Cross Field Smoothness
  - We can measure the smoothness via
  - $$E_\text{smooth} = \sum_{e_{ij} \in E} \left( \theta_i + \kappa_{ij} + \frac{\pi}{2} p_{ij} - \theta_j\right)^2$$
  - where $\theta_i$ wrt frame $j$ and $\kappa_{ij} \in (-\pi, \pi]$$. 

### Step 0 - Normalize Triangulation
- Turn any quads into triangles
- (current issue with the MD30P)

### Step 1 - Cross Field Generation
- Input: 
  - mesh $M = (V, E, F)$
  - $F_c \subset F$ with constrained directions $\theta_i = \hat{\theta}_i$
- Goal:
  - Minimize $E_\text{smooth}$
    - by finding integer $p_{ij}$ per edge and real $\theta_i$ per face
- Reduce the Search Space
  - to ensure a unique solution
    - fix one period jump per free triangle (e.g. zero) without changing the energy
    - no edge whose dual path connects two constrained faces should be fixed 
    - to construct a valid set of edges who period jump can be set to zero
      - use Dijkstra tress of the dual mesh
        - each constrained face $F_c$ is the root of a separate tree
        - on each step, each tree "conquers" a triangle
        - this results in tress with no closed loop and no intersections with other constrained faces
        - ![alt text](images/image.png)
          - the red arrows show the user specified constrained directions on faces $F_c$
          - the green forest show a set of valid edges to be set to zero without changing the energy
          - for each free face: $F \setminus F_c$, pick on of the valid edges and set to zero 
- Once a suitable set of edges to found:
  - for each edge in the set:
    - fix valid edge to zero
    - if we need a jump between constrained $f_i$ and $f_j$
      - $$p_{ij} = \text{round}\left( \frac{2}{\pi} \left( \hat{\theta}_j  - \hat{\theta}_i - \kappa_{ij} \right) \right).$$
- Now we wish to min $E_\text{smooth}$ which is a mixed-integer problem with: 
  - $| F \setminus F_C | \approx 2 |V|$ real valued variables
  - $|E|- | F \setminus F_C | \approx |V|$ integer valued variables $p_{ij}$.
- To use greedy-solver approach:
  - Taking the gradient of the $E_\text{smooth}$ and setting equal to zero gives the linear system:
    - $$\frac{\partial E_\text{smooth}}{\partial \theta_k}= \sum_{e_{kj} \in N(f_i)} 2 \left( \theta_k  + \kappa_{kj} + \frac{\pi}{2}p_{kj} - \theta_j \right) = 0$$
    - $$\frac{\partial E_\text{smooth}}{\partial p_{ij}}= \pi \left(\theta_i + \kappa_{ij} + \frac{\pi}{2} p_{ij} - \theta_j \right) = 0.$$
      - careful of antisymmetry on edges which lead to sign changes:
        - $p_{ij} = - p_{ji}$
        - $\kappa_{ij} = - \kappa_{ji}$
- after running the solver, the result:
  - smooth cross field 
  - integer valued period jumps correspond to singularities


### Step 2 - Parameterize Domain


### Step 3 - Extract Quadrangulation


### Greedy Mixed-Integer Solver
- decide to implement or find existing tool







---
Created: 2025-11-19
