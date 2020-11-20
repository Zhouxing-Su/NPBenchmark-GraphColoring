<p style="text-align: center; font-size: 32px; font-weight:bold;">
  Graph Coloring Models
</p>



[TOC]



# Graph Coloring Clustering Optimization Model

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                             |
| ---- | ---------------------- | ----------- | ------------ | ----------------- |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ |               |
| $E$ | **e**dge set | $[1, |N|^2]$ | $e, (n, m)$ |               |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |


## Decision

| Variable     | Description                        | Type | Domain     | Remark          |
| -------------- | -------------------------- | ---- | ------------- | --------------------- |
| $x_{nm}$ | nodes $n$ and $m$ are assigned to the same color | bool | $\{0, 1\}$ | $(n, m) \notin E, m \le n$ (lower triangle) |
| $y_{n}$ | node $n$ is the agent of its belonging cluster | bool | $\{0, 1\}$ | a color is used |


## Objective

### Minimize the Used Colors **OUC (used colors)**

the number of used colors.

$$
\min \sum_{n \in N} y_{n}
$$


## Constraint

all of the following constraints must be satisfied.
if $x_{nm}$ exists in an expression, $(n, m) \notin E$ and $m \le n$ must be respected.

- **HCB (cluster belonging)** each node belongs to a cluster whose containing nodes are assigned the same color. the agent of a cluster is the node with the smallest index.
  $$
  \sum_{m \in N} x_{nm} = 1, \quad \forall n \in N
  $$
  - it can be relaxed since belonging to more than one cluster never leads to better optima.
    $$
    \sum_{m \in N} x_{nm} \ge 1, \quad \forall n \in N
    $$

- **HAN (agent node)** node $m$ is the agent of cluster in which the nodes with the same color.
  $$
  \sum_{n \in N} x_{nm} \le |N| \cdot y_{m}, \quad \forall m \in N
  $$
  - it can be tighten if $n$ and $m$ are zero-based index, where $(|N| - m)$ is the upper bound of the summing term.
    $$
    \sum_{n \in N} x_{nm} \le (|N| - m) \cdot y_{m}, \quad \forall m \in N
    $$

- **HCP1 (conflict propagation)** adjacent nodes cannot share the same agent, i.e., be in the same cluster.
  $$
  x_{nm} + x_{n'm} \le 1, \quad \forall n, n', m \in N, (n, n') \in E
  $$



# Graph Coloring Crossover Decision Model

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                                                         |
| ---- | ---------------------- | ----------- | ------------ | ------------------------------------------------------------ |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ | $I_{n} \subset N$ is the set of independent sets containing node $n$ |
| $P$ | **p**arent set | $[2, 2000]$ | $p$ | $I_{p} \subset I$ is a set of related independent sets |
| $I$ | pseudo **i**ndependent sets | $[2, 2000]$ | $i$ | $N_{i} \subset N$ is a set of independent nodes |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |

### Constant

| Constant | Description                    | Type | Range       | Remark |
| -------- | ------------------------------ | ---- | ----------- | ------ |
| $S^{-}_{p}$ | min number of independent sets **s**elected from parent $p$ | int | $[0, |I_{p}|]$ |  |
| $S^{+}_{p}$ | max number of independent sets **s**elected from parent $p$ | int | $[0, |I_{p}|]$ |  |
| $L_{ni}$ | covering **l**evel of node $n$ by independent set $i$ | real | $[0, 1]$ | $L_{in} < 1$ if adjacent nodes of $n$ in $i$ exists |


## Decision

| Variable     | Description                                                | Type | Domain     | Remark                                                     |
| -------------- | ------------------------------------------------------------ | ---- | ------------- | ------------------------------------------------------------ |
| $s_{i}$ | independent node set $i$ is **s**elected | bool | $\{0, 1\}$ |  |
| $l_{n}$ | covering **l**evel of node $n$ | real | $[0, 1]$ | $l_{n} = 1$ means fully covered without conflict |
| $y_{ni}$ | covering node $n$ by independent set $i$ | real | $[0, 1]$ | auxiliary variable for linearization |


## Objective

### Maximize the Node Coverage **ONC (node coverage)**

cover nodes as many as possible.

$$
\max \sum_{n \in N} l_{n}
$$


## Constraint

all of the following constraints must be satisfied.

- **HNC (node coverage)** every node belongs to an independent set which will be assigned with the same color.
  use this if $\sum_{i} L_{in} < 1$ where $i$ belongs to the independent sets with conflicts on node $n$, otherwise, better use **HCL (coverage level)**.
  $$
  \sum_{i \in I_{n}} L_{ni} \cdot s_{i} \ge l_{n}, \quad \forall n \in N
  $$

- **HCL (coverage level)** every node belongs to an independent set which will be assigned with the same color.
  **HNC (node coverage)** is trivial when this constraint is enabled.
  $$
  \max_{i \in I_{n}} \{L_{ni} \cdot s_{i}\} \ge l_{n}, \quad \forall n \in N
  $$
  which can be linearized as
  $$
  y_{ni} + L_{ni} \cdot s_{i} \ge l_{n}, \quad \forall n \in N, \forall i \in I_{n}
  $$
  $$
  \sum_{i \in I_{n}} y_{ni} \le |I_{n}| - 1, \quad \forall n \in N
  $$

- **HSN (set number)** given number of independent sets will be selected.
  the $\ge$ can be replaced with $=$ to tighten the bound, but it may not improve the performance.
  $$
  \sum_{i \in I} s_{i} \ge |C|
  $$

- **HBS (balance selection)** the number of independent sets selected from the same parent is limited.
  $$
  S^{-}_{i} \le \sum_{i \in I_{p}} s_{i} \le S^{+}_{i}, \quad \forall p \in P
  $$



# Graph Coloring Column Generation Decision Model

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                        |
| ---- | ---------------------- | ----------- | ------------ | ------------------------- |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ | $I_{n} \subset N$ is the set of independent nodes containing node $n$ |
| $I$ | **i**ndependent sets | $[2, 2000]$ | $i$ | $N_{i} \subset N$ is a set of independent nodes |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |


## Decision

| Variable     | Description                                                | Type | Domain     | Remark            |
| -------------- | ------------------------------------------------------------ | ---- | ------------- | ------------- |
| $s_{i}$ | independent node set $i$ is **s**elected | bool | $\{0, 1\}$ |  |


## Constraint

all of the following constraints must be satisfied.

- **HNC (node coverage)** every node belongs to an independent set which will be assigned with the same color.
  $$
  \sum_{i \in I_{n}} s_{i} \ge 1, \quad \forall n \in N
  $$

- **HSN (set number)** given number of independent sets will be selected.
  the $\ge$ can be replaced with $=$ to tighten the bound, but it may not improve the performance.
  $$
  \sum_{i \in I} s_{i} \ge |C|
  $$



# Graph Coloring Boolean Decision Model

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                      |
| ---- | ---------------------- | ----------- | ------------ | --------------------- |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ |               |
| $E$ | **e**dge set | $[1, |N|^2]$ | $e, (n, m)$ |               |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |


## Decision

| Variable     | Description                       | Type | Domain     | Remark      |
| -------------- | ------------------------------ | ---- | ------------- | ------------------ |
| $x_{nc}$ | assign color $c$ to node $n$ | bool | $\{0, 1\}$ |  |

### Convention and Function

- define function $\textrm{A}(n)$ to be the set of adjacent nodes of node $n$.
- define function $\textrm{I}(N')$ to be the max independent set in node set $N'$.
- define function $\textrm{C}(N')$ to be the max clique in node set $N'$.
- define function $\textrm{G}(N', i)$ to get the $i^{\textrm{th}}$ node in node set $N'$.


## Constraint

all of the following constraints must be satisfied.

- **HSC (single color)** assign single color to each node only.
  $$
  \sum_{c \in C} x_{nc} = 1, \quad \forall n \in N
  $$

- **HCA (conflict avoidance)** the colors assigned to a pair of end nodes on single edge should be different.
  $$
  x_{nc} + x_{mc} \le 1, \quad \forall (n, m) \in E, \forall c \in C
  $$
  it can be reduced by combining the constraints $\forall n \in \textrm{A}(n)$.
  $$
  |\textrm{A}(n)| \cdot x_{nc} + \sum_{m \in \textrm{A}(n)} x_{mc} \le |\textrm{A}(n)|, \quad \forall n \in N, \forall c \in C
  $$
  it can be strengthened by the fact that the number of nodes with the same color never exceeds the size of the max independent set.
  $$
  |\textrm{I}(\textrm{A}(n))| \cdot x_{nc} + \sum_{m \in \textrm{A}(n)} x_{mc} \le |\textrm{I}(\textrm{A}(n))|, \quad \forall n \in N, \forall c \in C
  $$

- **HSO.O (size order)** the number of nodes assigned to certain color decreases monotonically as the color index increases.
  i.e., the number of nodes assigned to smaller color index should not exceed the number of nodes assigned to greater color index.
  it conflicts with **HEO.O (exploitation order)**.
  $$
  \sum_{n \in N} x_{nc'} \ge \sum_{n \in [0, n)} x_{nc}, \quad \forall c, c' \in C, c' = c - 1
  $$

- **HIO.O (index order)** each node cannot be assigned to the colors whose indices exceed the node's.
  $$
  x_{nc} = 0, \quad \forall n \in N, c \in C, n < c
  $$

- **HEO.O (exploitation order)** the min node index assigned to certain color increases monotonically as the color index increases.
  i.e., greater color index cannot be used if the smaller color index is not used by nodes with smaller index.
  it conflicts with **HSO.O (size order)**.
  $$
  \sum_{m \in [0, n)} x_{mc'} \ge x_{nc}, \quad \forall n \in N, \forall c, c' \in C, c' = c - 1
  $$
  if **HIO.O (index order)** is enabled, it can be reduced by substituting the 0-value items.
    $$
    \sum_{m \in [c', n)} x_{mc'} \ge x_{nc}, \quad \forall n \in N, \forall c, c' \in C, c' = c - 1, n \ge c
    $$

- **HFC.O (fixed clique)** the nodes in a clique should be assigned to different colors.
  if $\textrm{C}(N) \neq \{0, 1, 2, ..., |\textrm{C}(N)| - 1\}$, it conflicts with **HIO.O (index order)** and **HEO.O (exploitation order)**.
  $$
  x_{nc} = 1, \quad \forall n \in \textrm{C}(N), \forall c \in C, n = \textrm{G}(\textrm{C}(N), c)
  $$
  $$
  x_{nc} = 0, \quad \forall n \in \textrm{C}(N), \forall c \in C, n \neq \textrm{G}(\textrm{C}(N), c)
  $$



# Graph Coloring Relaxed Boolean Decision Model

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                             |
| ---- | ---------------------- | ----------- | ------------ | ----------------- |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ |               |
| $E$ | **e**dge set | $[1, |N|^2]$ | $e, (n, m)$ |               |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |


## Decision

| Variable     | Description                        | Type | Domain     | Remark          |
| -------------- | -------------------------- | ---- | ------------- | --------------------- |
| $x_{nc}$ | assign color $c$ to node $n$ | bool | $\{0, 1\}$ |  |
| $y_{e}$ | there is a conflict on edge $e$ | real | $[0, 1]$ |  |


## Objective

### Minimize the Conflict Number **OCN (conflict number)**

the number of conflicting edges.

$$
\min \sum_{e \in E} y_{e}
$$


## Constraint

all of the following constraints must be satisfied.

- **HSC (single color)** assign single color to each node only.
  $$
  \sum_{c \in C} x_{nc} = 1, \quad \forall n \in N
  $$

- **HCA (conflict avoidance)** the colors assigned to a pair of end nodes on single edge should be different.
  $$
  x_{nc} + x_{mc} \le 1 + y_{e}, \quad \forall e = (n, m) \in E, \forall c \in C
  $$



# Graph Coloring Integer Decision Model

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                 |
| ---- | ---------------------- | ----------- | ------------ | -------------------- |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ |               |
| $E$ | **e**dge set | $[1, |N|^2]$ | $e, (n, m)$ |               |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |


## Decision

| Variable     | Description                           | Type | Domain     | Remark                  |
| -------------- | ------------------------------------ | ---- | ------------- | ------------------- |
| $x_{n}$ | the color assigned to node $n$ | int | $C$ |  |
| $z_{e}, z_{nm}$ | the color index for node $n$ is less than the one for node $m$ | real | $[0, 1]$ | $e = (n, m) \in E$ |


## Constraint

all of the following constraints must be satisfied.

- **HCA (conflict avoidance)** the colors assigned to a pair of end nodes on single edge should be different.
  $$
  1 - |C| \cdot (1 - z_{e}) \le x_{n} - x_{m} \le |C| \cdot z_{e} - 1, \quad \forall e = (n, m) \in E
  $$



# Notation

there are several notations about the constraints.

- hard constraints
  - prefixed with **H** in the tag
  - constraints in the original problem that must be satisfied
- soft constraints
  - prefixed with **S** in the tag
  - objectives in the original problem
  - cut-off scores to the low-priority objectives in case the high-priority objectives are over dominant
- auxiliary constraints
  - prefixed with **A** in the tag
  - additional constraints for converting non-linear constraints into linear ones
- optional constraints
  - prefixed with **.O** in the tag
  - user cuts or constraints that are very likely to be ignored when changing requirement or implementation
- lazy constraints
  - suffixed with **.L** in the tag
  - constraints which are recommended to be added lazily
