<p style="text-align: center; font-size: 32px; font-weight:bold;">
  Graph Coloring Models
</p>



[TOC]



# Graph Coloring Crossover

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



# Graph Coloring Column Generation

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                                                         |
| ---- | ---------------------- | ----------- | ------------ | ------------------------------------------------------------ |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ | $I_{n} \subset N$ is the set of independent nodes containing node $n$ |
| $I$ | **i**ndependent sets | $[2, 2000]$ | $i$ | $N_{i} \subset N$ is a set of independent nodes |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |


## Decision

| Variable     | Description                                                | Type | Domain     | Remark                                                     |
| -------------- | ------------------------------------------------------------ | ---- | ------------- | ------------------------------------------------------------ |
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

| Set | Description                   | Size        | Element         | Remark                                                         |
| ---- | ---------------------- | ----------- | ------------ | ------------------------------------------------------------ |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ |               |
| $E$ | **e**dge set | $[1, |N|^2]$ | $e, (n, m)$ |               |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |


## Decision

| Variable     | Description                                                | Type | Domain     | Remark                                                     |
| -------------- | ------------------------------------------------------------ | ---- | ------------- | ------------------------------------------------------------ |
| $x_{nc}$ | assign color $c$ to node $n$ | bool | $\{0, 1\}$ |  |


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



# Graph Coloring Relaxed Boolean Decision Model

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                                                         |
| ---- | ---------------------- | ----------- | ------------ | ------------------------------------------------------------ |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ |               |
| $E$ | **e**dge set | $[1, |N|^2]$ | $e, (n, m)$ |               |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |


## Decision

| Variable     | Description                                                | Type | Domain     | Remark                                                     |
| -------------- | ------------------------------------------------------------ | ---- | ------------- | ------------------------------------------------------------ |
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

| Set | Description                   | Size        | Element         | Remark                                                         |
| ---- | ---------------------- | ----------- | ------------ | ------------------------------------------------------------ |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ |               |
| $E$ | **e**dge set | $[1, |N|^2]$ | $e, (n, m)$ |               |
| $C$ | **c**olor set | $[2, 2000]$ | $c$ |               |


## Decision

| Variable     | Description                                                | Type | Domain     | Remark                                                     |
| -------------- | ------------------------------------------------------------ | ---- | ------------- | ------------------------------------------------------------ |
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
