<p style="text-align: center; font-size: 32px; font-weight:bold;">
  Graph Coloring Crossover Model
</p>



[TOC]



# Graph Coloring Crossover

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                                                         |
| ---- | ---------------------- | ----------- | ------------ | ------------------------------------------------------------ |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ |               |
| $C$ | **c**olor set | $[2, 200]$ | $n$ |               |
| $P$ | **p**arent set | $[2, 20]$ | $p$ |               |

### Constant

| Constant | Description                    | Type | Range       | Remark |
| -------- | ------------------------------ | ---- | ----------- | ------ |
| $A_{pnc}$ | node $n$ is **a**ssigned with color $c$ in parent $p$ | bool | $\{0, 1\}$ |  |


## Decision

| Variable     | Description                                                | Type | Domain     | Remark                                                     |
| -------------- | ------------------------------------------------------------ | ---- | ------------- | ------------------------------------------------------------ |
| $s_{pc}$ | node set **s**elected | bool | $\{0, 1\}$ |  |
| $l_{n}$ | covering **l**evel of node $n$ | real | $[0, 1]$ | $l_{n} = 1$ means fully covered |


## Objective

### Maximize the Node Coverage **ONC (node coverage)**

cover nodes as many as possible.

$$
\max \sum_{n \in N} ...
$$


## Constraint

all of the following constraints must be satisfied.

- **HNC (node coverage)** ...
  $$
  ...
  $$



# Graph Coloring Boolean Decision Model

## Known

### Set

| Set | Description                   | Size        | Element         | Remark                                                         |
| ---- | ---------------------- | ----------- | ------------ | ------------------------------------------------------------ |
| $N$ | **n**ode set | $[2, 2000]$ | $n, m$ |               |
| $E$ | **e**dge set | $[1, |N|^2]$ | $e, (n, m)$ |               |
| $C$ | **c**olor set | $[2, 200]$ | $c$ |               |


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
| $C$ | **c**olor set | $[2, 200]$ | $c$ |               |


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
| $C$ | **c**olor set | $[2, 200]$ | $c$ |               |


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
