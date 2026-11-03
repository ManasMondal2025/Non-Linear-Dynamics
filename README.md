# Three-Species Non-Delayed Food Chain Model

Analysis of a tritrophic (prey, predator, top predator) system of nonlinear ODEs, covering equilibria, local stability, numerical simulation, and the emergence of chaos as parameters are varied.

Based on Section 13 of the M.Sc. dissertation *Non-Linear Dynamics* (Dept. of Mathematics, RKMVERI). Figure numbers below refer to that dissertation.

## Contents

1. [Model](#1-model)
2. [Equilibrium points](#2-equilibrium-points)
3. [Stability analysis](#3-stability-analysis)
4. [Numerical simulations](#4-numerical-simulations)
5. [Applications and implications](#5-applications-and-implications)
6. [Key findings](#6-key-findings)

---

## 1. Model

A food chain with three populations:

| Variable | Population |
|---|---|
| `X` | Prey |
| `Y` | Predator |
| `Z` | Top predator (superpredator) |

$$\frac{dX}{dt} = X\left[a_1\left(1-\frac{X}{K}\right) - \frac{wY}{\frac{X^2}{j}+X+a}\right]$$

$$\frac{dY}{dt} = Y\left[-a_2 + \frac{w_1X}{\frac{X^2}{j}+X+a} - \frac{w_2Z}{Y+D_2+bZ}\right]$$

$$\frac{dZ}{dt} = Z\left[-c + \frac{w_3Y}{Y+D_2+bZ}\right]$$

Parameters: `a1, K, w, j, a, a2, w1, w2, D2, b, c, w3`. Reading off the equations, `a1` is the prey growth rate and `K` its carrying capacity, while `a2` and `c` are the loss terms of the predator and top predator.

---

## 2. Equilibrium points

| Point | Coordinates | Notes |
|---|---|---|
| `E0` | `(0, 0, 0)` | Trivial equilibrium |
| `E1` | `(K, 0, 0)` | Prey-only equilibrium |
| `E2` | `(x2, y2, 0)` | Top predator absent |
| `E3` | `(x3, y3, z3)` | Coexistence equilibrium |

**E2**

$$x_2 = \frac{j}{2}\left[\left(\frac{w_1}{a_2}-1\right) \pm \sqrt{\left(\frac{w_1}{a_2}-1\right)^2-\frac{4a}{j}}\right],\qquad y_2 = \frac{a_1}{w}\left(1-\frac{x_2}{K}\right)\frac{w_1x_2}{a_2}$$

**E3**

$$y_3 = \frac{a_1}{w}\left(1-\frac{x_3}{K}\right)\left(\frac{x_3^2}{j}+x_3+a\right),\qquad z_3 = \frac{y_3(w_3-c)-cD_2}{bc}$$

`x3` is a positive root of the cubic

$$f(x_3) = B^{\ast}x_3^3 - \big(B^{\ast}K + (C^{\ast}-B^{\ast})j\big)x_3^2 + j\big((C^{\ast}-B^{\ast})K + aB^{\ast}\big)x_3 + jK(D_2 - B^{\ast}a) = 0$$

with

$$B^{\ast} = \frac{a_1}{w}\left(\frac{w_3}{c} + \frac{a_2w_3b}{w_2c} - 1\right),\qquad C^{\ast} = \frac{a_1w_1w_3b}{w\,w_2\,c}$$

**Existence conditions for E3**

- `x3 < K`
- `w3 > c`
- `y3 > c·D2 / (w3 − c)`
- `D2 + a·a1/w < (a·a1·w3 / (w·c)) · (1 + a2·b / w2)`

---

## 3. Stability analysis

Local stability is determined from the 3×3 Jacobian `J` of the system evaluated at each equilibrium.

### Stability of each equilibrium (as reported in the dissertation)

| Point | Result |
|---|---|
| `E0` | Eigenvalues `0, 0, c`, so non-hyperbolic |
| `E1` | Eigenvalues `0, ±c, ±a1`, so non-hyperbolic |
| `E2` | One eigenvalue is `−c + w3·y2 / (D2 + y2)`; the other two are given in closed form in the dissertation. Asymptotically stable if `c > w3·y2/(D2 + y2)` and the remaining two eigenvalues have negative real part |
| `E3` | Characteristic equation `λ³ + αλ² + βλ + γ = 0`; stable by Routh–Hurwitz (below) |

### Coexistence equilibrium E3

Characteristic polynomial coefficients:

$$\alpha = -(D+G+A),\qquad \beta = DG - EF - BC + AD + AG,\qquad \gamma = A(EF-DG) + BCG$$

**E3 is asymptotically stable when `α > 0`, `β > 0`, `γ > 0` and `αβ > γ`** (Routh–Hurwitz).

<details>
<summary>Jacobian entries A–G at E3</summary>

Let `Φ = x3²/j + x3 + a` and `Ψ = (y3 + D2 + b·z3)²`.

| Entry | Expression |
|---|---|
| `A` | `w·x3·y3·(2x3/j + 1) / Φ² − a1·x3 / K` |
| `B` | `−w·x3 / Φ` |
| `C` | `y3·w1·(a − x3²/j) / Φ²` |
| `D` | `w2·y3·z3 / Ψ` |
| `E` | `−(w2·y3² + D2·w2·y3) / Ψ` |
| `F` | `(D2·w3·z3 + b·w3·z3²) / Ψ` |
| `G` | `−c + (w3·y3² + D2·w3·y3) / Ψ` |

</details>

---

## 4. Numerical simulations

All simulations use the initial condition **(X, Y, Z) = (6, 5, 3)** unless stated otherwise.

### 4.1 Stable case

| a1 | K | w | j | a | a2 | w1 | w2 | D2 | b | c | w3 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 30 | 0.75 | 12 | 60 | 0.05 | 0.6 | 0.2 | 50 | 0.21 | 0.06 | 0.12 |

With these values the system settles to the coexistence equilibrium, i.e. it is asymptotically stable (Figs. 21–22: 3D phase portrait and time series).

### 4.2 Sensitivity to initial conditions

Two nearby initial conditions, **(6, 5, 3)** and **(6.21, 5.13, 3.03)**, start close together but their `X` time series diverge over time (Fig. 23).

- Both trajectories show frequent, abrupt, irregular drops in `X`.
- Both stay bounded (roughly between 0 and 30), so the system does not escape to infinity.
- Together this indicates sensitive dependence on initial conditions, a hallmark of chaos.

### 4.3 Case I: varying `a1`

Fixed parameters:

| K | w | j | a | a2 | w1 | w2 | D2 | b | w3 | c |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 0.75 | 12 | 20 | 0.08 | 0.89 | 0.5 | 50 | 0.21 | 0.12 | 0.01 |

Time series and phase portraits were computed for `a1 = 0.05, 0.1, 0.15, 0.51, 0.85` (Figs. 25–28, 30).

**Observations**

- The bifurcation diagram of the maxima of `x` versus `a1` (Fig. 24) shows a **period-doubling route to chaos** as `a1` increases.
- At higher `a1` the maxima are dense, scattered and irregular, with **windows of periodicity** inside the chaotic region.
- Sensitivity to initial conditions grows with `a1`, then saturates around `a1 ≈ 0.5`, suggesting a relatively stable chaotic regime.
- The **maximum Lyapunov exponent** (Fig. 31) is positive over the whole range, rising from about 0.04 at `a1 = 0.05` to about 0.275 at `a1 = 0.85`. The rise is steepest at small `a1` and flattens beyond roughly 0.3 to 0.5.

### 4.4 Case II: varying `w3`

Fixed parameters:

| a1 | K | w | j | a | a2 | w1 | w2 | D2 | b | c |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 30 | 0.75 | 12 | 20 | 0.08 | 0.89 | 0.5 | 50 | 0.20 | 0.12 |

Time series and phase portraits were computed for `w3 = 1.2, 0.75, 0.5, 0.25, 0.15, 0.1` (Figs. 32–37).

**Observations**

- Varying `w3` produces **chaotic windows interleaved with more regular behavior**, i.e. several transitions in the intensity of chaos.
- For large `w3` (Figs. 32–34) the dynamics are burst-like, with sparse large peaks in the top predator `Z` while `X` stays near its carrying capacity.
- For small `w3` (Figs. 35–37) the oscillations are denser and more irregular.
- The **maximum Lyapunov exponent** (Fig. 38) stays positive, between about 0.253 and 0.28, and is non-monotonic in `w3`, with local peaks near `w3 ≈ 0.75` and `w3 ≈ 1.1` and a dip near `w3 ≈ 0.9`.

---

## 5. Applications and implications

Mapping the model to a forest ecosystem:

| Variable | Example |
|---|---|
| `x(t)` | Prey (e.g. grass) |
| `y(t)` | Predator (e.g. deer) |
| `z(t)` | Top predator (e.g. tiger) |

Potential uses:

- Maintaining biodiversity and preventing degradation of forest ecosystems.
- Managing human–wildlife conflict.
- Predicting the ecological outcome of predator removal or overhunting (deer overpopulation, vegetation depletion).
- Simulating the long-term impact of ecological policies before they are implemented.

**Real-life illustration: Yellowstone National Park.** Deer and elk overpopulation led to overgrazing and damaged vegetation. The reintroduction of wolves in 1995 reduced elk numbers and let vegetation recover, benefiting species such as beavers and birds. The dissertation uses this as an example of sensitive dependence on initial conditions, where small ecological changes can have large, unpredictable effects.

---

## 6. Key findings

- The three-species model has a rich dynamical structure: stable, periodic and chaotic regimes all occur.
- With the Case I parameters the model is **chaotic for `a1` from 0.05 to 0.85**, and with the Case II parameters **chaotic for `w3` from 0.1 to 1.2**.
- Chaos is reached through **period-doubling** as `a1` increases, confirmed by the bifurcation diagram, phase portraits and positive Lyapunov exponents.
- Trajectories remain bounded while being sensitive to initial conditions.
