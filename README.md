# twoBinaryTrial (R package)

An R package for **probability calculations in clinical trials with two binary endpoints**
(e.g., **efficacy** and **safety**). The package builds the **joint distribution** of the
two endpoints using a **2×2 multinomial model** (four outcome cells) and an association
parameter defined by an **odds ratio** \(\phi\).

This package was created to support a research workflow where functions for:
- converting marginal probabilities \((p_E, p_S)\) and \(\phi\) into the four cell probabilities, and
- repeatedly computing joint/multinomial probabilities (e.g., for PCS calculations)

need to be reused consistently and reproducibly across simulations and manuscript results.

> **Note:** Replace `twoBinaryTrial` with your actual package name if different.

---

## Problem setup

Let the two binary endpoints be:

- **Efficacy** \(E \in \{0,1\}\)
- **Safety** \(S \in \{0,1\}\)

A subject falls into one of four joint outcomes:

| Cell | Meaning | Probability |
|---|---|---|
| \(p_{11}\) | \(E=1, S=1\) | \(P(E=1, S=1)\) |
| \(p_{10}\) | \(E=1, S=0\) | \(P(E=1, S=0)\) |
| \(p_{01}\) | \(E=0, S=1\) | \(P(E=0, S=1)\) |
| \(p_{00}\) | \(E=0, S=0\) | \(P(E=0, S=0)\) |

Given:
- marginal efficacy probability \(p_E = P(E=1)\),
- marginal safety probability \(p_S = P(S=1)\),
- an odds ratio (association) \(\phi\),

we solve for \((p_{11}, p_{10}, p_{01}, p_{00})\) such that:
- \(p_{10} = p_E - p_{11}\)
- \(p_{01} = p_S - p_{11}\)
- \(p_{00} = 1 - p_E - p_S + p_{11}\)
- \(\text{OR} = \dfrac{p_{11}p_{00}}{p_{10}p_{01}} = \phi\)

This yields a (typically) quadratic equation in \(p_{11}\). The valid solution is the root
that lies within:
\[
\max(0, p_E + p_S - 1) \le p_{11} \le \min(p_E, p_S)
\]
When \(\phi = 1\), the endpoints are independent and \(p_{11} = p_E p_S\).

Once the cell probabilities are obtained, the subject-level joint outcomes follow a
**Multinomial** distribution with 4 categories (a 2×2 table). For a sample size \(n\),
the count vector \((N_{11}, N_{10}, N_{01}, N_{00})\) satisfies:
\[
(N_{11}, N_{10}, N_{01}, N_{00}) \sim \text{Multinomial}\left(n; p_{11}, p_{10}, p_{01}, p_{00}\right).
\]

---

## Installation

### From source (local)
```r
# install.packages("devtools")
devtools::install_local("path/to/twoBinaryTrial")
