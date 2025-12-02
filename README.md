# SSLfmm

**Semisupervised Learning under a Mixed-Missingness Mechanism in Finite Mixture Models**

`SSLfmm` is an R package for semi-supervised learning in finite mixture models
under a *mixed label-missingness mechanism*. The package explicitly models both
missing completely at random (MCAR) and entropy-based missing at random (MAR)
processes through a logistic–entropy formulation. Parameter estimation is
conducted via an Expectation–Conditional Maximisation (ECM) algorithm with
robust initialisation for stable convergence.

The methodology is grounded in the statistical theory of informative
missingness and partially labelled mixture modelling as discussed in:

- Ahfock and McLachlan (2020), *Statistics and Computing*  
- Ahfock and McLachlan (2023), *Econometrics and Statistics*

---

## 📦 Installation

You can install the development version of `SSLfmm` directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("wujrtudou/SSLfmm")
