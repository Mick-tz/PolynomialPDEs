# PolynomialPDEs
An implementation of DTM for finding the solution of Boundary Value Problems in Calculus of Variations

We cite from "An analytic study on the Euler-Lagrange equation arising in calculus of variations"

Euler-lagrange equation plays an important role in the minimizing problems of the calculus of variations. Differential Transform Method (DTM), provides an anaylitical solution in the form of an infinite power series with with easily computable components. DTM is a semi-numerical-analytic technique that formalizes taylor series. This method gives exact values of the n-th of an analytic function at a point in terms of known and unknow boundary conditions in a fast manner.

Let $x(t)$ be analytic in a domain $D$ and let $t_i$ represent any point in $D$. Along this work, we will consider and ask for $t_0=0\in D$. The function $x(t)$ is then represented by one power series whose center is located at $t_i$. The $k$th differential transformation of the function $x(t)$ is defined as follows:

![f1]

\begin{equation}

\end{equation}

In this equation, $x(t)$ is the original function and $X(k)$ is the fransformed function. The inverse differential transformation is given by:

\begin{equation}
x(t)=\underset{k=0}{\overset{\infty}\sum}(t-t_i)^k X(k),\ \forall t\in D.
\end{equation}

In our case, since we are considering $t_0=0$ this simply becomes,

\begin{equation}
x(t)=\underset{k=0}{\overset{\infty}\sum}(t)^k X(k),\ \forall t\in D.
\end{equation}

where we evaluate the derivatives of $x$ in $t_0=0$ when computing $X(k)$.

Some of the fundamental mathematical operations performed by differential
transform method are listed in Table 1.

\begin{tabular}{ l r }
\hline
  Original function & Transformed function
\hline
 $x(t) = \alpha u(t) \pm \beta w(t)$  &  $X(k) = \alpha U(k) \pm \beta W (k)$\\
 $x(t) = u(t)v(t)$   &  $X(k) =\sum_{l=0}^k U(l)V (k − l)$ \\
 $x(t)=\frac{d^m}{dt^m}u(t)$   &   $X(k) = $\frac{(k+m)!}{k!} U(k+m)$\\
 $x(t) = t^m$   &   $X(k) = \delta_m(k − m)$\\
 $x(t) = e^{mt}$   &   $X(k) = \frac{m^k}{k!}$\\
 \hline
\end{tabular}

## Description of the problem



[f1]:http://chart.apis.google.com/chart?cht=tx&chl=$X(k)=\frac1{k!}\cdot\frac{d^kx(t)}{dt^k}\rvert_{t=t_i}.$
