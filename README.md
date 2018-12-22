# PolynomialPDEs
An implementation of DTM for finding the solution of Boundary Value Problems in Calculus of Variations

We cite from "An analytic study on the Euler-Lagrange equation arising in calculus of variations"

Euler-lagrange equation plays an important role in the minimizing problems of the calculus of variations. Differential Transform Method (DTM), provides an anaylitical solution in the form of an infinite power series with with easily computable components. DTM is a semi-numerical-analytic technique that formalizes taylor series. This method gives exact values of the n-th of an analytic function at a point in terms of known and unknow boundary conditions in a fast manner.

Let $x(t)$ be analytic in a domain $D$ and let $t_i$ represent any point in $D$. Along this work, we will consider and ask for $t_0=0\in D$. The function $x(t)$ is then represented by one power series whose center is located at $t_i$. The $k$th differential transformation of the function $x(t)$ is defined as follows:
\begin{equation}
 X(k)=\frac1{k!}\cdot \frac{d^k x(t)}{dt^k}\big|_{t=t_i}\ ,\ \forall t\in D.
\end{equation}



