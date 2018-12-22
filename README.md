  # PolynomialPDEs
\maketitle
An implementation of DTM for finding the solution of Boundary Value Problems in Calculus of Variations

We cite from \lq\lq An analytic study on the Euler-Lagrange equation arising in calculus of variations":

Euler-lagrange equation plays an important role in the minimizing problems of the calculus of variations. Differential Transform Method (DTM), provides an anaylitical solution in the form of an infinite power series with with easily computable components. DTM is a semi-numerical-analytic technique that formalizes taylor series. This method gives exact values of the n-th of an analytic function at a point in terms of known and unknow boundary conditions in a fast manner.

Let $x(t)$ be analytic in a domain $D$ and let $t_i$ represent any point in $D$. Along this work, we will consider and ask for $t_0=0\in D$. The function $x(t)$ is then represented by one power series whose center is located at $t_i$. The $k$th differential transformation of the function $x(t)$ is defined as follows:


\begin{equation}
X(k)=\frac1{k!}\cdot\frac{d^kx(t)}{dt^k}\Big\rvert_{t=t_i}.
\end{equation}

In this equation, $x(t)$ is the original function and $X(k)$ is the fransformed function. The inverse differential transformation is given by:

\begin{equation}
x(t)=\underset{k=0}{\overset{\infty}\sum}(t-t_i)^k X(k),\ \forall t\in D.
\end{equation}

In our case, since we are considering $t_0=0$ this simply becomes,

\begin{equation}
x(t)=\underset{k=0}{\overset{\infty}\sum}t^k X(k),\ \forall t\in D.
\end{equation}

where we evaluate the derivatives of $x$ in $t_0=0$ when computing $X(k)$.

\noindent Some of the fundamental mathematical operations (calculated when $t_0=0$) performed by differential
transform method are listed in Table 1.

\begin{center}
\begin{tabular}{ l r }
\hline
  Original function & Transformed function\\
\hline
 $x(t) = \alpha u(t) \pm \beta w(t)$  &  $X(k) = \alpha U(k) \pm \beta W (k)$\\
 $x(t) = u(t)v(t)$   &  $X(k)=\Sigma_{l=0}^k U(l)V(k - l)$ \\
 $x(t)=\frac{d^m}{dt^m}u(t)$   &   $X(k) = \frac{(k+m)!}{k!} U(k+m)$\\
 $x(t) = t^m$   &   $X(k) = \delta_m(k - m)$\\
 $x(t) = e^{mt}$   &   $X(k) = \frac{m^k}{k!}$\\
 \hline
\end{tabular}
\end{center}

\newpage
\section{Description of the problem}
Consider the problem of finding the extremum of the functional

\begin{equation}
J[y(t)]=\int_{t_0}^c L(t,y(t),y'(t))dt.
\end{equation}

The necessary condition for $y(t)$ to \lq\lq extremize" $J[y(t)]$ is that it should satisfy the
Euler-Lagrange equation

\begin{equation}\label{ele}
\frac{\partial L}{\partial y} - \frac{d}{dt}\Biggl(\frac{\partial L}{\partial y'}\Biggr) = 0,
\end{equation}

with appropriate boundary conditions. The boundary value problem does not
always have a solution and if the solution exists, it may not be unique. It is worth to
mention here that, if the solution of Euler-Lagrange satisfies the boundary conditions,
it is unique. The given boundary conditions will be expressed as:

\begin{equation}
y(0) = a
\end{equation}
\begin{equation}\label{bv1}
 y(1) = b 
\end{equation}
\begin{equation}
y'(0) = \gamma
\end{equation}

where for simplicity, we consider $c=1$ and $\gamma\in\mathbb R$ is an unknown parameter. Notice that $y(0)=Y(0)$ and $Y(1)=\gamma$.

\section{Solution with DTM}
We know that $y(t)$ can be obtained as

\begin{equation}
y(t)=\underset{k=0}{\overset{\infty}\sum}t^k Y(k),\ \forall t\in D.
\end{equation}

\noindent In real applications, the function y(t) can be approximated by a finite series as

\begin{equation}
y(t)=\underset{k=0}{\overset{N}\sum}t^k Y(k),\ \forall t\in D.
\end{equation}

Now we are only missing the factors $Y(k)$, this is fixed using equation \ref{ele}. Euler-Lagrange equation is in the form $E(t,y',y'',\dots) = 0$. When the highest order derivative involved in $E$ is $y''$ then we can use table 1 along $E$ to find $Y(k)$ in a recursive manner, depending on $\gamma$, for $k\geq2$ (we'll denote this as $Y(k,\gamma)$).

Finally, we can use $\ref{bv1}$ to aproximate $\gamma$. It's

\begin{equation}
b=y(1)=y(0)+\gamma+\sum_{k=2}^NY(k,\gamma).
\end{equation}

\noindent Solving for $\gamma$, we obtain the expected model.

\newpage
\section{An example}
