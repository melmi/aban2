# Single phase solver

## Derivative boundary conditions

In 1d single phaze case, derivatives at the starting boundary and starting cells can be written as:

$$ 
\begin{align} 
\frac{\partial \phi}{\partial x} \vert_b & = \frac{1}{\Delta x}2(\phi_0-\phi_b) = \frac{1}{2\Delta x}4(\phi_0-\phi_b) \\\\
\frac{\partial \phi}{\partial x} \vert_0 & = \frac{1}{2\Delta x}(\phi_0+\phi_1-2\phi_b) \\\\
\frac{\partial \phi}{\partial x} \vert_1 & = \frac{1}{2\Delta x}(\phi_2-\phi_0) \\\\
\frac{\partial \phi}{\partial x} \vert_2 & = \frac{1}{2\Delta x}(\phi_3-\phi_1) 
\end{align} 
$$

## Compact Laplace's operator boundary conditions

$$
\begin{align} 
\frac{\partial^2 \phi}{\partial x^2} \vert_0 
& = \frac{1}{\Delta x^2}[(\phi_1-\phi_0)-2(\phi_0-\phi_b)] \\\\
& = \frac{1}{\Delta x^2}(-3\phi_0+\phi_1+2\phi_b)
\end{align}
$$

## Extended Laplace's operator boundary conditions

$$
\begin{align} 
\frac{\partial^2 \phi}{\partial x^2} \vert_0 
& = \frac{1}{4\Delta x^2}[(\phi_0+\phi_1-2\phi_b+\phi_2-\phi_0)-8(\phi_0-\phi_b)] \\\\
& = \frac{1}{4\Delta x^2}(-8\phi_0+\phi_1+\phi_2+6\phi_b)
\end{align}
$$

$$
\begin{align} 
\frac{\partial^2 \phi}{\partial x^2} \vert_1
& = \frac{1}{4\Delta x^2}[(\phi_2-\phi_0+\phi_3-\phi_1)-(\phi_0+\phi_1-2\phi_b+\phi_2-\phi_0)] \\\\
& = \frac{1}{4\Delta x^2}(-\phi_0-2\phi_1+\phi_3+2\phi_b)
\end{align}
$$
