| EAPS 591 - Numerical Dynamics 	| Fall 2020 	| Prof. David Minton 	|
|:-----------------------------:	|:---------:	|:------------------:	|
# Module 3

## Prerequisites

I recommend creating a new VS Code workspace and storing it in the new `module3-[USERNAME]` directory that you cloned from using git.

Once you have created a new workspace, it is important to make sure that your Python environment is set correctly. Go to the *Command Pallette* from the *Settings* icon (or, alternatively use the keyboard shotcut <kbd>⇧</kbd><kbd>⌘</kbd><kbd>P</kbd>) and search for *Python: Select Interpreter*. Select the current workspace and then from the drop-down, select `('EAPS591-Dynamics-F2020': conda)`. 

Be sure to commit often, and be sure to push your commits up to GitHub at a miniumum prior to when we meet for class so that I can track your progress.


--- 
## Problem 1: Data processing
The `data` folder contains several files named `idXXXXXX-XV.csv`. Each file contains the cartesian position and velocity vectors for a single body in a 1 My n-body simulation. The units in these data files are AU-day. 

1. In Jupyter Notebook, read the data files into a Pandas dataframe.

2. Make a plot that contains the cartesian positions for all times in the X-Y plane. There must be a single plot with each separate body in a different color.

3. Repeat 2 but this time plot the X-Z plane.

## Problem 2: Coordinate conversions.
Write a Modern Fortran code that converts cartesian coordinate vectors ($\mathbf{r}, \mathbf{v}$) to orbital elements ($a$, $e$, $I$, $\Omega$, $\varpi$, $f$). 

1. Use the code to convert all of the data files from Problem 1 to orbital elements. Save the output files as `idXXXXXX-EL.csv`.

2. Make plots of $a$ vs $e$ and $a$ vs $i$ for each body (see Figs. 7.9-7.11 on pages 305-307 in Murray & Dermott as example of what they should look like).

3. On a single plot, plot $a$, $r_q$, and $r_Q$ vs. $t$. You may wish to use a semilog scale (i.e. $\log r$ vs $t$). In a Jupyter Notebook, plot it such that each body has its own color (using the same color scheme as in Problem 1), but $a$ is a solid bold line and the region between $r_q$ and $r_Q$ is shaded and partially transparent.

## Problem 3: Coordinate conversion the other way.
Write a Modern Fortran code that converts orbital elements ($a$, $e$, $I$, $\Omega$, $\varpi$, $f$) to cartesian coordinate vectors ($\mathbf{r}, \mathbf{v}$). 

1. Use the code to convert the data files you generated in Problem 2 *back* to cartesion position and velocity vectors. Save the output files as `idXXXXXX-XV-new.csv`

2. Verify that you have done the conversion correctly by loading the files into a Jupyter Notebook along with the original files. To do this, make two plots, one with $\left|\mathbf{r}_{orig}-\mathbf{r}_{new}\right|$ and $\left|\mathbf{v}_{orig}-\mathbf{v}_{new}\right|$ vs time for particle number 10 only.

## Problem 4: The Leapfrog method method for the 2 body problem.
Consider the dynamics of the two body problem in three dimensions. The two-body equation of motion can be written as:
$$\ddot{\mathbf{r}}=-\mu\frac{\mathbf{r}}{r^3}$$
Consider the \motion of a particle where:
$$\mathbf{r} = r_x\hat{x}+r_y\hat{y}+r_z\hat{z}$$
The equation of motion in vector form can now be written as three coupled scalar equations:
$$
\begin{aligned}
\ddot{r}_x&=-\mu\frac{r_x}{\left(r_x^2+r_y^2+r_z^2\right)^{3/2}}\\
\ddot{r}_y&=-\mu\frac{r_y}{\left(r_x^2+r_y^2+r_z^2\right)^{3/2}}.\\
\ddot{r}_z&=-\mu\frac{r_z}{\left(r_x^2+r_y^2+r_z^2\right)^{3/2}}.\\
\end{aligned}
$$
Ultimately, we'd like to use these equations solve problems that can be written as an *initial value problem* in which the position and velocity of particles is given at time $t_0$ and we wish to know what the position and velocity will be at some time in the future.  

For each of the dimensions ($x$, $y$, and $z$), the equation of motion is a 2$^{nd}$ order ODE. We can write it instead as a system of 1$^{st}$ order ODEs:
$$
\begin{aligned}
\dot{r}_x&=v_x\\
\dot{v}_x&=-\mu\frac{r_x}{\left(r_x^2+r_y^2+r_z^2\right)^{3/2}}.\\
\end{aligned}
$$

We can write our system of 1$^{st}$ order ODEs as vectors in the form:
$$
\mathbf{x}=\begin{bmatrix}
r_x\\
r_y\\
r_z\\
\end{bmatrix}
$$
$$
\mathbf{y}=\begin{bmatrix}
v_x\\
v_y\\
v_z\\
\end{bmatrix}
$$
with the following derivatives:
$$\dot{\mathbf{x}}=\mathbf{y}\quad\mathrm{and}\quad\dot{\mathbf{y}}=f\left(\mathbf{x}\right)$$
where
$$f\left(\mathbf{x}\right)=\begin{bmatrix}
-\mu\frac{r_x}{\left(r_x^2+r_y^2+r_z^2\right)^{3/2}}.\\
-\mu\frac{r_y}{\left(r_x^2+r_y^2+r_z^2\right)^{3/2}}.\\
-\mu\frac{r_z}{\left(r_x^2+r_y^2+r_z^2\right)^{3/2}}.\\
\end{bmatrix}
$$

We can choose to work in a unit system in which $\mu=1$. In this unit system, consider an elliptical orbit with semimajor axis $a=1$ and $e=0.3$, and choose initial conditions $\left(r_{x,0}, r_{y,0}, v_{x,0}, v_{y,0}\right)$ such that the body is at the periapsis of the orbit.

Given initial values for $\mathbf{x}_0$, and $\mathbf{y}_0$, as well as a step-size $h$ (aka time step, or $\Delta t$), the above system of equations can be solved numerically to find $\mathbf{x}_n$, and $\mathbf{y}_n$ where $t_n = t_0 + nh$. 

For this problem, you will implement the *Leapfrog method*. In the "drift-kick-drift" form, the Leapfrog method can be written as:

$$
\begin{aligned}
\mathbf{x}_{n+1/2} &= \mathbf{x}_n + \frac{1}{2}h \mathbf{y}_n\\
\mathbf{y}_{n+1} &= \mathbf{y}_n + h f\left(\mathbf{x}_{n+1/2}\right)\\
\mathbf{x}_{n+1} &= \mathbf{x}_n + \frac{1}{2}h \mathbf{y}_{n+1}\\
\end{aligned}
$$

Apply the Leapfrog method to integrate the nine bodies given in the data files analyzed previously. Assume that these bodies are "massless test particles." and be sure to use the correct units for $\mu$. Integrate the bodies for 1 My forward in time using a step size that is constrained by the innermost body in the system. You may need to experiment to determine a good step size. 

As your integration proceeds, monitor the total energy $\mathcal{E}=\frac{1}{2}v^2-\frac{1}{r}$ and angular momentum $\mathbf{h}=\mathbf{r}\times\mathbf{v}$ and plot $\Delta\mathcal{E}/\mathcal{E}_0$ and $\Delta h/h_0$ for each of the bodies in the system. 


## Problem 10: The Symplectic Integrator
Write a symplectic n-body integrator. You may use either Jacobi coordinates or democractic heliocentric coordinates, but be sure to use the correct form of the Hamiltonian depending on which method you use. Only implement the basic symplectic integrator, so do not worry about close encounters. Simulate the orbit of Pluto and Neptune for $10^5$ years (bodies 9 and 10). Produce plots of the total change in the system energy, $\Delta\mathcal{E}/\mathcal{E}_0$ and the resonance angle $\phi=3\lambda_P-2\lambda_N-\varpi_P$ vs. time. 








