# Lorenz_periodic_orbit
These are the key codes related to the periodic orbits of the Lorenz system, including the generation of the Lorenz attractor, similarity signature curves, Poincaré sections, quasi-periodic orbits, IntLab toolbox, and the Krawczyk operator. The parameters in the program can be adjusted to obtain different divisions of the quasi-periodic orbit.

​		Lorenz.m: Generate the Lorenz attractor. We can obtain a trajectory of shape 3 × iter, where each column represents a point on the trajectory. “iter” is the number of iterations.



​		compsig.m: Calculate the signature curves and similarity signature curves.  We can obtain a similarity signature curve of shape 3 × iter, where each column represents the similarity signature at a point on the trajectory of the Lorenz attractor.



​		quasi.m: Generate the quasi-periodic orbits. The parameters in the program can be adjusted again. Different effects can be observed by adjusting parameters such as initial point, window size, etc. We need to verify the points obtained in the quasi-periodic orbit one by one through the Krawczyk operator.



​		Poincare.m: Generate the Poincare section. The trajectory of the points on the Poincare section can be obtained step by step by changing the evolution time.



​		Intlab.zip: Interval analysis. The basic usage can be found in the literature. See: Hargreaves G I. Interval analysis in MATLAB[J]. Numerical Algorithms, 2002 (2009.1). IntLab V11 and higher versions can be more conducive to the calculation of interval vectors. In addition, the enclosure will be easily obtained by the Taylor method or the AWA method



​		Krawczyk.zip: Construct Krawczyk operator. We need to use the IntLab toolbox to generate interval vectors and then check the points one by one using the Krawczyk operator. All computations are carried out in interval arithmetic, which ensures that we obtain enclosure of the solution of the system. 

​				fun lorenz_rhs: Lorenz attractor function handles

​				fun lorenz_jacobian: Jacobian matrix of the Lorenz System

​				fun R: Poincare map

​				fun R_der: Calculating the Jacobi matrix for Poincare map

​				fun F: the existence of period-p orbits of R

​				fun Krawczyk:  Constructing the Krawczyk operator

​				fun rk4: Runge-Kutta methods



## Remark

When performing interval calculations, it is important to initialize the program variables with intervals. This can be achieved by creating intervals using the interval data type in the initialization section of the program. For instance, for a variable named x, an interval can be initialized as:

```matlab
x = intval([a, b]) 
```

or

```matlab
x = intval(zeros(a,b));
```

Once all the variables have been initialized with intervals, interval operators such as +, -, *, and / can be used to perform interval calculations. These operators will use the range of the intervals to compute the result and return a new interval containing the resulting range.

It is important to keep in mind that numerical errors in interval calculations can cause the resulting interval to be wider, requiring a more conservative analysis of results. 

