Strange fitting problem
-----------------------

The file `parampoints.txt` contains a point set with parameters (with the format `u v x y z`),
where the parameters are in [0, 1]. There is a larger dataset in `parampoints2.txt`.

The starting B-spline surface is contained in `surface.bss`. The first control row,
associated with `u=0`, is retained, and the remaining control points are used to approximate
the given points. Even though the input is approximately symmetric, the output is not.

![image](problem.png?raw=true "The fitted B-spline is not symmetric")

(The program generates two files: `points.obj`, which is just the input points, and `result.obj`,
 that contains the control points of the fitted B-spline surface.)
 
Experiments:

1. Make the data really symmetric (`type=1` in `readParamPoints`):
   The original data is approximately symmetric to the `y` axis;
   here this is enforced by mirroring everything with a negative `y` coordinate,
   and deleting the points with positive `y` coordinates.
   The result is substantially different (!), but nowhere near symmetric.
2. Use an algebraic function (`type=2` in `readParamPoints`):
   The input points are replaced by a parabolic function of the parameters, 
   that runs close to the original data.
   The result is nice and symmetric.
3. Add more control points (using `surface2.bss`).
   There are 10 control points to use in the fitting now, and this also works nicely.
4. Use a different library (`test2` instead of `test`):
   Implement the same algorithm with [libgeom](http://github.com/salvipeter/libgeom)
   instead of the original *Sketches* library. No change.

So what is the problem?

- It cannot be the input surface, because of experiment #2.
- It cannot be the input points, because of experiment #3.
- It cannot be the library, because of experiment #4.
