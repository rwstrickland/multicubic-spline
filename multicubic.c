#include <stdio.h>
#include <stdlib.h>

#define NATSPLN 1.e30

int
csp1d(double x[], double y[], int n, double ypl, double ypr, double y2[])
/*
 * Purpose
 *
 *   Caclulate cubic spline of a tabulated smooth function of one variable
 *
 * Reference
 *
 *   Adapted from http://oceancolor.gsfc.nasa.gov/staff/norman/
 *      seawifs_image_cookbook/faux_shuttle/spline.c
 *      by Norman Kuring 31 Mar 1999:
 *
 * Comments from Kuring
 *
 *   This code is based on the cubic spline interpolation code presented in
 *   Numerical Recipes in C: The Art of Scientific Computing by
 *   William H. Press, Brian P. Flannery, Saul A. Teukolsky, and
 *   William T. Vetterling, 1988 (and 1992 for the 2nd edition).
 *
 *   This version assumes zero-offset arrays instead of the unit-offset arrays
 *   suggested by the authors.  You may style me rebel or conformist
 *   depending on your point of view.
 *
 * Revised
 *
 *   Robert Strickland September 21, 2015: Added comments and return values
 *
 * Input
 *
 *   double x[]  independent variable x coordinates, strictly increasing order
 *   double y[]  value of function at each x
 *   int n       Length of x and y
 *   ypl         First derivative at left end, x[0], or >= 1.e30 to signal
 *               to use natural spline (2nd derivative zero)
 *   ypr         First derivative at right end, x[n-1], or >= 1.e30 to
 *               signal to use natural spline
 *
 * Output
 *
 *   y2[]        The second derivative of spline at the nodes. See below.
 *               This changes linearly between the nodes.
 *
 * Returns
 *
 *   0  no error
 *   1  memory allocation error
 *   2  x[] not in strictly increasing order
 *
 * Discussion
 *
 *   Calculate cubic spline versus x through n points (x[i], y[i]). The
 *   independent variable must be in strictly increasing order,
 *   x[i+1] > x[i].  The first and second derivatives of the resulting
 *   spline will be continous. To evaluate the spline at xx, find the
 *   interval i such that x[i] <= xx < x[i+1]. Then let
 *
 *     dx = x[i+1] - x[i]
 *     a = (x[i+1] - xx)/h
 *     b = 1. - a
 *
 *   The value of the cubic spline s(xx) at the point xx is
 *
 *     f(xx) = a*y[i] + b*y[i+1] +
 *            ((a*a*a - a)*y2[i] + (b*b*b - b)*y2[i+1])*(dx*dx)/6.
 */
{
  int i, k;
  double p, qn, sig, un, *u, dx;

  u = (double *) malloc((n-1)*sizeof(double));
  if (u == NULL)
    return 1;

  if (ypl > 0.99e30) { /* Use natural spline at left end */
    y2[0] = 0.;
    u[0] = 0.;
  }
  else { /* User provided 1st derivative at left end */
    y2[0] = -0.5;
    u[0] = (3./(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-ypl);
  }
  for (i = 1; i < n-1; i++) {
    dx = x[i] - x[i-1];
    if (dx <= 0.) {
      free(u);
      return 2;
    }
    sig = dx/(x[i+1] - x[i-1]);
    p = sig*y2[i-1] + 2.;
    y2[i] = (sig - 1.)/p;
    u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/dx;
    u[i] = (6.*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
  }
  if (ypr > 0.99e30) { /* Use natural spline at right end */
    qn = 0.;
    un = 0.;
  }
  else { /* User provided 1st derivative at right end */
    qn = 0.5;
    un = (3./(x[n-1] - x[n-2]))*(ypr - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
  }
  y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.);
  for (k = n-2; k >= 0; k--)
    y2[k] = y2[k]*y2[k+1] + u[k];
  free(u);
  return 0;
}

void
csp2d_close(double **splines)
/*
 * Free memory allocated by csp2d_open()
 */
{
  int i;

  if (splines == NULL)
    return;
  for (i = 0; i < 3; i++) {
    if (splines[i] != NULL) {
      free(splines[i]);
      splines[i] = NULL;
    }
  }
  free(splines);
  splines = NULL;
}

double **
csp2d_open(double *xv[], int nv[], double *f)
/*
 * Purpose
 *
 *   Calculate 2-D cubic spline coefficients for a smooth function f(x, y)
 *   of two independent variables tabulated on a 2-D rectilinear grid.
 *   To evaluate this 2-D spline, use csp2d_eval(), one call per point.
 *
 * Input
 *
 *   double *xv[]  Array of 2 arrays of independent variables;
 *                 xv[] is the array of {x[], y[]} Both x[] and y[]
 *                 must be in strictly-increasing order
 *   int nv[]      Array of length 2 containing lengths of xv[];
 *                 nv[] = {nx, ny}
 *   double f[]    f(x[i], y[j]) packed into a one-dimensional
 *                 array of length nx*ny. The y[j] index varies the
 *                 fastest. That is, instead of f[i][j], we use f[i*ny+j]
 * Returns
 *
 *   double **splines  On success, returns spline coefficient array. This is
 *                     an array of length 3 of arrays of length nx*ny. The
 *                     user passes this to csp2d_eval(). Call
 *                     sp2d_close(splines) to free this memory.
 *
 *                     On error, returns NULL. Errors are memory
 *                     allocation or xv[][] not in strictly increasing order.
 *
 * Method
 *
 *   This calls csp1d() to calculate a 1d spline along lines through
 *   every grid point parallel to each axis in turn. The csp1d() calculates
 *   the second derivative versus the independent variable. These are
 *   stored in an array the same size as f[]. Then cascaded cubic splines of
 *   these spline coefficients are calculated and stored. See csp2d_eval() to
 *   see how these are used to evaluate the spline.
 *
 * By
 *
 *   Robert W. Strickland, September 21, 2015
 */
{
  int nw;
  double *w, *s;
  int i, j;
  double **splines;
  double *sx, *sy, *sxy;
  int nf;
  double *x, *y;
  int nx, ny;
  int nerr;

  x = xv[0];
  y = xv[1];
  nx = nv[0];
  ny = nv[1];

  nw = (nx > ny) ? nx : ny;
  s = (double *) malloc(nw * sizeof(double));
  w = (double *) malloc(nw * sizeof(double));
  if (s == NULL || w == NULL)
    return NULL;

  nf = nx*ny;

  splines = (double **) malloc(3*sizeof(double *));
  if (splines == NULL)
    return NULL;
  for (i = 0; i < 3; i++)
    splines[i] = NULL;
  for (i = 0; i < 3; i++) {
    splines[i] = (double *) malloc(nf*sizeof(double));
    if (splines[i] == NULL) {
      csp2d_close(splines);
      return NULL;
    }
  }
    
  /* This must match the order in csp2d_eval(); */
  sx = splines[0];
  sy = splines[1];
  sxy = splines[2];

  nerr = 0;

  /*
   * Spline of f vs x -> sx
   */
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++)
      w[i] = f[i+nx*j];
    nerr += csp1d(x, w, nx, NATSPLN, NATSPLN, s);
    for (i = 0; i < nx; i++)
      sx[i+nx*j] = s[i];
  }
  /* Splines of f vs y -> sy
   * Cascaded spline of sx vs y -> sxy
   */
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++)
      w[j] = f[i+nx*j];
    nerr += csp1d(y, w, ny, NATSPLN, NATSPLN, s);
    for (j = 0; j < ny; j++)
      sy[i+nx*j] = s[j];
    for (j = 0; j < ny; j++)
      w[j] = sx[i+nx*j];
    nerr += csp1d(y, w, ny, NATSPLN, NATSPLN, s);
    for (j = 0; j < ny; j++)
      sxy[i+nx*j] = s[j];
  }
  free(s);
  free(w);
  if (nerr != 0) {
    csp2d_close(splines);
    splines = NULL;
  }
  return splines;
}

double
csp2d_eval(double p[], double *xv[], int nv[], double f[],
    double **splines, int interval[])
/*
 * Purpose
 *
 *   Interpolate a smooth function of two variables tabulated on a rectilinear
 *   grid. This uses the coefficients calculated by csp2d_open().
 *
 * Input
 *
 *   double p[2]  The point P(x, y) where we want f interpolated
 *   double xv[2] The array of {x[], y[]} where f is tabulated.
 *   int nv[2]    The array of {nx, ny}, the lengths of the 
 *                {x[], y[]} arrays
 *   double f[]   f(x[i], y[j]) packed into a one-dimensional
 *                array of length nx*ny. The y[j] index varies the
 *                fastest. That is, instead of f[i][j], we use f[i*ny+j]
 *   double **splines
 *                Spline coefficients returned by csp2d_open()
 *   int interval[2]
 *                Work array of length 2 to keep track of interval. If
 *                there are several parallel calls at the same value of p[],
 *                they should all share the same interval[] array.
 *
 * Method
 *   
 *   We start by finding the rectangle in our grid that contains P. Find
 *   i and j such that
 *
 *     x[i] < xp < x[i+1]
 *     y[j] < yp < y[j+1]
 *
 *   If the point falls off the edges of the table, the function will
 *   silently extrapolate. Results will probably be poor.
 *
 *   Interpolate along the 2 edges of the rectangle parallel to the y axis
 *   to get a line segment, eliminating the y dependence. Finally, we
 *   interpolate along the line segment to get a single point and to eliminate
 *   the x dependence. This can be considered the process of repeatedly
 *   slicing perpendicular to one of the axes, reducing the dimension by one
 *   each time until a single point remains.
 *
 *   To get a better understanding, of the underlying process it might be
 *   best to read this code from the bottom up. We'll use the following
 *   scheme to index the points:
 *
 *      y[j+1] -01-----+---11
 *               |     |    |
 *          yy --U-----P----V
 *               |     |    |
 *               |     |    |
 *        y[j] -00-----+---10
 *               |     |    |
 *             x[i]   xx  x[i+1]
 *   
 * By
 *
 *   Robert W. Strickland, September 21, 2015
 */
{
  double a, b, g, c, d;
  int i0, i1, j0, j1;
  int k00, k01, k10, k11;
  double dx, dy;
  double f0, f1, s0, s1;
  int n;
  double *sx, *sy, *sxy;
  double *x, *y;
  int nx;
  double xx, yy;

  x = xv[0];
  y = xv[1];
  nx = nv[0];
  xx = p[0];
  yy = p[1];

  /* This must match the order in csp2d_open() */

  sx = splines[0];
  sy = splines[1];
  sxy = splines[2];

  /* Linear search for interval containing each independent variable */

  for (n = 0; n < 2; n++) {
    if (interval[n] < 0 || interval[n] > nv[n]-2)
      interval[n] = nv[n]/2;
    while (interval[n] > 0 && p[n] < xv[n][interval[n]])
      interval[n]--;
    while (interval[n] < nv[n]-2 && p[n] > xv[n][interval[n]+1])
      interval[n]++;
  }

  i0 = interval[0];
  i1 = i0+1;
  j0 = interval[1];
  j1 = j0+1;

  dx = x[i1] - x[i0];
  if (dx <= 0.)
    return 1.e10;

  dy = y[j1] - y[j0];
  if (dy <= 0.)
    return 1.e10;

  k00 = i0 + nx*j0;
  k01 = i0 + nx*j1;
  k10 = i1 + nx*j0;
  k11 = i1 + nx*j1;

  /* Interpolation coefficients for y direction */

  b = (yy - y[j0])/dy;
  a = 1.-b;
  g = dy*dy/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  /* Interpolate in y direction for f(x[], y) to get f at points U and V */

  f0 = a*f[k00] + b*f[k01] + c*sy[k00] + d*sy[k01];
  f1 = a*f[k10] + b*f[k11] + c*sy[k10] + d*sy[k11];

  /* Interpolate in y direction for 2nd derivative of f(x[], y) at U and V */

  s0 = a*sx[k00] + b*sx[k01] + c*sxy[k00] + d*sxy[k01];
  s1 = a*sx[k10] + b*sx[k11] + c*sxy[k10] + d*sxy[k11];

  /* Interpolate in x direction for value of function at our point P */

  b = (xx - x[i0])/dx;
  a = 1.-b;
  g = dx*dx/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  return a*f0 + b*f1 + c*s0 + d*s1;
}

/*
 * Free memory allocated by csp3d_open()
 */
void csp3d_close(double **splines)
{
  int i;

  if (splines == NULL)
    return;
  for (i = 0; i < 7; i++) {
    if (splines[i] != NULL) {
      free(splines[i]);
      splines[i] = NULL;
    }
  }
  free(splines);
  splines = NULL;
}

double **
csp3d_open(double *xv[], int nv[], double f[])
/*
 * Purpose
 *
 *   Calculate 3-D cubic spline coefficients for a function f(x, y, z)
 *   of three independent variables tabulated on a 3-D rectilinear grid.
 *
 *   To evaluate this 3-D spline, use csp3d_eval()
 *
 * Input
 *
 *   double *xv[]  Array of 3 arrays of independent variables;
 *                 xv[] is the array of {x[], y[], z[]}
 *   int nv[]      Array of length 3 containing lengths of xv[];
 *                 nv[] = {nx, ny, nz}
 *   double f[]    f(x[i], y[j], z[k]) packed into a one-dimensional
 *                 array of length nx*ny*nz. The z[k] index varies the
 *                 fastest. That is, instead of f[i][j][k], we use
 *                 f[(i*ny+j)*nz+k]
 * Returns
 *
 *   double **splines  Spline coefficient array. This is an array of length
 *                     7 of arrays of length nx*ny*nz. The user passes
 *                     this to csp3d_eval(). Use csp3d_close(splines) to free
 *                     this memory.
 *
 *                     On error, returns NULL. Errors is either memory
 *                     allocation or xv[][] not in strictly increasing order.
 *
 * Method
 *
 *   This calls csp1d() to calculate a 1-D cubic spline through every
 *   grid point parallel to each axis in turn. The csp1d() calculates the
 *   second derivative at each input point. These spline coefficients are
 *   placed in arrays the same size as f[]. Cascaded cubic splines of these
 *   spline coefficeints are calculated versus directions perpendicular to
 *   the original independent variable. See csp3d_eval() to see how these are
 *   used to evaluate the spline.
 *
 * By
 *
 *   Robert W. Strickland, September 21, 2015
 */
{
  int nw;
  double *w, *s;
  int i, j, k;
  double *sx, *sy, *sz, *sxy, *sxz, *syz, *sxyz;
  double **splines;
  int nf;
  double *x, *y, *z;
  int nx, ny, nz;
  int nerr;

  x = xv[0];
  y = xv[1];
  z = xv[2];

  nx = nv[0];
  ny = nv[1];
  nz = nv[2];
 
  nf = nx*ny*nz;
  splines = (double **) malloc(7*sizeof(double *));
  if (splines == NULL)
    return NULL;
  for (i = 0; i < 7; i++) {
    splines[i] = (double *) malloc(nf*sizeof(double));
    if (splines[i] == NULL)
      return NULL;
  }
  sx   = splines[0];
  sy   = splines[1];
  sz   = splines[2];
  sxy  = splines[3];
  sxz  = splines[4];
  syz  = splines[5];
  sxyz = splines[6];

  /*
   * Allocate work arrays that are the length of largest of nx, ny, nz
   */
  nw = (nx > ny) ? nx : ny;
  nw = (nz > nw) ? nz : nw;
  w = (double *) malloc(nw * sizeof(double));
  s = (double *) malloc(nw * sizeof(double));
  if (w == NULL || s == NULL) {
    for (i = 0; i < 7; i++) {
      free(splines[i]);
      free(splines);
    }
    return NULL;
  }

  nerr = 0;

  /* Spline of f vs x -> sx */

  for (j = 0; j < ny; j++) {
    for (k = 0; k < nz; k++) {
      for (i = 0; i < nx; i++)
        w[i] = f[(i*ny+j)*nz+k];
      nerr += csp1d(x, w, nx, NATSPLN, NATSPLN, s);
      for (i = 0; i < nx; i++)
        sx[(i*ny+j)*nz+k] = s[i];
    }
  }

  /* Spline of: f vs y -> sy
   * Cascaded spline of sx vs y -> sxy
   */

  for (i = 0; i < nx; i++) {
    for (k = 0; k < nz; k++) {
      for (j = 0; j < ny; j++)
        w[j] = f[(i*ny+j)*nz+k];
      nerr += csp1d(y, w, ny, NATSPLN, NATSPLN, s);
      for (j = 0; j < ny; j++)
        sy[(i*ny+j)*nz+k] = s[j];

      for (j = 0; j < ny; j++)
        w[j] = sx[(i*ny+j)*nz+k];
      nerr += csp1d(y, w, ny, NATSPLN, NATSPLN, s);
      for (j = 0; j < ny; j++)
        sxy[(i*ny+j)*nz+k] = s[j];
    }
  }

  /* Spline of f vs. z -> sz
   * Cascaded splines of: sx vs z -> sxz
   *                      sy vs z -> syz
   *                      sxy vs z -> sxyz
   */

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++)
        w[k] = f[(i*ny+j)*nz+k];
      nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
      for (k = 0; k < nz; k++)
        sz[(i*ny+j)*nz+k] = s[k];

      for (k = 0; k < nz; k++)
        w[k] = sx[(i*ny+j)*nz+k];
      nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
      for (k = 0; k < nz; k++)
        sxz[(i*ny+j)*nz+k] = s[k];

      for (k = 0; k < nz; k++)
        w[k] = sy[(i*ny+j)*nz+k];
      nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
      for (k = 0; k < nz; k++)
        syz[(i*ny+j)*nz+k] = s[k];

      for (k = 0; k < nz; k++)
        w[k] = sxy[(i*ny+j)*nz+k];
      nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
      for (k = 0; k < nz; k++)
        sxyz[(i*ny+j)*nz+k] = s[k];
    }
  }

  free(s);
  free(w);
  if (nerr != 0) {
    csp3d_close(splines);
    splines = NULL;
  }
  return splines;
}

double
csp3d_eval(double p[], double *xv[], int nv[], double f[],
    double **splines, int interval[])
/*
 * Purpose
 *
 *   Interpolate a smooth function of three variables tabulated on a
 *   rectilinear grid. Requires csp3d_open() to calculate coefficients
 *   stored in splines. There is just one call to csp3d_open() ahead of
 *   time and then multiple calls to csp3d_eval() to evaluate the spline
 *   at each point.
 *
 * Input
 *
 *   double p[3]  The point P(x, y, z) where we want f interpolated
 *   double xv[3] The array of {x[], y[], z[]} where f is tabulated.
 *   int nv[3]    The array of {nx, ny, nz}, the lengths of the 
 *                {x[], y[], z[]} arrays
 *   double f[]   f(x[i], y[j], z[k]) packed into a one-dimensional
 *                array of length nx*ny*nz. The z[k] index varies the
 *                fastest. That is, instead of f[i][j][k], we use
 *                f[(i*ny+j)*nz+k]
 *   double **splines
 *                Spline coefficients returned by csp3d_open()
 *   int interval[3]  Work array of length 3 to keep track of interval. If
 *                there are several parallel calls at the same value of p[],
 *
 * Method
 *   
 *   We start by finding the 3-D box in our grid that contains P(x,y,z).
 *   Find i, j, k such that
 *
 *     x[i] < x < x[i+1]
 *     y[j] < y < y[j+1]
 *     z[k] < z < z[k+1]
 *
 *   Next we interpolate along the 4 edges of the 3-D box parallel to the
 *   z axis to find a rectangle, eliminating the z dependence. Next we
 *   interpolate along the 2 edges of the rectangle parallel to the y axis
 *   to get a line segment, eliminating the y dependence. Finally, we
 *   interpolate along the line segment to get a single point and to eliminate
 *   the x dependence. This can be considered the process of repeatedly
 *   slicing perpendicular to one of the axes, reducing the dimension by one
 *   each time until a single point remains.
 *
 *   To get a better understanding of the underlying process it might be
 *   best to read this code from the bottom up.
 *
 * By
 *
 *   Robert W. Strickland, September 21, 2015
 */
{
  double g, a, b, c, d;
  int i0, i1, j0, j1, k0, k1;
  int i8[8];
  double dx, dy, dz;
  double f4[4], sx4[4], sy4[4], sxy4[4], f0, f1, sx0, sx1;
  int n, n0, n1;
  double *sx, *sy, *sz, *sxy, *sxz, *syz, *sxyz;
  double *x, *y, *z, xx, yy, zz;
  int ny, nz;
  
  xx = p[0];
  yy = p[1];
  zz = p[2];

  ny = nv[1];
  nz = nv[2];

  sx   = splines[0];
  sy   = splines[1];
  sz   = splines[2];
  sxy  = splines[3];
  sxz  = splines[4];
  syz  = splines[5];
  sxyz = splines[6];

  /* Linear search for interval containing each independent variable */

  for (n = 0; n < 3; n++) {
    if (interval[n] < 0 || interval[n] > nv[n]-2)
      interval[n] = nv[n]/2;
    while (interval[n] > 0 && p[n] < xv[n][interval[n]])
      interval[n]--;
    while (interval[n] < nv[n]-2 && p[n] > xv[n][interval[n]+1])
      interval[n]++;
  }

  i0 = interval[0];
  i1 = i0+1;
  j0 = interval[1];
  j1 = j0+1;
  k0 = interval[2];
  k1 = k0+1;

  /* Get indices for 8 corners of box containing (xx, yy, zz). The index
   * that goes from 0-7 can be considered a 3-bit field where the most sig-
   * nificant bit corresponsds to the x coordinate and the least significant
   * bit corresponds to the z coordinate
   */

  i8[0] = (i0*ny+j0)*nz+k0;
  i8[1] = (i0*ny+j0)*nz+k1;
  i8[2] = (i0*ny+j1)*nz+k0;
  i8[3] = (i0*ny+j1)*nz+k1;
  i8[4] = (i1*ny+j0)*nz+k0;
  i8[5] = (i1*ny+j0)*nz+k1;
  i8[6] = (i1*ny+j1)*nz+k0;
  i8[7] = (i1*ny+j1)*nz+k1;

  /* Interpolation coefficients for z direction */

  z = xv[2];
  dz = z[k1]-z[k0];
  if (dz <= 0.)
    return 1.e30;
  b = (zz - z[k0])/dz;
  a = 1.-b;
  g = dz*dz/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  /* Interpolate in z direction for corners of a rectangular slice through
   * z=zz plane containing point P(xx,yy,zz). We'll find it convenient to
   * label the points on this slice as follows:
   *
   *  y[j+1] -01-----+---11
   *           |     |    |
   *      yy --U-----P----V
   *           |     |    |
   *           |     |    |
   *    y[j] -00-----+---10
   *           |     |    |
   *         x[i]   xx  x[i+1]
   *
   * The loop from 0 to 3 goes through these 4 points in binary counting order
   */

  for (n = 0; n < 4; n++) {
    n0 = i8[2*n];
    n1 = i8[2*n+1];
    f4  [n] = a*  f[n0] + b*  f[n1] + c*  sz[n0] + d*  sz[n1];
    sx4 [n] = a* sx[n0] + b* sx[n1] + c* sxz[n0] + d* sxz[n1];
    sy4 [n] = a* sy[n0] + b* sy[n1] + c* syz[n0] + d* syz[n1];
    sxy4[n] = a*sxy[n0] + b*sxy[n1] + c*sxyz[n0] + d*sxyz[n1];
  }

  /* Interpolation coefficients for y direction */

  y = xv[1];
  dy = y[j1]-y[j0];
  if (dy <= 0.)
    return 1.e30;
  b = (yy - y[j0])/dy;
  a = 1.-b;
  g = dy*dy/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  /* Interpolate in y direction for f(U) and f(V) */

  f0  = a*f4[0] + b*f4[1] + c*sy4[0] + d*sy4[1];
  f1  = a*f4[2] + b*f4[3] + c*sy4[2] + d*sy4[3];

  /* Interpolate in y direction for sx(U) and sx(V) */

  sx0 = a*sx4[0] + b*sx4[1] + c*sxy4[0] + d*sxy4[1];
  sx1 = a*sx4[2] + b*sx4[3] + c*sxy4[2] + d*sxy4[3];

  /* Interpolate in x direction for value of function at P=(xx, yy) */

  x = xv[0];
  dx = x[i1]-x[i0];
  if (dx <= 0.)
    return 1.e30;
  b = (xx - x[i0])/dx;
  a = 1.-b;
  g = dx*dx/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  return a*f0 + b*f1 + c*sx0 + d*sx1;
}


/*
 * Free memory allocated by csp4d_open()
 */
void
csp4d_close(double **splines)
{
  int i;
  for (i = 0; i < 15; i++) {
    if (splines[i] != NULL) {
      free(splines[i]);
      splines[i] = NULL;
    }
  }
  free(splines);
  splines = NULL;
  return;
}

double **
csp4d_open(double *xv[], int nv[], double f[])
/*
 * Purpose
 *
 *   Calculate 4-D cubic spline coefficients for a smooth function
 *   f(x, y, z, t) four independent variables tabulated on a 4-D
 *   rectilinear grid.
 *
 *   To evaluate this 4-D spline, use csp4d_eval()
 *
 * Input
 *
 *   double *xv[]  Array of 4 arrays of independent variables;
 *                 xv[] is the array of {x[], y[], z[], t[]}
 *   int nv[]      Array of length 4 containing lengths of xv[];
 *                 nv[] = {nx, ny, nz, nt}
 *   double f[]    f(x[i], y[j], z[k], t[l]) packed into a one-dimensional
 *                 array of length nx*ny*nz*nt. The t[l] index varies the
 *                 fastest. That is, instead of f[i][j][k][l], we use
 *                 f[((i*ny+j)*nz+k)*nt+l]
 * Returns
 *
 *   double **splines  Spline coefficient array. This is an array of length
 *                     15 of arrays of length nx*ny*nz*nt. The user passes
 *                     this to csp4d_eval(). Use csp4d_close(splines) to free
 *                     this memory.
 * Method
 *
 *   This calls csp1d() to calculate a 1d spline along lines through
 *   every grid point parallel to in every axis in turn. The csp1d()
 *   calcualtes the second derivative at each grid point. These spline
 *   coefficients are stored in arrays the same size as f[]. Then cascades
 *   of splines are calculated in directions perpendicular to the original
 *   splines. These coefficients are stored. See csp34_eval() to see how
 *   these are used to evaluate the spline.
 *
 * By
 *
 *   Robert W. Strickland, September 21, 2015
 */
{
  int nx, ny, nz, nt, nw;
  double *w, *s;
  int i, j, k, l;
  int nf;
  /*
   * Here are our 2nd derivative 1-D cubic spline coefficients.
   * For cascaded splines, we'll use the following notation. The rightmost
   * character indicates the independent variable for the spline. A spline of
   * a spline adds the independent variable to the end. Thus, szxy would be
   *
   * spline of {
   *   spline of {
   *     spline of {
   *       f vs z
   *     }
   *   }
   *   vs x
   * }
   * vs y
   */
  double *sx, *sy, *sz, *st, *szt, *sxt, *sxy, *syt, *syz, *syzt, *sxz,
         *sxzt, *sxyt, *sxyz, *sxyzt;
  double **splines;
  double *x, *y, *z, *t;
  int nerr;

  nx = nv[0];
  ny = nv[1];
  nz = nv[2];
  nt = nv[3];

  x = xv[0];
  y = xv[1];
  z = xv[2];
  t = xv[3];

  /* Allocate memory for 4-D spline coefficients. This is an array of
   * 15 D arrays
   */

  nf = nx*ny*nz*nt;

  splines = (double **) malloc(15*sizeof(double*));
  if (splines == NULL)
    return NULL;
  for (i = 0; i < 15; i++)
    splines[i] = NULL;
  for (i = 0; i < 15; i++) {
    splines[i] = (double *) malloc(nf*sizeof(double));
    if (splines[i] == NULL) {
      csp4d_close(splines);
      return NULL;
    }
  }

  /* The following order must match that in csp4d_eval() */

  st    = splines[0];
  sz    = splines[1];
  sy    = splines[2];
  sx    = splines[3];
  sxt   = splines[4];
  sxz   = splines[5];
  sxy   = splines[6];
  syz   = splines[7];
  syt   = splines[8];
  szt   = splines[9];
  syzt  = splines[10];
  sxzt  = splines[11];
  sxyt  = splines[12];
  sxyz  = splines[13];
  sxyzt = splines[14];

  /*
   * Allocate work arrays that are as the length of largest of nx, ny,
   * nz, nt. These hold 1-D cubic spline coefficients temporarily until
   * they can be copied into the splines array.
   */
  nw = nv[0];
  for (i = 1; i < 4; i++)
    if (nv[i] > nw)
      nw = nv[i];
  w = (double *) malloc(nw * sizeof(double));
  s = (double *) malloc(nw * sizeof(double));
  if (w == NULL || s == NULL) {
    csp4d_close(splines);
    return NULL;
  }

  nerr = 0;

  /* Cubic spline of f vs. x -> sx */

  for (j = 0; j < ny; j++) {
    for (k = 0; k < nz; k++) {
      for (l = 0; l < nt; l++) {
        for (i = 0; i < nx; i++)
          w[i] = f[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(x, w, nx, NATSPLN, NATSPLN, s);
        for (i = 0; i < nx; i++)
          sx[((i*ny+j)*nz+k)*nt+l] = s[i];
      }
    }
  }

  /*    Cubic spline of  f vs. y -> sy
   * cascaded spline of sx vs. y -> sxy
   */

  for (i = 0; i < nx; i++) {
    for (k = 0; k < nz; k++) {
      for (l = 0; l < nt; l++) {
        for (j = 0; j < ny; j++)
          w[j] = f[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(y, w, ny, NATSPLN, NATSPLN, s);
        for (j = 0; j < ny; j++)
          sy[((i*ny+j)*nz+k)*nt+l] = s[j];
        for (j = 0; j < ny; j++)
          w[j] = sx[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(y, w, ny, NATSPLN, NATSPLN, s);
        for (j = 0; j < ny; j++)
          sxy[((i*ny+j)*nz+k)*nt+l] = s[j];
      }
    }
  }

  /* 
   *     Cubic spline of:  f vs. z -> sz
   * Cascaded splines of: sy vs. z -> syz
   *                      sx vs. z -> sxz
   *                     sxy vs. z -> sxyz
   */
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (l = 0; l < nt; l++) {
        for (k = 0; k < nz; k++)
          w[k] = f[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
        for (k = 0; k < nz; k++)
          sz[((i*ny+j)*nz+k)*nt+l] = s[k];
        for (k = 0; k < nz; k++)
          w[k] = sy[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
        for (k = 0; k < nz; k++)
          syz[((i*ny+j)*nz+k)*nt+l] = s[k];
        for (k = 0; k < nz; k++)
          w[k] = sx[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
        for (k = 0; k < nz; k++)
          sxz[((i*ny+j)*nz+k)*nt+l] = s[k];
        for (k = 0; k < nz; k++)
          w[k] = sxy[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
        for (k = 0; k < nz; k++)
          sxyz[((i*ny+j)*nz+k)*nt+l] = s[k];
      }
    }
  }

  /*     Cubic spline of:   f vs. t -> st
   * Cascaded splines of:  sx vs. t -> sxt
   *                       sy vs. t -> syt
   *                       sz vs. t -> szt
   *                      sxy vs. t -> sxyt
   *                      syz vs. t -> syzt
   *                      sxz vs. t -> sxzt
   *                     sxyz vs. t -> sxyzt
   */
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        for (l = 0; l < nt; l++)
          w[l] = f[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
        for (l = 0; l < nt; l++)
          st[((i*ny+j)*nz+k)*nt+l] = s[l];
        for (l = 0; l < nt; l++)
          w[l] = sx[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
        for (l = 0; l < nt; l++)
          sxt[((i*ny+j)*nz+k)*nt+l] = s[l];
        for (l = 0; l < nt; l++)
          w[l] = sy[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
        for (l = 0; l < nt; l++)
          syt[((i*ny+j)*nz+k)*nt+l] = s[l];
        for (l = 0; l < nt; l++)
          w[l] = sz[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
        for (l = 0; l < nt; l++)
          szt[((i*ny+j)*nz+k)*nt+l] = s[l];
        for (l = 0; l < nt; l++)
          w[l] = sxy[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s); 
        for (l = 0; l < nt; l++)
          sxyt[((i*ny+j)*nz+k)*nt+l] = s[l];
        for (l = 0; l < nt; l++)
          w[l] = syz[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s); 
        for (l = 0; l < nt; l++)
          syzt[((i*ny+j)*nz+k)*nt+l] = s[l];
        for (l = 0; l < nt; l++)
          w[l] = sxz[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
        for (l = 0; l < nt; l++)
          sxzt[((i*ny+j)*nz+k)*nt+l] = s[l];
        for (l = 0; l < nt; l++)
          w[l] = sxyz[((i*ny+j)*nz+k)*nt+l];
        nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
        for (l = 0; l < nt; l++)
          sxyzt[((i*ny+j)*nz+k)*nt+l] = s[l];
      }
    }
  }
  free(w);
  free(s);
  if (nerr != 0) {
    csp4d_close(splines);
    splines = NULL;
  }
  return splines;
}

double
csp4d_eval(double p[], double *xv[], int nv[], double f[],
    double **splines, int interval[])
/*
 * Purpose
 *
 *   Interpolate a smooth function of four variables tabulated on a
 *   rectilinear grid. We have a function f(x, y, z, t) tabulated at fixed
 *   values of x, y, z and t. Requires csp4d_open() to calculate coefficients
 *   stored in splines. There is just one call to csp4d_open() ahead of time
 *   and then multiple calls to csp4d_eval() to evaluate the spline at each
 *   point.
 *
 * Input
 *
 *   double p[4]  The point P(x, y, z, t) where we want f interpolated
 *   double xv[4] The array of {x[], y[], z[], t[]} where f is tabulated.
 *   int nv[4]    The array of {nx, ny, nz, nt}, the lengths of the 
 *                {x[], y[], z[], t[]} arrays
 *   double f[]   f(x[i], y[j], z[k], t[l]) packed into a one-dimensional
 *                array of length nx*ny*nz*nt. The t[l] index varies the
 *                fastest. That is,
 *   double **splines
 *                Spline coefficients returned by csp4d_open()
 *   int interval[4]
 *                Work array of length 4 to keep track of interval. If
 *                there are several parallel calls at the same value of p[],
 *                all should share the same interval[] array.
 *
 * Method
 *   
 *   We start by finding the 4-D box in our grid that contains P.
 *
 *     x[i] < x < x[i+1]
 *     y[j] < y < y[j+1]
 *     z[k] < z < z[k+1]
 *     t[l] < t < t[l+1]
 *
 *   Then we interpolate our 1-D splines along the 8 edges of the 4-D box
 *   parallel to the t axis to eliminate the t dependence. This gives us
 *   a 3-D box. Next we interpolate along the 4 edges of the 3-D box
 *   parallel to the z axis to find a rectangle, eliminating the z
 *   dependence. Next we interpolate along the 2 edges of the rectangle
 *   parallel to the y axis to get a line segment, eliminating the y
 *   dependence. Finally, we interpolate along the line segment to get a
 *   single point and to eliminate the x dependence. This can be considered
 *   the process of repeatedly slicing perpendicular to one of the axes,
 *   reducing the dimension by one each time.
 *
 *   To get a better understanding, of the underlying process it might be
 *   best to read this code from the bottom up.
 *
 * By
 *
 *   Robert W. Strickland, September 21, 2015
 */
{
  int i16[16];
  double f8[8], sy8[8], sx8[8], sxy8[8], sz8[8], sxz8[8], syz8[8], sxyz8[8];
  double f4[4], sy4[4], sx4[4], sxy4[4];
  double f0, f1, s0x, s1x;

  double g, a, b, c, d;
  int i0, i1, j0, j1, k0, k1, l0, l1;
  double dx, dy, dz, dt;
  double *sx, *sy, *sz, *st, *szt, *sxt, *sxy, *syt, *syz, *syzt, *sxz,
         *sxzt, *sxyt, *sxyz, *sxyzt;
  int n, n0, n1;
  int ny, nz, nt;
  double *x, *y, *z, *t;

  /* Linear search for interval containing each independent variable */

  for (n = 0; n < 4; n++) {
    if (interval[n] < 0 || interval[n] > nv[n]-2)
      interval[n] = nv[n]/2;
    while (interval[n] > 0 && p[n] < xv[n][interval[n]])
      interval[n]--;
    while (interval[n] < nv[n]-2 && p[n] > xv[n][interval[n]+1])
      interval[n]++;
  }

  /* Get indices for of corners of 4-D box containing (xx, yy, zz, tt) */

  i0 = interval[0];
  i1 = i0+1;
  j0 = interval[1];
  j1 = j0+1;
  k0 = interval[2];
  k1 = k0+1;
  l0 = interval[3];
  l1 = l0+1;

  ny = nv[1];
  nz = nv[2];
  nt = nv[3];

  i16[ 0] = ((i0*ny+j0)*nz+k0)*nt+l0;
  i16[ 1] = ((i0*ny+j0)*nz+k0)*nt+l1;
  i16[ 2] = ((i0*ny+j0)*nz+k1)*nt+l0;
  i16[ 3] = ((i0*ny+j0)*nz+k1)*nt+l1;
  i16[ 4] = ((i0*ny+j1)*nz+k0)*nt+l0;
  i16[ 5] = ((i0*ny+j1)*nz+k0)*nt+l1;
  i16[ 6] = ((i0*ny+j1)*nz+k1)*nt+l0;
  i16[ 7] = ((i0*ny+j1)*nz+k1)*nt+l1;
  i16[ 8] = ((i1*ny+j0)*nz+k0)*nt+l0;
  i16[ 9] = ((i1*ny+j0)*nz+k0)*nt+l1;
  i16[10] = ((i1*ny+j0)*nz+k1)*nt+l0;
  i16[11] = ((i1*ny+j0)*nz+k1)*nt+l1;
  i16[12] = ((i1*ny+j1)*nz+k0)*nt+l0;
  i16[13] = ((i1*ny+j1)*nz+k0)*nt+l1;
  i16[14] = ((i1*ny+j1)*nz+k1)*nt+l0;
  i16[15] = ((i1*ny+j1)*nz+k1)*nt+l1;

  /* Retrieve the spline pointers from splines array. The order here must
   * match that in csp4d_open()
   */
  st    = splines[0];
  sz    = splines[1];
  sy    = splines[2];
  sx    = splines[3];
  sxt   = splines[4];
  sxz   = splines[5];
  sxy   = splines[6];
  syz   = splines[7];
  syt   = splines[8];
  szt   = splines[9];
  syzt  = splines[10];
  sxzt  = splines[11];
  sxyt  = splines[12];
  sxyz  = splines[13];
  sxyzt = splines[14];

  /* Interpolate in the t direction to get all 8 corners of the 3-D box
   * Our point P=(xx, yy, zz, tt) lies in the 4-D box
   *
   *   x[i] < xx < x[i+1]
   *   y[j] < yy < y[j+1]
   *   z[k] < zz < z[k+1]
   *   t[l] < tt < t[l+1]
   *
   * We can index the 16 corners of this box with a 4-bit binary index
   * where the most significant bit refers to the x coordinate and the
   * least significant bit refers to the t coordinate. For instance, index
   * 0110 refers to the point (x[i],y[j+1],z[k+1],t[l]).
   *
   * First we'll interpolate in the t direction. We'll interpolate all
   * of the values of f and all of the splines we'll need for later steps.
   * After we interpolate for t, we will have values on the corners of
   * a 3-D box. We can drop the least significant bit (the t index) and
   * use a 3-bit index.
   *
   *     011+--------+------------+111
   *        |\   z   |\           |\
   *        | \  ^   | \          | \
   *        |  +-|------+------------+
   *        |  |\|   |  |\        |  |\
   *        | 001+--------+------------+101
   *        +--|-|---+--|-|-------+  | |
   *        |\.|.|....\.|.|........\ | |
   *        | \|.|.....\|.|.........\| |
   *        |  U-|------P-|----------V |
   *        |  |\|.......\|...........\|
   *      y |  | +--------+------------+
   *       \|  | |   |  | |       |  | |
   *     010+--|-|---+--|-|-------+110 |
   *         \ | |    \ | |        \ | |
   *          \| |     \| |         \| |
   *           +-|------+-|----------+ |
   *            \|       \|           \|
   *          000+--------+------------+100->x
   *
   * Again, the msb refers to the x coordinate.
   */

  t = xv[3];
  dt = t[l1]-t[l0];
  if (dt <= 0.)
    return 1.e30;
  b = (p[3] - t[l0])/dt;
  a = 1.-b;
  g = dt*dt/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  /* The n index goes through all 8 corners of the box from binary 000 to 111
   * (decimal 7) in binary counting order.
   *
   * n0 is the index into the f[] array (or one of its splines) with
   * tt=t[l]. It is the 4-bit index to a corner of the the 4-D box
   * containing our point.
   *
   * n1 is the index into the f[] array (or one of its splines) with
   * tt=t[l+1]. It, too is a 4-bit corner index.
   */
  for (n = 0; n < 8; n++) {
    n0 = i16[2*n];
    n1 = i16[2*n+1];
    f8   [n] = a*   f[n0] + b*   f[n1] + c*   st[n0] + d*   st[n1];
    sx8  [n] = a*  sx[n0] + b*  sx[n1] + c*  sxt[n0] + d*  sxt[n1];
    sy8  [n] = a*  sy[n0] + b*  sy[n1] + c*  syt[n0] + d*  syt[n1];
    sz8  [n] = a*  sz[n0] + b*  sz[n1] + c*  szt[n0] + d*  szt[n1];
    sxy8 [n] = a* sxy[n0] + b* sxy[n1] + c* sxyt[n0] + d* sxyt[n1];
    syz8 [n] = a* syz[n0] + b* syz[n1] + c* syzt[n0] + d* syzt[n1];
    sxz8 [n] = a* sxz[n0] + b* sxz[n1] + c* sxzt[n0] + d* sxzt[n1];
    sxyz8[n] = a*sxyz[n0] + b*sxyz[n1] + c*sxyzt[n0] + d*sxyzt[n1];
  }

  /* Interpolate in the z direction along the edges of the 3-box to get
   * corners of the rectangle parallel to the x-z plane at zz (the shaded
   * plane above)
   *
   *  y[j+1] -01-----+---11
   *           |     |    |
   *      yy --U-----P----V
   *           |     |    |
   *           |     |    |
   *    y[j] -00-----+---10
   *           |     |    |
   *         x[i]   xx  x[i+1]
   */

  z = xv[2];
  dz = z[k1]-z[k0];
  if (dz <= 0.)
    return 1.e30;
  b = (p[2] - z[k0])/dz;
  a = 1.-b;
  g = dz*dz/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  /* The index n refers to the corners of the rectangle. Here we have dropped
   * the z index, and n takes on all the values from binary 00 to 11. The
   * msb still refers to x.
   *
   * n0 refers to the bottom z[k] plane of our 3-D box, and
   * n1 refers to the top z[k+1] plane.
   */
  for (n = 0; n < 4; n++) {
    n0 = 2*n;
    n1 = 2*n+1;
    f4  [n] = a*  f8[n0] + b*  f8[n1] + c*  sz8[n0] + d*  sz8[n1];
    sx4 [n] = a* sx8[n0] + b* sx8[n1] + c* sxz8[n0] + d* sxz8[n1];
    sy4 [n] = a* sy8[n0] + b* sy8[n1] + c* syz8[n0] + d* syz8[n1];
    sxy4[n] = a*sxy8[n0] + b*sxy8[n1] + c*sxyz8[n0] + d*sxyz8[n1];
  }

  /* Interpolate in the y direction across the rectangle to get function
   * and spline at the endpoints of our line segment from points U to V
   * (here indexed 0 and 1)
   */
  y = xv[1];
  dy = y[j1]-y[j0];
  if (dy <= 0.)
    return 1.e30;
  b = (p[1] - y[j0])/dy;
  a = 1.-b;
  g = dy*dy/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  f0  = a* f4[0] + b* f4[1] + c* sy4[0] + d* sy4[1];
  f1  = a* f4[2] + b* f4[3] + c* sy4[2] + d* sy4[3];
  s0x = a*sx4[0] + b*sx4[1] + c*sxy4[0] + d*sxy4[1];
  s1x = a*sx4[2] + b*sx4[3] + c*sxy4[2] + d*sxy4[3];

  /* Interpolate in the x direction along the line segment to get f(P) */

  x = xv[0];
  dx = x[i1]-x[i0];
  if (dx <= 0.)
    return 1.e30;
  b = (p[0] - x[i0])/dx;
  a = 1.-b;
  g = dx*dx/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  return a*f0 + b*f1 + c*s0x + d*s1x;
}

/*
 * Free memory allocated by csp5d_open()
 */
void
csp5d_close(double **splines)
{
  int i;
  for (i = 0; i < 29; i++) {
    if (splines[i] != NULL) {
      free(splines[i]);
      splines[i] = NULL;
    }
  }
  free(splines);
  splines = NULL;
  return;
}

double **
csp5d_open(double *xv[], int nv[], double f[])
/*
 * Purpose
 *
 *   Calculate 5-D cubic spline coefficients for a smooth function
 *   f(x, y, z, t) four independent variables tabulated on a 5-D
 *   rectilinear grid.
 *
 *   To evaluate this 5-D spline, use csp5d_eval()
 *
 * Input
 *
 *   double *xv[]  Array of 5 arrays of independent variables;
 *                 xv[] is the array of {x[], y[], z[], t[]}
 *   int nv[]      Array of length 5 containing lengths of xv[];
 *                 nv[] = {nx, ny, nz, nt, nu}
 *   double f[]    f(x[i], y[j], z[k], t[l], u[m]) packed into a one-
 *                 dimensional array of length nx*ny*nz*nt*nu. The u[m]
 *                 index varies the fastest. That is, instead of
 *                 f[i][j][k][l][m], we use
 *                 f[(((i*ny+j)*nz+k)*nt+l)*nu+m]
 * Returns
 *
 *   double **splines  Spline coefficient array. This is an array of length
 *                     29 of arrays of length nx*ny*nz*nt*nu. The user passes
 *                     this to csp5d_eval(). Use csp5d_close(splines) to free
 *                     this memory.
 * Method
 *
 *   This calls csp1d() to calculate a 1d spline along lines through
 *   every grid point parallel to in every axis in turn. The csp1d()
 *   calcualtes the second derivative at each grid point. These spline
 *   coefficients are stored in arrays the same size as f[]. Then cascades
 *   of splines are calculated in directions perpendicular to the original
 *   splines. These coefficients are stored. See csp5d_eval() to see how
 *   these are used to evaluate the spline.
 *
 * By
 *
 *   Robert W. Strickland, September 29, 2015
 */
{
  int nx, ny, nz, nt, nu, nw;
  double *w, *s;
  int i, j, k, l, m;
  int nf;
  /*
   * Here are our 2nd derivative 1-D cubic spline coefficients.
   * For cascaded splines, we'll use the following notation. The rightmost
   * character indicates the independent variable for the spline. A spline of
   * a spline adds the independent variable to the end. Thus, szxy would be
   *
   * spline of {
   *   spline of {
   *     spline of {
   *       f vs z
   *     }
   *   }
   *   vs x
   * }
   * vs y
   */
  double *sx;
  double *sy, *sxy;
  double *sz, *sxz, *syz, *sxyz;
  double *st, *syt, *sxyt, *szt, *sxzt, *syzt, *sxyzt;
  double *su, *sxu, *syu, *sxyu, *szu, *sxzu, *syzu, *sxyzu,
         *stu, *sytu, *sxytu, *sztu, *sxztu, *syztu, *sxyztu;
  double *x, *y, *z, *t, *u;
  double **splines;
  int nerr;

  nx = nv[0];
  ny = nv[1];
  nz = nv[2];
  nt = nv[3];
  nu = nv[4];

  x = xv[0];
  y = xv[1];
  z = xv[2];
  t = xv[3];
  u = xv[4];

  /* Allocate memory for 5-D spline coefficients. This is an array of
   * 29 arrays
   */

  nf = nx*ny*nz*nt*nu;

  splines = (double **) malloc(29*sizeof(double*));
  if (splines == NULL)
    return NULL;
  for (i = 0; i < 29; i++)
    splines[i] = NULL;
  for (i = 0; i < 29; i++) {
    splines[i] = (double *) malloc(nf*sizeof(double));
    if (splines[i] == NULL) {
      csp5d_close(splines);
      return NULL;
    }
  }

  /* Retrieve the spline pointers from splines array. The order here must
   * match that in csp5d_eval()
   */
  sx = splines[0];
  sy = splines[1];
  sxy = splines[2];
  sz = splines[3];
  syz = splines[4];
  sxz = splines[5];
  sxyz = splines[6];
  st = splines[7];
  syt = splines[8];
  szt = splines[9];
  sxyt = splines[10];
  syzt = splines[11];
  sxzt = splines[12];
  sxyzt = splines[13];
  su = splines[14];
  stu = splines[15];
  sxu = splines[16];
  syu = splines[17];
  sytu = splines[18];
  szu = splines[19];
  sztu = splines[20];
  sxyu = splines[21];
  sxytu = splines[22];
  syzu = splines[23];
  syztu = splines[24];
  sxzu = splines[25];
  sxztu = splines[26];
  sxyzu = splines[27];
  sxyztu = splines[28];

  /*
   * Allocate work arrays that are as the length of largest of nx, ny,
   * nz, nt, nu. These hold 1-D cubic spline coefficients temporarily until
   * they can be copied into the splines array.
   */
  nw = nv[0];
  for (i = 1; i < 5; i++)
    if (nv[i] > nw)
      nw = nv[i];
  w = (double *) malloc(nw * sizeof(double));
  s = (double *) malloc(nw * sizeof(double));
  if (w == NULL || s == NULL) {
    csp5d_close(splines);
    return NULL;
  }

  nerr = 0;

  /* Cubic spline of f vs. x -> sx */

  for (j = 0; j < ny; j++) {
    for (k = 0; k < nz; k++) {
      for (l = 0; l < nt; l++) {
        for (m = 0; m < nu; m++) {
          for (i = 0; i < nx; i++)
            w[i] = f[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(x, w, nx, NATSPLN, NATSPLN, s);
          for (i = 0; i < nx; i++)
            sx[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[i];
        }
      }
    }
  }

  /*    Cubic spline of:  f vs. y -> sy
   * cascaded spline of: sx vs. y -> sxy
   */

  for (i = 0; i < nx; i++) {
    for (k = 0; k < nz; k++) {
      for (l = 0; l < nt; l++) {
        for (m = 0; m < nu; m++) {
          for (j = 0; j < ny; j++)
            w[j] = f[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(y, w, ny, NATSPLN, NATSPLN, s);
          for (j = 0; j < ny; j++)
            sy[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[j];
          for (j = 0; j < ny; j++)
            w[j] = sx[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(y, w, ny, NATSPLN, NATSPLN, s);
          for (j = 0; j < ny; j++)
            sxy[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[j];
        }
      }
    }
  }

  /* 
   *     Cubic spline of:  f vs. z -> sz
   * Cascaded splines of: sy vs. z -> syz
   *                      sx vs. z -> sxz
   *                     sxy vs. z -> sxyz
   */
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (l = 0; l < nt; l++) {
        for (m = 0; m < nu; m++) {
          for (k = 0; k < nz; k++)
            w[k] = f[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
          for (k = 0; k < nz; k++)
            sz[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[k];
          for (k = 0; k < nz; k++)
            w[k] = sy[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
          for (k = 0; k < nz; k++)
            syz[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[k];
          for (k = 0; k < nz; k++)
            w[k] = sx[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
          for (k = 0; k < nz; k++)
            sxz[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[k];
          for (k = 0; k < nz; k++)
            w[k] = sxy[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(z, w, nz, NATSPLN, NATSPLN, s);
          for (k = 0; k < nz; k++)
            sxyz[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[k];
        }
      }
    }
  }

  /*     Cubic spline of:   f vs. t -> st
   *                       sy vs. t -> syt
   *                       sz vs. t -> szt
   *                      sxy vs. t -> sxyt
   *                      syz vs. t -> syzt
   *                      sxz vs. t -> sxzt
   *                     sxyz vs. t -> sxyzt
   */
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        for (m = 0; m < nu; m++) {
          for (l = 0; l < nt; l++)
            w[l] = f[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
          for (l = 0; l < nt; l++)
            st[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[l];
          for (l = 0; l < nt; l++)
            w[l] = sy[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
          for (l = 0; l < nt; l++)
            syt[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[l];
          for (l = 0; l < nt; l++)
            w[l] = sz[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
          for (l = 0; l < nt; l++)
            szt[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[l];
          for (l = 0; l < nt; l++)
            w[l] = sxy[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s); 
          for (l = 0; l < nt; l++)
            sxyt[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[l];
          for (l = 0; l < nt; l++)
            w[l] = syz[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s); 
          for (l = 0; l < nt; l++)
            syzt[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[l];
          for (l = 0; l < nt; l++)
            w[l] = sxz[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
          for (l = 0; l < nt; l++)
            sxzt[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[l];
          for (l = 0; l < nt; l++)
            w[l] = sxyz[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(t, w, nt, NATSPLN, NATSPLN, s);
          for (l = 0; l < nt; l++)
            sxyzt[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[l];
        }
      }
    }
  }

  /*     Cubic spline of:     f vs. u -> su
   * cascaded splines of:    st vs. u -> stu
   *                         sy vs. u -> syu
   *                         sx vs. u -> sxu
   *                        sxy vs. u -> sxyu
   *                         sz vs. u -> szu
   *                        syz vs. u -> sxyu
   *                        sxz vs. u -> sxzu
   *                       sxyz vs. u -> sxyzu
   *                         st vs. u -> stu
   *                        syt vs. u -> sytu
   *                        szt vs. u -> sztu
   *                       sxyt vs. u -> sxytu
   *                       syzt vs. u -> syztu
   *                       sxzt vs. u -> sxztu
   *                      sxyzt vs. u -> sxyztu
   */
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        for (l = 0; l < nt; l++) {
          for (m = 0; m < nu; m++)
            w[m] = f[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            su[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = st[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            stu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = sy[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            syu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = sx[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            sxu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = sxy[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            sxyu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = sz[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            szu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = syz[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            syzu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = sxz[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            sxzu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = sxyz[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            sxyzu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = st[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            stu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = syt[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            sytu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = szt[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            sztu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = sxyt[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            sxytu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = syzt[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            syztu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = sxzt[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            sxztu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
          for (m = 0; m < nu; m++)
            w[m] = sxyzt[(((i*ny+j)*nz+k)*nt+l)*nu+m];
          nerr += csp1d(u, w, nu, NATSPLN, NATSPLN, s);
          for (m = 0; m < nu; m++)
            sxyztu[(((i*ny+j)*nz+k)*nt+l)*nu+m] = s[m];
        }
      }
    }
  }
  free(w);
  free(s);
  if (nerr != 0) {
    csp5d_close(splines);
    splines = NULL;
  }
  return splines;
}

double
csp5d_eval(double p[], double *xv[], int nv[], double f[],
    double **splines, int interval[])
/*
 * Purpose
 *
 *   Interpolate a smooth function of five variables tabulated on a
 *   rectilinear grid. We have a function f(x, y, z, t, u) tabulated at fixed
 *   values of x, y, z, t and u. Requires csp5d_open() to calculate coefficients
 *   stored in splines. There is just one call to csp5d_open() ahead of time
 *   and then multiple calls to csp5d_eval() to evaluate the spline at each
 *   point.
 *
 * Input
 *
 *   double p[5]  The point P(x, y, z, t, u) where we want f interpolated
 *   double xv[5] The array of {x[], y[], z[], t[], u[]} where f is tabulated.
 *   int nv[5]    The array of {nx, ny, nz, nt, nu}, the lengths of the 
 *                {x[], y[], z[], t[], u[]} arrays
 *   double f[]   f(x[i], y[j], z[k], t[l], u[l]) packed into a one-dimensional
 *                array of length nx*ny*nz*nt. The u[l] index varies the
 *                fastest.
 *   double **splines
 *                Spline coefficients returned by csp5d_open()
 *   int interval[5]
 *                Work array of length 5 to keep track of interval. If
 *                there are several parallel calls at the same value of p[],
 *                all should share the same interval[] array.
 *
 * Method
 *   
 *   We start by finding the 5-D box in our grid that contains P.
 *
 *     x[i] < x < x[i+1]
 *     y[j] < y < y[j+1]
 *     z[k] < z < z[k+1]
 *     t[l] < t < t[l+1]
 *     u[m] < t < u[m+1]
 *
 *   Then we interpolate our 1-D splines along the 8 edges of the 5-D box
 *   parallel to the u axis to eliminate the u dependence. This gives us
 *   a 4-D box. Next we interpolate along the 16 edges of the 4-D box parallel
 *   to the t axis to eliminate the t dependence. This gives us a 3-D box.
 *   Next we interpolate along the 8 edges of the 3-D box parallel to the
 *   z axis to find a rectangle, eliminating the z dependence. Next we
 *   interpolate along the 2 edges of the rectangle parallel to the y axis to
 *   get a line segment, eliminating the y dependence. Finally, we interpolate
 *   along the line segment parallel to the x axis to elimate the x dependence
 *   to get a single point. This can be considered the process of repeatedly
 *   slicing perpendicular to one of the axes, reducing the dimension by one
 *   each time.
 *
 *   To get a better understanding, of the underlying process it might be
 *   best to read this code from the bottom up.
 *
 * By
 *
 *   Robert W. Strickland, September 29, 2015
 */
{
  int i32[32];
  double f16[16], sy16[16], sx16[16], sxy16[16], sz16[16], sxz16[16],
         syz16[16], sxyz16[16];
  double st16[16], syt16[16], sxt16[16], sxyt16[16], szt16[16], sxzt16[16],
         syzt16[16], sxyzt16[16];
  double f8[8], sy8[8], sx8[8], sxy8[8], sz8[8], sxz8[8], syz8[8], sxyz8[8];
  double f4[4], sy4[4], sx4[4], sxy4[4];
  double f0, f1, s0x, s1x;

  double g, a, b, c, d;
  int i0, i1, j0, j1, k0, k1, l0, l1, m0, m1;
  double dx, dy, dz, dt, du;

  double *sx;
  double *sy, *sxy;
  double *sz, *sxz, *syz, *sxyz;
  double *st, *syt, *sxyt, *szt, *sxzt, *syzt, *sxyzt;
  double *su, *sxu, *syu, *sxyu, *szu, *sxzu, *syzu, *sxyzu,
         *stu, *sytu, *sxytu, *sztu, *sxztu, *syztu, *sxyztu;

  int n, n0, n1;
  int ny, nz, nt, nu;
  double *x, *y, *z, *t, *u;

  /* Linear search for interval containing each independent variable */

  for (n = 0; n < 5; n++) {
    if (interval[n] < 0 || interval[n] > nv[n]-2)
      interval[n] = nv[n]/2;
    while (interval[n] > 0 && p[n] < xv[n][interval[n]])
      interval[n]--;
    while (interval[n] < nv[n]-2 && p[n] > xv[n][interval[n]+1])
      interval[n]++;
  }

  /* Get indices for of corners of 5-D box containing (xx, yy, zz, tt, uu) */

  i0 = interval[0];
  i1 = i0+1;
  j0 = interval[1];
  j1 = j0+1;
  k0 = interval[2];
  k1 = k0+1;
  l0 = interval[3];
  l1 = l0+1;
  m0 = interval[4];
  m1 = m0+1;

  ny = nv[1];
  nz = nv[2];
  nt = nv[3];
  nu = nv[4];

  i32[ 0] = (((i0*ny+j0)*nz+k0)*nt+l0)*nu+m0;
  i32[ 1] = (((i0*ny+j0)*nz+k0)*nt+l0)*nu+m1;
  i32[ 2] = (((i0*ny+j0)*nz+k0)*nt+l1)*nu+m0;
  i32[ 3] = (((i0*ny+j0)*nz+k0)*nt+l1)*nu+m1;
  i32[ 4] = (((i0*ny+j0)*nz+k1)*nt+l0)*nu+m0;
  i32[ 5] = (((i0*ny+j0)*nz+k1)*nt+l0)*nu+m1;
  i32[ 6] = (((i0*ny+j0)*nz+k1)*nt+l1)*nu+m0;
  i32[ 7] = (((i0*ny+j0)*nz+k1)*nt+l1)*nu+m1;
  i32[ 8] = (((i0*ny+j1)*nz+k0)*nt+l0)*nu+m0;
  i32[ 9] = (((i0*ny+j1)*nz+k0)*nt+l0)*nu+m1;
  i32[10] = (((i0*ny+j1)*nz+k0)*nt+l1)*nu+m0;
  i32[11] = (((i0*ny+j1)*nz+k0)*nt+l1)*nu+m1;
  i32[12] = (((i0*ny+j1)*nz+k1)*nt+l0)*nu+m0;
  i32[13] = (((i0*ny+j1)*nz+k1)*nt+l0)*nu+m1;
  i32[14] = (((i0*ny+j1)*nz+k1)*nt+l1)*nu+m0;
  i32[15] = (((i0*ny+j1)*nz+k1)*nt+l1)*nu+m1;
  i32[16] = (((i1*ny+j0)*nz+k0)*nt+l0)*nu+m0;
  i32[17] = (((i1*ny+j0)*nz+k0)*nt+l0)*nu+m1;
  i32[18] = (((i1*ny+j0)*nz+k0)*nt+l1)*nu+m0;
  i32[19] = (((i1*ny+j0)*nz+k0)*nt+l1)*nu+m1;
  i32[20] = (((i1*ny+j0)*nz+k1)*nt+l0)*nu+m0;
  i32[21] = (((i1*ny+j0)*nz+k1)*nt+l0)*nu+m1;
  i32[22] = (((i1*ny+j0)*nz+k1)*nt+l1)*nu+m0;
  i32[23] = (((i1*ny+j0)*nz+k1)*nt+l1)*nu+m1;
  i32[24] = (((i1*ny+j1)*nz+k0)*nt+l0)*nu+m0;
  i32[25] = (((i1*ny+j1)*nz+k0)*nt+l0)*nu+m1;
  i32[26] = (((i1*ny+j1)*nz+k0)*nt+l1)*nu+m0;
  i32[27] = (((i1*ny+j1)*nz+k0)*nt+l1)*nu+m1;
  i32[28] = (((i1*ny+j1)*nz+k1)*nt+l0)*nu+m0;
  i32[29] = (((i1*ny+j1)*nz+k1)*nt+l0)*nu+m1;
  i32[30] = (((i1*ny+j1)*nz+k1)*nt+l1)*nu+m0;
  i32[31] = (((i1*ny+j1)*nz+k1)*nt+l1)*nu+m1;

  /* Retrieve the spline pointers from splines array. The order here must
   * match that in csp5d_open()
   */
  sx = splines[0];
  sy = splines[1];
  sxy = splines[2];
  sz = splines[3];
  syz = splines[4];
  sxz = splines[5];
  sxyz = splines[6];
  st = splines[7];
  syt = splines[8];
  szt = splines[9];
  sxyt = splines[10];
  syzt = splines[11];
  sxzt = splines[12];
  sxyzt = splines[13];
  su = splines[14];
  stu = splines[15];
  sxu = splines[16];
  syu = splines[17];
  sytu = splines[18];
  szu = splines[19];
  sztu = splines[20];
  sxyu = splines[21];
  sxytu = splines[22];
  syzu = splines[23];
  syztu = splines[24];
  sxzu = splines[25];
  sxztu = splines[26];
  sxyzu = splines[27];
  sxyztu = splines[28];

  /* Interpolate in the u direction to get all 16 corners of the 4-D box.
   *
   * Our point P=(xx, yy, zz, tt, uu) lies in the 5-D box
   *
   *   x[i] < xx < x[i+1]
   *   y[j] < yy < y[j+1]
   *   z[k] < zz < z[k+1]
   *   t[l] < tt < t[l+1]
   *   u[m] < uu < u[m+1]
   *
   * We can index the 32 corners of this box with a 5-bit binary index
   * where the most significant bit refers to the x coordinate and the
   * least significant bit refers to the u coordinate. For instance, index
   * 01101 refers to the point (x[i],y[j+1],z[k+1],t[l],u[m+1]).
   *
   * First we'll interpolate in the u direction. We'll interpolate all
   * of the values of f and all of the splines we'll need for later steps.
   * After we interpolate for u, we will have values on the corners of
   * a 4-D box. We can drop the least significant bit (the u index) and
   * use a 4-bit index.
   */

  u = xv[4];
  du = u[m1]-u[m0];
  if (du <= 0.)
    return 1.e30;
  b = (p[4] - u[m0])/du;
  a = 1.-b;
  g = du*du/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  for (n = 0; n < 16; n++) {
    n0 = i32[2*n];
    n1 = i32[2*n+1];
    f16    [n] = a*    f[n0] + b*    f[n1] + c*    su[n0] + d*    su[n1];
    st16   [n] = a*   st[n0] + b*   st[n1] + c*   stu[n0] + d*   stu[n1];
    sx16   [n] = a*   sx[n0] + b*   sx[n1] + c*   sxu[n0] + d*   sxu[n1];
    sy16   [n] = a*   sy[n0] + b*   sy[n1] + c*   syu[n0] + d*   syu[n1];
    syt16  [n] = a*  syt[n0] + b*  syt[n1] + c*  sytu[n0] + d*  sytu[n1];
    sz16   [n] = a*   sz[n0] + b*   sz[n1] + c*   szu[n0] + d*   szu[n1];
    szt16  [n] = a*  szt[n0] + b*  szt[n1] + c*  sztu[n0] + d*  sztu[n1];
    sxy16  [n] = a*  sxy[n0] + b*  sxy[n1] + c*  sxyu[n0] + d*  sxyu[n1];
    sxyt16 [n] = a* sxyt[n0] + b* sxyt[n1] + c* sxytu[n0] + d* sxytu[n1];
    syz16  [n] = a*  syz[n0] + b*  syz[n1] + c*  syzu[n0] + d*  syzu[n1];
    syzt16 [n] = a* syzt[n0] + b* syzt[n1] + c* syztu[n0] + d* syztu[n1];
    sxz16  [n] = a*  sxz[n0] + b*  sxz[n1] + c*  sxzu[n0] + d*  sxzu[n1];
    sxzt16 [n] = a* sxzt[n0] + b* sxzt[n1] + c* sxztu[n0] + d* sxztu[n1];
    sxyz16 [n] = a* sxyz[n0] + b* sxyz[n1] + c* sxyzu[n0] + d* sxyzu[n1];
    sxyzt16[n] = a*sxyzt[n0] + b*sxyzt[n1] + c*sxyztu[n0] + d*sxyztu[n1];
  }

  /* Interpolate in the t direction to get all 8 corners of the 3-D box
   *
   * We can index the 16 corners of this box with a 4-bit binary index
   * where the most significant bit refers to the x coordinate and the
   * least significant bit refers to the t coordinate. For instance, index
   * 0110 refers to the point (x[i],y[j+1],z[k+1],t[l]).
   *
   * Next we'll interpolate in the t direction. We'll interpolate all
   * of the values of f and all of the splines we'll need for later steps.
   * After we interpolate for t, we will have values on the corners of
   * a 3-D box. We can drop the least significant bit (the t index) and
   * use a 3-bit index.
   *
   *     011+--------+------------+111
   *        |\   z   |\           |\
   *        | \  ^   | \          | \
   *        |  +-|------+------------+
   *        |  |\|   |  |\        |  |\
   *        | 001+--------+------------+101
   *        +--|-|---+--|-|-------+  | |
   *        |\.|.|....\.|.|........\ | |
   *        | \|.|.....\|.|.........\| |
   *        |  U-|------P-|----------V |
   *        |  |\|.......\|...........\|
   *      y |  | +--------+------------+
   *       \|  | |   |  | |       |  | |
   *     010+--|-|---+--|-|-------+110 |
   *         \ | |    \ | |        \ | |
   *          \| |     \| |         \| |
   *           +-|------+-|----------+ |
   *            \|       \|           \|
   *          000+--------+------------+100->x
   *
   * Again, the msb refers to the x coordinate.
   */

  t = xv[3];
  dt = t[l1]-t[l0];
  if (dt <= 0.)
    return 1.e30;
  b = (p[3] - t[l0])/dt;
  a = 1.-b;
  g = dt*dt/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  /* The n index goes through all 8 corners of the box from binary 000 to 111
   * (decimal 7) in binary counting order.
   *
   * n0 is the index into the f[] array (or one of its splines) with
   * tt=t[l]. It is the 4-bit index to a corner of the the 4-D box
   * containing our point.
   *
   * n1 is the index into the f[] array (or one of its splines) with
   * tt=t[l+1]. It, too is a 4-bit corner index.
   */
  for (n = 0; n < 8; n++) {
    n0 = 2*n;
    n1 = 2*n+1;
    f8   [n] = a*   f16[n0] + b*   f16[n1] + c*   st16[n0] + d*   st16[n1];
    sx8  [n] = a*  sx16[n0] + b*  sx16[n1] + c*  sxt16[n0] + d*  sxt16[n1];
    sy8  [n] = a*  sy16[n0] + b*  sy16[n1] + c*  syt16[n0] + d*  syt16[n1];
    sz8  [n] = a*  sz16[n0] + b*  sz16[n1] + c*  szt16[n0] + d*  szt16[n1];
    sxy8 [n] = a* sxy16[n0] + b* sxy16[n1] + c* sxyt16[n0] + d* sxyt16[n1];
    syz8 [n] = a* syz16[n0] + b* syz16[n1] + c* syzt16[n0] + d* syzt16[n1];
    sxz8 [n] = a* sxz16[n0] + b* sxz16[n1] + c* sxzt16[n0] + d* sxzt16[n1];
    sxyz8[n] = a*sxyz16[n0] + b*sxyz16[n1] + c*sxyzt16[n0] + d*sxyzt16[n1];
  }

  /* Interpolate in the z direction along the edges of the 3-box to get
   * corners of the rectangle parallel to the x-z plane at zz (the shaded
   * plane above)
   *
   *  y[j+1] -01-----+---11
   *           |     |    |
   *      yy --U-----P----V
   *           |     |    |
   *           |     |    |
   *    y[j] -00-----+---10
   *           |     |    |
   *         x[i]   xx  x[i+1]
   */

  z = xv[2];
  dz = z[k1]-z[k0];
  if (dz <= 0.)
    return 1.e30;
  b = (p[2] - z[k0])/dz;
  a = 1.-b;
  g = dz*dz/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  /* The index n refers to the corners of the rectangle. Here we have dropped
   * the z index, and n takes on all the values from binary 00 to 11. The
   * msb still refers to x.
   *
   * n0 refers to the bottom z[k] plane of our 3-D box, and
   * n1 refers to the top z[k+1] plane.
   */
  for (n = 0; n < 4; n++) {
    n0 = 2*n;
    n1 = 2*n+1;
    f4  [n] = a*  f8[n0] + b*  f8[n1] + c*  sz8[n0] + d*  sz8[n1];
    sx4 [n] = a* sx8[n0] + b* sx8[n1] + c* sxz8[n0] + d* sxz8[n1];
    sy4 [n] = a* sy8[n0] + b* sy8[n1] + c* syz8[n0] + d* syz8[n1];
    sxy4[n] = a*sxy8[n0] + b*sxy8[n1] + c*sxyz8[n0] + d*sxyz8[n1];
  }

  /* Interpolate in the y direction across the rectangle to get function
   * and spline at the endpoints of our line segment from points U to V
   * (here indexed 0 and 1)
   */
  y = xv[1];
  dy = y[j1]-y[j0];
  if (dy <= 0.)
    return 1.e30;
  b = (p[1] - y[j0])/dy;
  a = 1.-b;
  g = dy*dy/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  f0  = a* f4[0] + b* f4[1] + c* sy4[0] + d* sy4[1];
  f1  = a* f4[2] + b* f4[3] + c* sy4[2] + d* sy4[3];
  s0x = a*sx4[0] + b*sx4[1] + c*sxy4[0] + d*sxy4[1];
  s1x = a*sx4[2] + b*sx4[3] + c*sxy4[2] + d*sxy4[3];

  /* Interpolate in the x direction along the line segment to get f(P) */

  x = xv[0];
  dx = x[i1]-x[i0];
  if (dx <= 0.)
    return 1.e30;
  b = (p[0] - x[i0])/dx;
  a = 1.-b;
  g = dx*dx/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  return a*f0 + b*f1 + c*s0x + d*s1x;
}

#include <math.h>
#include <time.h>

int testmain5d(void)
{
  double *f;
  int nx, ny, nz, nt, nu, nf;
  int i, j, k, l, m;
  double xx, yy, zz, tt, uu;
  double *x, *y, *z, *t, *u;
  double r;
  double **splines;
  int interval[5];
  double *xv[5];
  int nv[5];
  double p[5];
  double elapsed;
  clock_t tm;

  nx = 12;
  ny = 12;
  nz = 12;
  nt = 12;
  nu = 12;

  x = (double *) malloc(nx*sizeof(double));
  y = (double *) malloc(ny*sizeof(double));
  z = (double *) malloc(nz*sizeof(double));
  t = (double *) malloc(nt*sizeof(double));
  u = (double *) malloc(nu*sizeof(double));

  xv[0] = x;
  xv[1] = y;
  xv[2] = z;
  xv[3] = t;
  xv[4] = u;

  x[0] = 0.;
  y[0] = 0.;
  z[0] = 0.;
  t[0] = 0.;
  u[0] = 0.;
  /*
   * To get proper first derivative at center of splash, use a point close
   * to edge
   */
  x[1] = 0.1;
  y[1] = 0.1;
  z[1] = 0.1;
  t[1] = 0.1;
  u[1] = 0.1;
  for (i = 2; i < nx; i++) {
    x[i] = i-1.;
    y[i] = i-1.;
    z[i] = i-1.;
    t[i] = i-1.;
    u[i] = i-1.;
  }

  nf = nx*ny*nz*nt*nu;
  nv[0] = nx;
  nv[1] = ny;
  nv[2] = nz;
  nv[3] = nt;
  nv[4] = nu;

  f =  (double *) malloc(nf*sizeof(double));

  nf = 0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        for (l = 0; l < nt; l++) {
          for (m = 0; m < nu; m++) {
            r = x[i]*x[i] + y[j]*y[j] + z[k]*z[k] + t[l]*t[l] + u[m]*u[m];
            r = M_PI*sqrt(r)/2.;
            yy = (r != 0.) ? sin(r)/r : 1.;
            int q = (((i*ny+j)*nz+k)*nt+l)*nu+m;
            if (q != nf) {
              fprintf(stderr, "Mismatch\n");
              return 1;
            }
            f[q] = yy;
            nf++;
          }
        }
      }
    }
  }
  fprintf(stderr, "Function values calculated\n");
  splines = csp5d_open(xv, nv, f);
  if (splines == NULL) {
    fprintf(stderr, "failed\n");
    return 1;
  }
  fprintf(stderr, "Splines calculated\n");
  for (i = 0; i < 5; i++)
    interval[i] = 0;
  tm = clock();
  for (i = 0; i < 5000; i++) {
    double e, spl;
    double ex[5];

    if (i == 0) {
      ex[0] = 1.;
      ex[1] = 0.;
      ex[2] = 0.;
      ex[3] = 0.;
      ex[4] = 0.;
    }
    else if (i == 1) {
      ex[0] = 1.;
      ex[1] = 1.;
      ex[2] = 1.;
      ex[3] = 1.;
      ex[4] = 1.;
    }
    else {
      for (j = 0; j < 5; j++) {
        ex[j] = rand()/(double) RAND_MAX;
      }
      ex[rand() % 5] = 1.; /* force to surface of hyperbox */
    }
    for (j = 0; j < 201; j++) {
      xx = j * ex[0] / 20.;
      yy = j * ex[1] / 20.;
      zz = j * ex[2] / 20.;
      tt = j * ex[3] / 20.;
      uu = j * ex[4] / 20.;
      p[0] = xx;
      p[1] = yy;
      p[2] = zz;
      p[3] = tt;
      p[4] = uu;
      spl = csp5d_eval(p, xv, nv, f, splines, interval);
      r = M_PI*sqrt(xx*xx + yy*yy + zz*zz + tt*tt + uu*uu)/2.;
      e = (r != 0.) ? sin(r)/r : 1.;
      printf("%8.5f %9.6f %9.6f\n", r, spl, e);
    }
  }
  tm = clock()-tm;
  free(f);
  csp5d_close(splines);
  elapsed = (double) tm / (double) CLOCKS_PER_SEC;
  fprintf(stderr, "5000 points in %.2f seconds = %.2f points per second\n",
      elapsed, 5000./elapsed);
  return 0;
}

int testmain4d(void)
{
  double *f;
  int nx, ny, nz, nt, nf;
  int i, j, k, l;
  double xx, yy, zz, tt;
  double *x, *y, *z, *t;
  double r;
  double **splines;
  int interval[4];
  double *xv[4];
  int nv[4];
  double p[4];

  nx = 12;
  ny = 12;
  nz = 12;
  nt = 12;

  x = (double *) malloc(nx*sizeof(double));
  y = (double *) malloc(ny*sizeof(double));
  z = (double *) malloc(nz*sizeof(double));
  t = (double *) malloc(nt*sizeof(double));

  xv[0] = x;
  xv[1] = y;
  xv[2] = z;
  xv[3] = t;

  x[0] = 0.;
  y[0] = 0.;
  z[0] = 0.;
  t[0] = 0.;
  /*
   * To get proper first derivative at center of splash, use a point close
   * to edge
   */
  x[1] = 0.1;
  y[1] = 0.1;
  z[1] = 0.1;
  t[1] = 0.1;
  for (i = 2; i < nx; i++) {
    x[i] = i-1.;
    y[i] = i-1.;
    z[i] = i-1.;
    t[i] = i-1.;
  }

  nf = nx*ny*nz*nt;
  nv[0] = nx;
  nv[1] = ny;
  nv[2] = nz;
  nv[3] = nt;

  f =  (double *) malloc(nf*sizeof(double));

  nf = 0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        for (l = 0; l < nt; l++) {
          r = M_PI*sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k] + t[l]*t[l])/2.;
          yy = (r != 0.) ? sin(r)/r : 1.;
          int m = ((i*ny+j)*nz+k)*nt+l;
          if (m != nf) {
            printf("Mismatch\n");
            return 1;
          }
          f[m] = yy;
          nf++;
        }
      }
    }
  }
  splines = csp4d_open(xv, nv, f);
  for (i = 0; i < 4; i++)
    interval[i] = 0;
  clock_t tm = clock();
  for (i = 0; i < 10000; i++) {
    double e, spl;
    double ex[4];

    if (i == 0) {
      ex[0] = 1.;
      ex[1] = 0.;
      ex[2] = 0.;
      ex[3] = 0.;
    }
    else if (i == 1) {
      ex[0] = 1.;
      ex[1] = 1.;
      ex[2] = 1.;
      ex[3] = 1.;
    }
    else {
      for (j = 0; j < 4; j++) {
        ex[j] = rand()/(double) RAND_MAX;
      }
      ex[rand() % 4] = 1.; /* force to surface of hyperbox */
    }
    for (j = 0; j < 201; j++) {
      xx = j * ex[0] / 20.;
      yy = j * ex[1] / 20.;
      zz = j * ex[2] / 20.;
      tt = j * ex[3] / 20.;
      p[0] = xx;
      p[1] = yy;
      p[2] = zz;
      p[3] = tt;
      spl = csp4d_eval(p, xv, nv, f, splines, interval);
      r = M_PI*sqrt(xx*xx + yy*yy + zz*zz + tt*tt)/2.;
      e = (r != 0.) ? sin(r)/r : 1.;
      printf("%8.5f %9.6f %9.6f\n", r, spl, e);
    }
  }
  tm = clock()-tm;
  double elapsed = (double) tm / (double) CLOCKS_PER_SEC;
  fprintf(stderr, "10000 points in %.2f seconds = %.2f points per second\n",
      elapsed, 10000./elapsed);
  free(f);
  csp4d_close(splines);
  return 0;
}

int testmain3d(void)
{
  double *f;
  int nx, ny, nz, nf;
  int i, j, k;
  double xx, yy, zz;
  double *x, *y, *z;
  double r;
  double **splines;
  int interval[3];
  double *xv[3];
  int nv[3];
  double p[3];

  nx = 12;
  ny = 12;
  nz = 12;

  x = (double *) malloc(nx*sizeof(double));
  y = (double *) malloc(ny*sizeof(double));
  z = (double *) malloc(nz*sizeof(double));

  xv[0] = x;
  xv[1] = y;
  xv[2] = z;

  x[0] = 0.;
  y[0] = 0.;
  z[0] = 0.;
  /*
   * To get proper first derivative at center of splash, use a point close
   * to edge
   */
  x[1] = 0.1;
  y[1] = 0.1;
  z[1] = 0.1;
  for (i = 2; i < nx; i++) {
    x[i] = i-1.;
    y[i] = i-1.;
    z[i] = i-1.;
  }

  nf = nx*ny*nz;
  nv[0] = nx;
  nv[1] = ny;
  nv[2] = nz;

  f =  (double *) malloc(nf*sizeof(double));

  nf = 0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        r = M_PI*sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])/2.;
        yy = (r != 0.) ? sin(r)/r : 1.;
        int m = (i*ny+j)*nz+k;
        if (m != nf) {
          printf("Mismatch\n");
          return 1;
        }
        f[m] = yy;
        nf++;
      }
    }
  }
  splines = csp3d_open(xv, nv, f);
  for (i = 0; i < 3; i++)
    interval[i] = 0;
  for (i = 0; i < 1002; i++) {
    double e, spl;
    double ex[3];

    if (i == 0) {
      ex[0] = 1.;
      ex[1] = 0.;
      ex[2] = 0.;
    }
    else if (i == 1) {
      ex[0] = 1.;
      ex[1] = 1.;
      ex[2] = 1.;
    }
    else {
      for (j = 0; j < 3; j++) {
        ex[j] = rand()/(double) RAND_MAX;
      }
      ex[rand() % 3] = 1.; /* force to surface of hyperbox */
    }
    for (j = 0; j < 101; j++) {
      xx = j * ex[0] / 10.;
      yy = j * ex[1] / 10.;
      zz = j * ex[2] / 10.;
      p[0] = xx;
      p[1] = yy;
      p[2] = zz;
      spl = csp3d_eval(p, xv, nv, f, splines, interval);
      r = M_PI*sqrt(xx*xx + yy*yy + zz*zz)/2.;
      e = (r != 0.) ? sin(r)/r : 1.;
      printf("%8.5f %9.6f %9.6f\n", r, spl, e);
    }
  }
  free(f);
  csp3d_close(splines);
  return 0;
}

int testmain2d(void)
{
  double *f;
  int nx, ny, nf;
  int i, j;
  double xx, yy;
  double *x, *y;
  double r;
  double **splines;
  int interval[2];
  double *xv[2];
  int nv[2];
  double p[2];

  nx = 12;
  ny = 12;

  x = (double *) malloc(nx*sizeof(double));
  y = (double *) malloc(ny*sizeof(double));

  xv[0] = x;
  xv[1] = y;

  x[0] = 0.;
  y[0] = 0.;
  /*
   * To get proper first derivative at center of splash, use a point close
   * to edge
   */
  x[1] = 0.1;
  y[1] = 0.1;
  for (i = 2; i < nx; i++) {
    x[i] = i-1.;
    y[i] = i-1.;
  }

  nf = nx*ny;
  nv[0] = nx;
  nv[1] = ny;

  f =  (double *) malloc(nf*sizeof(double));

  nf = 0;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      r = M_PI*sqrt(x[i]*x[i] + y[j]*y[j])/2.;
      yy = (r != 0.) ? sin(r)/r : 1.;
      int m = i*ny+j;
      if (m != nf) {
        printf("Mismatch\n");
        return 1;
      }
      f[m] = yy;
      nf++;
    }
  }
  splines = csp2d_open(xv, nv, f);
  for (i = 0; i < 2; i++)
    interval[i] = 0;
  for (i = 0; i < 50; i++) {
    double e, spl;
    double ex[2];

    if (i == 0) {
      ex[0] = 1.;
      ex[1] = 0.;
    }
    else if (i == 1) {
      ex[0] = 1.;
      ex[1] = 1.;
    }
    else {
      for (j = 0; j < 2; j++) {
        ex[j] = rand()/(double) RAND_MAX;
      }
      ex[rand() % 2] = 1.; /* force to surface of hyperbox */
    }
    for (j = 0; j < 201; j++) {
      xx = j * ex[0] / 20.;
      yy = j * ex[1] / 20.;
      p[0] = xx;
      p[1] = yy;
      spl = csp2d_eval(p, xv, nv, f, splines, interval);
      r = M_PI*sqrt(xx*xx + yy*yy)/2.;
      e = (r != 0.) ? sin(r)/r : 1.;
      printf("%8.5f %9.6f %9.6f %8.5f %8.5f\n", r, spl, e, xx, yy);
    }
  }
  free(f);
  csp2d_close(splines);
  return 0;
}

int testmain1d(void)
{
  double *x, *f, *s, xx, yy, r, spl, e;
  int i, j, k, n;
  double dx, a, b, c, d, g;

  n = 12;

  x = (double *) malloc(n*sizeof(double));
  f = (double *) malloc(n*sizeof(double));
  s = (double *) malloc(n*sizeof(double));

  x[0] = 0.;
  /*
   * To get proper first derivative at center of splash, use a point close
   * to edge
   */
  x[1] = 0.1;
  for (i = 2; i < n; i++) {
    x[i] = i-1.;
  }

  for (i = 0; i < n; i++) {
    r = M_PI*sqrt(x[i]*x[i])/2.;
    yy = (r != 0.) ? sin(r)/r : 1.;
    f[i] = yy;
    printf("%8.5f %9.6f\n", r, yy);
  }
  if (csp1d(x, f, n, NATSPLN, NATSPLN, s) != 0)
    return 1;
  k = 0;
  for (j = 0; j < 101; j++) {
    xx = j / 10.;
    while (k > 0 && xx < x[k])
      k--;
    while (k < n-2 && xx > x[k+1])
      k++;
    dx = x[k+1]-x[k];
    b = (xx - x[k])/dx;
    a = 1.-b;
    g = dx*dx/6.;
    c = g*(a*a*a - a);
    d = g*(b*b*b - b);
    spl = a*f[k] + b*f[k+1] + c*s[k] + d*s[k+1];
    r = M_PI*sqrt(xx*xx)/2.;
    e = (r != 0.) ? sin(r)/r : 1.;
    printf("%8.5f %9.6f %9.6f\n", r, spl, e);
  }
  free(f);
  free(s);
  free(x);
  return 0;
}

/*
int
main(void)
{
  return testmain5d();
}
*/
