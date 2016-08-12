/**
 * Main.cc
 *
 *  @date Created on: Jun 20, 2016
 *  @author Keith D Brauss
 */

/**
 * @mainpage
 *
 * @section section_toc Table of Contents
 * <ul>
 *   <li> @ref introduction_section </li>
 *   <li> @ref complex_analysis_section </li>
 *     <ul>
 *       <li> @ref complex_plane_subsection </li>
 *       <li> @ref complex_exponential_function_subsection </li>
 *       <li> @ref complex_logarithmic_function_subsection </li>
 *     </ul>
 *   <li> @ref fmm_2d_section </li>
 *     <ul>
 *       <li> @ref approximating_the_potential_subsection </li>
 *       <li> @ref field_expansion_subsection </li>
 *         <ul>
 *           <li> @ref far_field_expansion_subsubsection </li>
 *           <li> @ref near_field_expansion_subsubsection </li>
 *         </ul>
 *       <li> @ref translation_subsection </li>
 *         <ul>
 *           <li> @ref local_to_local_translation_subsubsection </li>
 *           <li> @ref far_to_near_translation_subsubsection </li>
 *           <li> @ref far_to_far_translation_subsubsection </li>
 *         </ul>
 *       <li> @ref algorithm_subsection </li>
 *         <ul>
 *           <li> @ref quadtree_partition_subsubsection </li>
 *           <li> @ref far_field_expansion_subsubsection </li>
 *           <li> @ref upward_pass_subsubsection </li>
 *           <li> @ref downward_pass_subsubsection </li>
 *         </ul>
 *     </ul>
 * </ul>
 *
 * @section introduction_section Introduction
 *
 * This code implements the Fast Multipole Method (FMM) in two dimensions (2D)
 * for the calculation of the elecctric potential experienced by a particle \f$y_j\f$ from
 * the presence of charged particles \f$x_i\f$.  The potential function describing
 * the action (from a distance) of \f$x_i\f$ on \f$y_j\f$ is given by the formula
 * \f[\phi_{ji} = q_i \ln ( \|y_j - x_i\| )\f]
 * where \f$x_i\f$ and \f$y_j\f$ are the coordinates for particles in two-dimensional space
 * and \f$q_i\f$ is the charge of particle \f$x_i\f$.
 * <ul>
 *   <li> \f$x_i\f$ and \f$y_j\f$ are complex numbers . </li>
 *   <li> \f$y_j\f$ are target particles (points). </li>
 *   <li> \f$x_i\f$ are source particles </li>
 * </ul>
 *
 *  For a target particle \f$y_j\f$ in the complex plane that encounters action from \f$N_s\f$ multiple source
 *  particles \f$x_i\f$ located in the complex plane (we can think of the xy-plane \f$\mathbb{R}^2\f$ since there is a
 *  one-to-one correspondece between the complex plane and the xy-plane), the potential is calculated
 *  using superposition.
 *  \f[ \sum_{i=1}^{N_s} \phi_{ji} = \sum_{i=1}^{N_s} \ln (\|y_j - x_i\|) q_i
 *                                 = \begin{bmatrix} \ln(\|y_j - x_1\|) & \ln(\|y_j - x_2\|) & \hdots & \ln(\|y_j-x_{N_s}\|)
 *                                   \end{bmatrix}
 *                                   \begin{bmatrix} q_1 \\ q_2 \\ \vdots \\ q_{N_s} \end{bmatrix} \f]
 *
 *  <ol type="1">
 *    <li> The summation above can be viewed as a vector-vector product having an operation count \f$O(N_s)\f$. </li>
 *    <li> We are interested in performing this operation for a number \f$N_t\f$ of target particles \f$y_j\f$. </li>
 *    <li> Performing the \f$N_t\f$ summation operations can be viewed as a matrix-vector product. </li>
 *  </ol>
 *  The goal of the code and FMM is to minimize the \f$O(N^2)\f$ operation count (assuming \f$N = N_t = N_s\f$)
 *  required to perform the matrix-vector product.  To discuss the FMM for two dimensions and for the electrostatic
 *  potential function, we need to review some complex numbers and complex analysis.
 *
 *  @section complex_analysis_section Complex Analysis Review
 *
 *  For a more in-depth review of complex analysis see Chapter 13 of \cite kreyszig2010.
 *  Recall the complex plane above.
 *
 *  @subsection complex_plane_subsection Complex Plane
 *
 *
 *
 *  @image html polar.png
 *
 *  Consider the complex plane shown above.
 *  <ul>
 *    <li> the x-axis corresponds to the real axis </li>
 *    <li> the y-axis is the imaginary axis </li>
 *    <li> the complex number \f$z = x + iy\f$ is the point in the complex plane with coordinates \f$(x,y)\f$
 *    <li> \f$\theta\f$ is the angle from the positive x-axis to the complex number (point)
 *         \f$z = x + iy\f$ </li>
 *    <li> \f$\theta\f$ is called the argument of z, denoted \f$arg(z)\f$ (\f$\theta = arg(z)\f$) </li>
 *    <li> from the Pythagorean Theorem, the hypotenuse \f$r = \|z\| = \sqrt{x^2 + y^2}\f$.
 *    <li> using trig, we have \f$\tan (\theta) = \frac{y}{x}\f$ </li>
 *    <li> therefore \f$z = x + iy = r \cos (\theta) + i r \sin(\theta) = r( \cos(\theta) + i \sin (\theta) \f$ </li>
 *    <li> using addition formulas \f$z_1z_2 = r_1r_2(cos(\theta_1+\theta_2)+i\sin(\theta_1+\theta_2))\f$ where
 *              \f$z_i = r_i(\cos(\theta_i)+i\sin(\theta_i))\f$ </li>
 *  </ul>
 *
 *  @subsection complex_exponential_function_subsection Complex Exponential Function
 *
 *   Let \f$f(z) \colon \mathbb{C} \to \mathbb{C}\f$ be the complex-valued function
 *   defined by w = f(z) = u(x,y) + iv(x,y).  From the Cauchy-Riemann equations, we have
 *   the following formulas for the complex derivative of \f$f\f$
 *            \f[ f'(z) = \lim_{\Delta z \to 0} \frac{f(z+\Delta z) - f(z)}{\Delta z}
 *                         = u_x + iv_x = -iu_y + v_y \f]
 *
 *  To retain the properties enjoyed by the real-value exponential function, we define
 *  the complex exponential function \f$f(z) = e^z\f$ for \f$z = x + iy\f$ by the formula
 *   \f[e^z = e^x (cos(y) + isin(y))\f]
 *  The following properties result
 *  <ul>
 *    <li>  \f$(e^z)' = e^xcos(y))_x + i(e^xsin(y))_x = e^xcos(y) + ie^xsin(y) = e^z\f$ </li>
 *    <li>  \f$e^{z_1 + z_2} = e^{z_1}e^{z_2}\f$ </li>
 *    <li>  \f$e^z = e^x (cos(0)+isin(0) = e^x\f$ if \f$z = x \in \mathbb{R}\f$ </li>
 *    <li>  \f$e^z = e^0(\cos(y)+i\sin(y)) = e^{iy} = \cos(y)+i\sin(y)\f$ if \f$z = iy\f$ is purely imaginary </li>
 *  </ul>
 *
 *  Some important properties of \f$z = x + iy\f$ and \f$e^z\f$ are
 *  <ul>
 *    <li> \f$z = r(\cos(\theta) + i\sin(\theta)) = re^{i \theta} \f$  (so \f$y = \theta = arg(z)\f$) </li>
 *    <li> \f$e^{2 \pi i} = cos(2 \pi) + i sin(2 \pi) = 1\f$ </li>
 *    <li> \f$e^z = e^xe^{iy} = e^xe^{iy}e^{2 \pi i} = e^{z + 2 \pi i} = e^x e^(y + 2 \pi i)\f$ </li>
 *    <li> \f$arg(e^z) = y +/- 2 n \pi\f$ </li>
 *    <li> \f$\|e^{iy}\| = \|cos(y) + isin(y)\| = cos^2(y) + sin^2(y) = 1\f$ </li>
 *    <li> \f$\|e^z\| = \|e^xe^iy\| = \|e^x\| \|e^iy\| = \|e^x\| = e^x\f$ </li>
 *  </ul>
 *
 *  @subsection complex_logarithmic_function_subsection Complex Logarithmic Function
 *
 *  We define the complex logarithmic function \f$w = \ln(z)\f$ for \f$z = x + iy\f$
 *  as the inverse of the complex exponential function
 *  \f[ f(z) = e^z \quad \mbox{and} \quad  f^{-1}(z) = \ln(z) \f]
 *  by way of the formula
 *  \f[ w = \ln(z)  \quad \iff \quad e^w = z \f]
 *
 *  Let the complex numbers \f$w = u + iv\f$ and \f$z = re^{i \theta}\f$.
 *  Recall that \f$r = \|z\|\f$ and \f$arg(z) = \theta\f$.  Then
 *  <ul>
 *    <li> \f$z = e^w \implies e^{u + iv} = re^{i \theta}\f$ and \f$e^u = r\f$ and \f$v = \theta\f$ </li>
 *    <li> since \f$e^u = r\f$ and \f$r \in \mathbb{R}\f$, \f$u = \ln(r)\f$ (where \f$\ln\f$ is the real natural logarithmic function) </li>
 *    <li> therefore \f$ln(z) = w = u + iv = ln(r) + i \theta\f$ </li>
 *    <li> recall that \f$arg(z) = \theta\f$ is determined only up to integer multiples of \f$2 \pi\f$ </li>
 *    <li> therefore, \f$ln(z)\f$ has infinitely many values and is an infinitely many-valued 'function' </li>
 *    <li> that is, \f$ln(z) = ln(r) + i(\theta + 2n \pi)\f$ for all \f$n \in \mathbb{Z}\f$ </li>
 *    <li> note \f$ln(z) = ln(r) + i(\theta + 2n \pi)\f$ implies that \f$\ln(r) = \mbox{Re}(\ln(z))\f$ </li>
 *    <li> therefore since \f$r = \|z\|\f$ we have \f$ln(\|z\|) = ln(r) = \mbox{Re} (ln(z))\f$ </li>
 *  </ul>
 *
 *  @section fmm_2d_section FMM 2D
 *
 *  The Fast Multipole Method is 'fast' since it has a lower floating point operation count in comparison to
 *  the matrix-vector product
 *
 *  \f{eqnarray*}{
 *   \begin{bmatrix} v_1 \\ v_2 \\ v_3 \\ \vdots \\ v_{N_t} \end{bmatrix} =
 *   \begin{bmatrix} \ln(\|y_1-x_1\|) &  \ln(\|y_1-x_1\|) & \ln(\|y_1-x_1\|) & \hdots & \ln(\|y_1-x_{N_s}\|) \\
 *                   \ln(\|y_2-x_1\|) &  \ln(\|y_2-x_1\|) & \ln(\|y_2-x_1\|) & \hdots & \ln(\|y_2-x_{N_s}\|) \\
 *                   \ln(\|y_3-x_1\|) &  \ln(\|y_3-x_1\|) & \ln(\|y_3-x_1\|) & \hdots & \ln(\|y_3-x_{N_s}\|) \\
 *                   \vdots           *  \vdots           & \vdots           &        & \vdots \\
 *                   \ln(\|y_{N_t}-x_1\|) &  \ln(\|y_{N_t}-x_1\|) & \ln(\|y_{N_t}-x_1\|) & \hdots & \ln(\|y_{N_t}-x_{N_s}\|) \\
 *   \end{bmatrix}
 *   \begin{bmatrix} q_1 \\ q_2 \\ q_3 \\ \vdots \\ q_{N_s} \end{bmatrix}
 *  \f}
 *
 *  The key to this 'speed' is approximation.  The Fast Multipole Method uses finite series approximations of operators.
 *  We use Taylor series representations (below) to represent these operations.  Our work gives exact
 *  representations first (with infinite series).  Afterwards, we truncate the infinite series to obatin approximations to the
 *  operators.
 *
 *  The first function that we approximate is the natural logarithm.  We look at two different Taylor series representations
 *  for the function (actually they are the same Taylor series with different centers - one at \f$z\f$ and the other at \f$t\f$).
 *  The reason for using two different series is the radius of convergence.  The radius of convergence for the first series will be
 *  for values \f$y\f$ within a disc of radius \f$R_*\f$.  The second series will be convergent for \f$y\f$-values outside the
 *  disc of radius \f$R_*\f$.  This allows us to take care of all scenarios (except for values on the boundary of the disc,
 *  but we assume that (and somehow will have to make sure that) this does not happen for our problems.
 *
 *  The next series approximations are of translation operations that translate a series center from one location
 *  to another location.
 *  There are three different types of these translations of center that we encounter in the Fast Multipole Method: (i)
 *  near-field to near-field translations (local-to-local translations) called the \f$R|R\f$ translation
 *  (ii) far-field to near-field translations (far-to-near) called the \f$S|R\f$ translation, and (iii) far-field
 *  to far-field translations (far-to-far) called the \f$S|S\f$ translation.
 *
 *  We develop the \f$R|R\f$ translation by writing each power \f$R_m(y-x_*+t) = (y-x_*+t)^m\f$ where \f$m \geq 0\f$
 *  as a Taylor series.  Eventually the derivatives of the function \f$R_m(y-x_*+t)\f$ will be zero and the Taylor series
 *  representation will be finite for each \f$R_m(y-x_*+t)\f$ (a nice feature of this translation).
 *
 *  The \f$S|R\f$ and \f$S|S\f$ translations are developed from the same Taylor series representation but having two different
 *  centers (as mentioned with the \f$R\f$ and \f$S\f$ expansions above).  The Taylor series representations that result
 *  are binomial series and the combination notation for binomial series is used.
 *
 *  We start with the \f$R|R\f$ translation.
 *
 *  @subsection approximating_the_potential_subsection Approximating the Potential Function
 *
 *  From our work above, we see that we can write \f$\phi_{ji} = log( \|y_j - x_i\| )\f$ as
 *  \f$ \phi_{ji} = Re( log( y_j - x_i ) ). \f$
 *  Therefore, we just consider the potential function as
 *  \f[ \phi_{ji} = log( y_j - x_i ). \f]
 *  We approximate the potential function by its Taylor series (see Chapter 15 of \cite kreyszig2010 for a review).
 *  Recall the formula for the Taylor series centered at the complex number \f$z_0\f$.
 *
 *  \f[ f(z) = \sum_{n=0}^{\infty} \frac{f^{(n)}(z_0)}{n!} (z-z_0)^n \f]
 *
 *  Let \f$z = y_j - x_i\f$ and consider the Taylor expansion for \f$f(z) = \log(z) = \log(y_j - x_i)\f$
 *
 *  \f{eqnarray*}{
 *      f(z) & = & \sum_{n=0}^{\infty} \frac{ f^{(n)}(z_0) }{ n! } (z - z_0)^{n} \\
 *           & = & f(z_0) + \frac{f'(z_0)}{1!} (z - z_0)^1 + \frac{1}{2!} f''(z_0)(z-z_0)^2
 *           + \frac{f'''(z_0)}{3!} (z-z_0)^3 + \frac{f^{(4)}(z_0)}{4!} (z-z_0)^4 + ... + \frac{f^{(n)}(z_0)}{n!} (z-z_0)^n + ... \\
 *           & = & ln(z_0) + \frac{1}{1!} \frac{z-z_0}{z_0} + \frac{(-1)^1}{2} \frac{z-z_0}{z_0}
 *                  + \frac{(-1)^2}{3} \frac{(z-z_0)^3}{z_0^3} + \frac{(-1)^3}{4} \frac{(z-z_0)^4}{z_0^4} + ...
 *                  + \frac{(-1)^{n-1}}{n} \frac{(z-z_0)^n}{z_0^n} + ...
 *  \f}
 *
 *  where we have obtained the derivarives from (ln(z))' = 1/z and the power rule applied to complex functions
 *  \f{eqnarray*}{
 *    f(z) &=& ln(z)                                 \quad \quad \implies       f(z_0) = ln(z_0) \\
 *    f'(z) &=& \frac{1}{z} = \frac{(-1)^0 0!}{z^1}    \quad \implies       f'(z_0) = \frac{1}{z_0} \\
 *    f''(z)               &=& \frac{(-1)^1 1!}{z^2} \quad \quad \implies      f''(z_0) = \frac{(-1)^1 1!}{z_0^2} \\
 *    f'''(z)              &=& \frac{(-1)^2 2!}{z^3} \quad \quad \implies      f'''(z_0) = \frac{(-1)^2 2!}{z_0^3} \\
 *    f^{(4)}(z)           &=& \frac{(-1)^3 3!}{z^4} \quad \quad \implies      f^{(4)}(z_0) = \frac{(-1)^3 3!}{z_0^4} \\
 *    &...& \\
 *    f^{(n)}(z) &=& \frac{(-1)^{n-1} (n-1)!}{z^n}   \quad \implies      f^{(n)}(z_0) = \frac{(-1)^{n-1} (n-1)!}{z_0^n}
 *  \f}
 *
 *  From the Ratio Test, we know that the series above \f$\sum_{n=0}^{\infty} a_n\f$ will converge (absolutely)
 *  if
 *  \f[ lim_{n \to \infty} \| \frac{a_{n+1}}{a_{n}} \| = \lim_{n \to \infty} \| \frac{n+1}{n} \frac{z-z_0}{z_0} \| < 1 \implies
 *            \frac{\|z-z_0\|}{\|z_0\|} < 1 \f]
 *
 *  If we take \f$z = y - x_i\f$ and \f$z_0 = y - x_*\f$, then \f$z - z_0 = x_* - x_i\f$ and the Taylor series
 *  \f[ f(z) = ln(z_0) + \sum_{n=1}^{\infty} \frac{ (-1)^{n-1} }{ n } \frac{(z - z_0)^n}{z_0^n} \f]
 *  becomes (where we have changed the letter of the index from \f$n\f$ to \f$m\f$ to correspond with \cite wang2005)
 *  \f{eqnarray*}{
 *    f(y - x_i) &=& ln(y - x_*) + \sum_{n=1}^{\infty} \frac{ (-1)^{n-1} }{ n } \frac{(x_* - x_i)^n}{(y - x_*)^n} \\
 *               &=& ln(y - x_*)
 *                        + \sum_{m=1}^{\infty} \frac{ (-1)^{m-1} * (-1)^{m} }{ m } \frac{(x_i - x_*)^m}{(y - x_*)^m} \\
 *               &=& ln(y - x_*) + \sum_{m=1}^{\infty} \frac{-1}{m} \frac{(x_i-x_*)^m}{(y-x_*)^m}
 *  \f}
 *
 *  The Ratio Test implies that
 *  \f[\frac{ \|z - z_0\| }{ \|z_0\| } = \frac{ \|x_i - x_*\| }{ \|y - x_*\| } < 1. \f]
 *  That is \f$\|x_i - x_*\| < \|y - x_*\|\f$.  Geomtrically, this means the distance between \f$y\f$ and \f$x_*\f$
 *  must be greater than \f$x_i\f$ and \f$x_*\f$, or that \f$x_i\f$ must live inside the disc with center \f$x_*\f$ and radius
 *  \f$\|y - x_*\|\f$.
 *
 *  If instead we let \f$z = y - x_i\f$ and \f$z_0 = x_* - x_i\f$, then \f$z - z_0 = y - x_*\f$ and
 *  \f{eqnarray*}{
 *    f(z) &=& ln(z_0) + \sum_{n=1}^{\infty} \frac{ (-1)^{n-1} }{ n } \frac{(z - z_0)^n}{z_0^n} \\
 *    f(y-x_i) &=& ln(x_*-x_i) + \sum_{n=1}^{\infty} \frac{(-1)^{n-1}}{n} \frac{(y-x_*)^n}{(x_*-x_i)^n} \\
 *             &=& ln(x_*-x_i) + \sum_{m=1}^{\infty} \frac{-1}{m} \frac{(y-x_*)^m}{(x_i-x_*)^m} \\
 *  \f}
 *
 *  The series converges absolutely if \f$\lim_{n \to \infty} \| \frac{a_n+1}{a_n} \| = \| \frac{y-x_*}{x_i-x_*} \| < 1\f$.
 *  That is,
 *  \f[ \|y-x_*\| < \|x_i-x_*\|. \f]
 *  Therefore, for convergence, it is necessary that \f$x_i\f$ live outside the disc with center \f$x_*\f$ and
 *  radius \f$\|y - x_*\|\f$
 *
 *  @subsection field_expansion_subsection S and R Expansions
 *
 *  Far-Field and Near-Field Expansions
 *
 *  @image html fields3.png
 *
 *  @subsubsection far_field_expansion_subsubsection S-Expansion (Far-Field Expansion)
 *
 *  We call the first Taylor series expansion shown above
 *  \f{eqnarray*}{
 *   \ln(y-x_i) &=& \ln(y-x_*) + \sum_{m=1}^{\infty} \frac{-1}{m} \frac{(x_i-x_*)^m}{(y-x_*)^m} \\
 *              &=& \sum_{m=0}^{\infty} b_m(x_i,x_*) S_m(y-x_*)
 *  \f}
 *  the S-expansion or far-field expansion,
 *  where
 *  \f{eqnarray*}{
 *   b_m(x_i,x_*) &=&  \left\{ \begin{matrix}
 *                           1                        & \mbox{if} & m = 0 \\
 *                           \frac{-1}{m} (x_i-x_*)^m & \mbox{if} & m \geq 1
 *                          \end{matrix}
 *                     \right. \\
 *   S_m(y-x_*) &=&    \left\{ \begin{matrix}
 *                               \ln(y-x_*)           & \mbox{if} & m = 0 \\
 *                               \frac{1}{(y-x_*)^m}  & \mbox{if} & m \geq 1
 *                             \end{matrix}
 *                     \right.
 *  \f}
 *
 *  Indeed, from the Ratio Test we have seen that
 *  the Taylor series converges if \f$y\f$ is outside (or far from ) the disc of radius \f$R_* = \|x_i-x_*\|\f$ centered
 *  at \f$x_*\f$.
 *
 *  @subsubsection near_field_expansion_subsubsection R-Expansion (Near-Field Expansion)
 *
 *  We call the second Taylor series expansion above
 *  \f{eqnarray*}{
 *   \ln(y-x_i) &=& \ln(x_*-x_i) + \sum_{m=1}^{\infty} \frac{-1}{m} \frac{(y-x_*)^m}{(x_i-x_*)^m} \\
 *              &=& \sum_{m=0}^{\infty} a_m(x_i,x_*) R_m(y-x_*)
 *  \f}
 *  the R-expansion or near-field expansion (or local expansion),
 *  where
 *  \f{eqnarray*}{
 *   a_m(x_i,x_*) &=&  \left\{ \begin{matrix}
 *                           \ln(x_*-x_i)             & \mbox{if} & m = 0 \\
 *                           \frac{-1}{m (x_i-x_*)^m} & \mbox{if} & m \geq 1
 *                          \end{matrix}
 *                     \right. \\
 *   R_m(y-x_*) &=&    \left\{ \begin{matrix}
 *                               1                    & \mbox{if} & m = 0 \\
 *                               (y-x_*)^m            & \mbox{if} & m \geq 1
 *                             \end{matrix}
 *                     \right.
 *  \f}
 *
 *  And, from the Ratio Test we have seen that
 *  the Taylor series converges if \f$y\f$ lives inside (near or local to) the disc of radius \f$r_* = \|x_i-x_*\|\f$ centered
 *  at \f$x_*\f$.
 *
 *  @subsection translation_subsection Translations
 *
 *  @subsubsection local_to_local_translation_subsubsection (R|R) Translation (local-to-local translation)
 *
 *  Let \f$z = y-x_*\f$ and consider the function \f$f(z) = z^m\f$.  We wish to write the function
 *  \f$f(z+t) = (z+t)^m\f$ as a Taylor series expansion about \f$t = 0\f$.  Let \f$\tilde{t} = z+t\f$
 *  be our complex variable.  Then the function we wish to expand is
 *  \f[ f(\tilde{t}) = (\tilde{t})^m = (z+t)^m = f(z+t) \f]
 *  about \f$\tilde{t} = z\f$. Using the chain rule and the power rule for complex functions on
 *  \f{eqnarray*}{
 *    \frac{d}{d \tilde{t}} f(\tilde{t}) &=& f'(\tilde{t}) = m(\tilde{t})^{m-1} = \frac{m!}{(m-1)!} (\tilde{t})^{m-1}
 *         \implies \quad \frac{df}{d\tilde{t}} \big|_{\tilde{t}=z} =  \frac{m!}{(m-1)!} z^{m-1} \\
 *    \frac{d^2}{dt^2} f(\tilde{t}) &=& f''(z+t) = m(m-1)(\tilde{t})^{m-2} = \frac{m!}{(m-2)!} (\tilde{t})^{m-2}
 *         \quad \implies \quad \frac{d^2f}{d \tilde{t}^2} \big|_{\tilde{t}=z} = \frac{m!}{(m-2)!} z^{m-2} \\
 *    f'''(\tilde{t}) &=& m(m-1)(m-2)(\tilde{t})^{m-3} = \frac{m!}{(m-3)!} (\tilde{t})^{m-3}
 *         \quad \implies \quad f'''(z) = \frac{m!}{(m-3)!} z^{m-3} \\
 *    f^{(4)}(\tilde{t}) &=& m(m-1)(m-2)(m-3)(\tilde{t})^{m-4} = \frac{m!}{(m-4)!} (\tilde{t})^{m-4}
 *         \quad \implies \quad f^{(4)}(z) = \frac{m!}{(m-4)!} z^{m-4} \\
 *    & ... & \\
 *    f^{(m-2)}(\tilde{t}) &=& \frac{m!}{(m-(m-2))!} (\tilde{t})^{m-(m-2)} = \frac{m!}{2!} (\tilde{t})^{2}
 *         \quad \implies \quad f^{(m-2)}(z) = \frac{m!}{2!} z^2 \\
 *    f^{(m-1)}(\tilde{t}) &=& \frac{m!}{(m-(m-1))!} (\tilde{t})^{m-(m-1)} = \frac{m!}{1!} (\tilde{t})^{1}
 *         \quad \implies \quad f^{(m-1)}(z) = \frac{m!}{1!} z^1 = m! z \\
 *    f^{(m)}(\tilde{t})   &=& \frac{m!}{(m-m)!} (\tilde{t})^(m-m) = \frac{m!}{0!} (\tilde{t})^0 = m!
 *         \quad \implies \quad f^{(m)}(z) = m! \\
 *    f^{(m+1)}(\tilde{t}) &=& \frac{d}{dt} m! = 0
 *         \quad \implies \quad f^{(m+1)}(z) = 0
 *  \f}
 *
 *  Therefore, the expansion is
 *  \f[ (y-(x_*-t))^m = (y-x_*+t)^m = (z+t)^m = (\tilde{t})^m = f(\tilde{t}) = \sum_{n=0}^{\infty} \frac{f^{(n)}(z)}{n!} (\tilde{t} - z)^n
 *              = \sum_{n=0}^{m} \frac{m!}{(m-n)!n!} z^{m-n} t^n
 *              = \sum_{n=0}^{m} \begin{pmatrix} m \\ n \end{pmatrix} z^{m-n} t^n
 *  \f]
 *
 *  Letting the expansion being centered at \f$t\f$ and noting that \f$\tilde{t}-t = z\f$, then \f$(y-(x_*-t))^m\f$ can be
 *  written as
 *  \f[ (y-(x_*-t))^m = \sum_{n=0}^{m} t^{m-n} z^n = \sum_{n=0}^{m} (y-x_*)^n t^{m-n} \f]
 *
 *  Note \f$R_m(y-x_*+t) = (y-x_*+t)^m\f$.  Therefore,
 *  \f{eqnarray*}{
 *    R_0(y-x_*+t) & = & \sum_{n=0}^{0} \begin{pmatrix} 0 \\ n \end{pmatrix} (y - x_*)^n t^{m-n} = 1 = R_0(y-x_*) \\
 *    R_1(y-x_*+t) & = & \sum_{n=0}^{1} \begin{pmatrix} 1 \\ n \end{pmatrix} (y-x_*)^n t^{m-n}
 *                   =   \begin{pmatrix} 1 \\ 0 \end{pmatrix} (y-x_*)^0 t^1
 *                       + \begin{pmatrix} 1 \\ 1 \end{pmatrix} (y-x_*)^1 t^0 \\
 *                 & = & t + (y-x_*) \\
 *                 & = & tR_0(y-x_*) + R_1(y-x_*) \\
 *    R_2(y-x_*+t) & = & \sum_{n=0}^{2} \begin{pmatrix} 2 \\ n \end{pmatrix} (y - x_*)^n t^{m-n} \\
 *                 & = & \begin{pmatrix} 2 \\ 0 \end{pmatrix} (y-x_*)^0 t^2
 *                       + \begin{pmatrix} 2 \\ 1 \end{pmatrix} (y-x_*)^1 t^1
 *                       + \begin{pmatrix} 2 \\ 2 \end{pmatrix} (y-x_*)^2 t^0 \\
 *                 & = & t^2 + 2(y-x_*)t + (y-x_*)^2 \\
 *                 & = & t^2R_0(y-x_*) + 2tR_1(y-x_*) + R_2(y-x_*) \\
 *    R_3(y-x_*+t) & = & \sum_{n=0}^{3} \begin{pmatrix} 3 \\ n \end{pmatrix} (y - x_*)^n t^{m-n} \\
 *                 & = & \begin{pmatrix} 3 \\ 0 \end{pmatrix} (y-x_*)^0 t^3
 *                       + \begin{pmatrix} 3 \\ 1 \end{pmatrix} (y-x_*)^1 t^2
 *                       + \begin{pmatrix} 3 \\ 2 \end{pmatrix} (y-x_*)^2 t^1
 *                       + \begin{pmatrix} 3 \\ 3 \end{pmatrix} (y-x_*)^3 t^0 \\
 *                 & = & t^3 + 3(y-x_*)t^2 + 3(y-x_*)^2t + (y-x_*)^3 \\
 *                 & = & t^3R_0(y-x_*) + 3t^2R_1(y-x_*) + 3tR_2(y-x_*) + R_3(y-x_*) \\
 *    R_4(y-x_*+t) & = & \sum_{n=0}^{4} \begin{pmatrix} 4 \\ n \end{pmatrix} (y - x_*)^n t^{m-n} \\
 *                 & = & \begin{pmatrix} 4 \\ 0 \end{pmatrix} (y-x_*)^0 t^4
 *                       + \begin{pmatrix} 4 \\ 1 \end{pmatrix} (y-x_*)^1 t^3
 *                       + \begin{pmatrix} 4 \\ 2 \end{pmatrix} (y-x_*)^2 t^2
 *                       + \begin{pmatrix} 4 \\ 3 \end{pmatrix} (y-x_*)^3 t^1
 *                       + \begin{pmatrix} 4 \\ 4 \end{pmatrix} (y-x_*)^4 t^0 \\
 *                 & = & t^4 + 4(y-x_*)t^3 + 6(y-x_*)^2t^2 + 4(y-x_*)^3t + (y-x_*)^4 \\
 *                 & = & t^4R_0(y-x_*) + 4t^3R_1(y-x_*) + 6t^2R_2(y-x_*) + 4tR_3(y-x_*) + R_4(y-x_*) \\
 *    R_5(y-x_*+t) & = & \sum_{n=0}^{5} \begin{pmatrix} 5 \\ n \end{pmatrix} (y - x_*)^n t^{m-n} \\
 *                 & = & \begin{pmatrix} 5 \\ 0 \end{pmatrix} (y-x_*)^0 t^5
 *                       + \begin{pmatrix} 5 \\ 1 \end{pmatrix} (y-x_*)^1 t^4
 *                       + \begin{pmatrix} 5 \\ 2 \end{pmatrix} (y-x_*)^2 t^3
 *                       + \begin{pmatrix} 5 \\ 3 \end{pmatrix} (y-x_*)^3 t^2
 *                       + \begin{pmatrix} 5 \\ 4 \end{pmatrix} (y-x_*)^4 t^1
 *                       + \begin{pmatrix} 5 \\ 5 \end{pmatrix} (y-x_*)^5 t^0 \\
 *                 & = & t^5 + 5(y-x_*)t^4 + 10(y-x_*)^2t^3 + 10(y-x_*)^3t^2 + 5(y-x_*)^4t + (y-x_*)^5 \\
 *                 & = & t^5R_0(y-x_*) + 5t^4R_1(y-x_*) + 10t^3R_2(y-x_*) + 10t^2R_3(y-x_*) + 5tR_4(y-x_*) + R_5(y-x_*) \\
 *                 & \vdots &
 *  \f}
 *
 *  As a transformation, this results in an infinite matrix
 *
 *  \f[  \begin{bmatrix} R_0(y-x_*+t) \\ R_1(y-x_*+t) \\ R_2(y-x_*+t) \\ R_3(y-x_*+t) \\ R_4(y-x_*+t) \\ R_5(y-x_*+t) \\
 *                       R_6(y-x_*+t) \\ \vdots
 *       \end{bmatrix}
 *       =
 *      \begin{bmatrix} 1      & 0      & 0      & 0      & 0      & 0      & 0      & \hdots & 0      & \hdots \\
 *                      t      & 1      & 0      & 0      & 0      & 0      & 0      & \hdots & 0      & \hdots  \\
 *                      t^2    & 2t     & 1      & 0      & 0      & 0      & 0      & \hdots & 0      & \hdots  \\
 *                      t^3    & 3t^2   & 3t     & 1      & 0      & 0      & 0      & \hdots & 0      & \hdots  \\
 *                      t^4    & 4t^3   & 6t^2   & 4t     & 1      & 0      & 0      & \hdots & 0      & \hdots  \\
 *                      t^5    & 5t^4   & 10t^3  & 10t^2  & 5t     & 1      & 0      & \hdots & 0      & \hdots  \\
 *                      \vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \ddots & \ddots & \vdots & \vdots  \\
 *      \end{bmatrix}
 *      \begin{bmatrix} R_0(y-x_*) \\ R_1(y-x_*) \\ R_2(y-x_*) \\ R_3(y-x_*) \\ R_4(y-x_*) \\ R_5(y-x_*) \\
 *                      R_6(y-x_*) \\ \vdots  \end{bmatrix}
 *  \f]
 *
 *  To obtain the formula for the translation from the center \f$x_*-t\f$ to \f$x_*\f$, we multiply each equation by its corresponding
 *  coefficient \f$a_m(x_i,x_*-t)\f$ and sum the equations.  The left-hand side is the Taylor series expansion at \f$(x_*-t)\f$
 *  \f[ \sum_{m=0}^{\infty} a_m(x_i,x_*-t) R_m(y-x_*+t) \f]
 *  and the right-hand side is the formula for the translation
 *
 *  \f{eqnarray*}{
 *   \sum_{m=0}^{\infty} a_m(x_i,x_*-t) R_m(y-x_*+t)
 *     & = & \sum_{k=0}^{\infty} a_k(x_i,x_*-t) \begin{pmatrix} k \\ 0 \end{pmatrix} t^k R_0(y-x_*)
 *           + \sum_{k=1}^{\infty} a_{k}(x_i,x_*-t) \begin{pmatrix} k \\ 1 \end{pmatrix}  t^{k-1} R_1(y-x_*) \\
 *     &&    + \sum_{k=2}^{\infty} a_{k}(x_i,x_*-t) \begin{pmatrix} k \\ 2 \end{pmatrix}  t^{k-2} R_2(y-x_*)
 *           + \sum_{k=3}^{\infty} a_{k}(x_i,x_*-t) \begin{pmatrix} k \\ 3 \end{pmatrix}  t^{k-3} R_3(y-x_*) \\
 *     &&    + \hdots
 *           + \sum_{k=l-1}^{\infty} a_{k}(x_i,x_*-t) \begin{pmatrix} k \\ l-1 \end{pmatrix}  t^{k-(l-1)} R_{l-1}(y-x_*)
 *           + \sum_{k=l}^{\infty} a_{k}(x_i,x_*-t) \begin{pmatrix} k \\ l \end{pmatrix}  t^{k-l} R_l(y-x_*) + \hdots \\
 *     & = & \sum_{m=0}^{\infty} \left( \sum_{k=m}^{\infty} a_k(x_i,x_*-t) \begin{pmatrix} k \\ m \end{pmatrix} t^{k-m} \right) R_m(y-x_*)
 *  \f}
 *
 *  Therefore, we have the formula for the translation
 *  \f[ \sum_{m=0}^{\infty} a_m(x_i,x_*-t) R_m(y-x_*+t) = \sum_{m=0}^{\infty} a_m(x_i,x_*) R_m(y-x_*) \f]
 *  where
 *  \f[ a_m(x_i,x_*) = \sum_{k=m}^{\infty} a_k(x_i,x_*-t) \begin{pmatrix} k \\ m \end{pmatrix} t^{k-m} \f]
 *
 *
 *  <b> Approximation </b>
 *
 *  We will appproximate the RR-translation by truncating the series to \f$p\f$ terms.
 *  This results in a \f$p\f$-by-\f$p\f$ translation matrix.  Note that the indexing
 *  for \f$p\f$ starts on one.
 *
 *
 *  \f{eqnarray*}{
 *      \begin{bmatrix} 1       & 0      & 0      & 0      & 0      & 0      & \hdots & 0  \\
 *                      t       & 1      & 0      & 0      & 0      & 0      & \hdots & 0  \\
 *                      t^2     & 2t     & 1      & 0      & 0      & 0      & \hdots & 0  \\
 *                      t^3     & 3t^2   & 3t     & 1      & 0      & 0      & \hdots & 0  \\
 *                      t^4     & 4t^3   & 6t^2   & 4t     & 1      & 0      & \hdots & 0  \\
 *                      t^5     & 5t^4   & 10t^3  & 10t^2  & 5t     & 1      & \ddots & \vdots  \\
 *                      \vdots  & \vdots & \vdots & \vdots & \vdots & \ddots & \ddots & 0  \\
 *                      t^{p-1} & \begin{pmatrix} p-1 \\ 1 \end{pmatrix} t^{p-2}
 *                                       & \begin{pmatrix} p-1 \\ 2 \end{pmatrix} t^{p-3}
 *                                                & \begin{pmatrix} p-1 \\ 3 \end{pmatrix} t^{p-4}
 *                                                         & \begin{pmatrix} p-1 \\ 4 \end{pmatrix} t^{p-5}
 *                                                                  & \hdots
 *                                                                           & \begin{pmatrix} p-1 \\ p-2 \end{pmatrix} t^{p-(p-1)}
 *                                                                                    & 1
 *      \end{bmatrix}
 *  \f}
 *
 *  We work with the transpose of this approximation in the code (getRR member function in Potential.cc).
 *
 *
 *
 *
 *  @subsubsection far_to_near_translation_subsubsection (S|R) Translation (far-to-near translation)
 *
 *  Recall the Far-Field Expansion
 *  \f{eqnarray*}{
 *   \ln(y-x_i) &=& ln(y-x_*) + \sum_{m=1}^{\infty} \frac{-1}{m} \frac{(x_i-x_*)^m}{(y-x_*)^m} \\
 *              & = & \sum_{m=0}^{\infty} b_m(x_i,x_*) S_m(y-x_*)
 *  \f}
 *
 *  Let \f$z = y-x_*\f$ and consider the function \f$f(z) = \frac{1}{z^m} \f$.
 *  Let \f$\tilde{t} = z+t and \f$k = -m\f$, then
 *  \f[ f(\tilde{t}) = (\tilde{t})^k = (z+t)^{-m} = f(z+t). \f]
 *  We wish to determine the Taylor series centered at \f$\tilde{t} = z\f$.
 *  Using the chain rule and the power rule for complex functions on
 *
 *  \f{eqnarray*}{
 *    f (\tilde{t}) &=& \tilde{t}^{-m}
 *          \quad \implies \quad
 *          f(z) = z^{-m} \\
 *    f'  (\tilde{t}) &=& -m \tilde{t}^{-m-1} = -m \tilde{t}^{-m-1}
 *          \quad \implies \quad  f'(z) = -m z^{-m-1} \\
 *    f'' (\tilde{t}) &=& -m(-m-1) \tilde{t}^{-m-2} = -m(-m-1) \tilde{t}^{-m-2} = \frac{(-m)!}{(-m-2)!} \tilde{t}^{-m-2}
 *          \quad \implies \quad  f''(z) = \frac{(-m)!}{(-m-2)!} z^{-m-2} \\
 *    f'''(\tilde{t}) &=& -m(-m-1)(-m-2) \tilde{t}^{-m-3} = \frac{(-m)!}{(-m-3)!} \tilde{t}^{-m-3}
 *          \quad \implies \quad  f'''(z) = \frac{(-m)!}{(-m-3)!} z^{-m-3} = \\
 *    f^{(4)} (\tilde{t}) & = & -m(-m-1)(-m-2)(-m-3) \tilde{t}^{-m-4} = \frac{(-m)!}{(-m-4)!} z^{-m-4}
 *           \quad \implies \quad  f^{(4)}(z) = (-1)^4 \frac{(-m)!}{(-m-4)!} z^{-m-4} \\
 *    f^{(5)} (\tilde{t}) & = & -m(-m-1)(-m-2)(-m-3)(-m-4) \tilde{t}^{-(m+5)} = (-1)^5 \frac{(-m)!}{(-m-5)!} \tilde{t}^{-m-5}
 *          \quad \implies \quad f^{(5)}(z) = \frac{(-m)!}{(-m-5)!} z^{-m-5} \\
 *          & \vdots &
 *  \f}
 *
 *  Realizing that we are dealing with a binomial series (a special Taylor series), let \f$k = -m\f$ and write the expansion as
 *
 *  \f{eqnarray*}{
 *    \frac{1}{\tilde{t}^m} = \tilde{t}^k = f(\tilde{t}) & = & \sum_{n=0}^{\infty} \frac{f^{(n)}(z)}{n!} (\tilde{t} - z)^n
 *                   = \sum_{n=0}^{\infty}  \frac{k!}{(k-n)!n!} z^{k-n}t^n
 *                   = \sum_{n=0}^{\infty} \begin{pmatrix} k \\ n \end{pmatrix} z^{k-n} t^n
 *                   = \sum_{n=0}^{\infty} \begin{pmatrix} -m \\ n \end{pmatrix} z^{-m-n} t^n \\
 *                 & = & \sum_{n=0}^{\infty} \begin{pmatrix} -m \\ n \end{pmatrix} z^n t^{-m-n}
 *  \f}
 *
 *  where the last equality results from the center of the expansion around \f$t\f$ where \f$\tilde{t}-t = z\f$.
 *  We use the combination notation loosely in the sense that factorials of negative integers are not
 *  formally defined here.  However, the notation is the same as that used for the binomial series
 *  \f[ \begin{pmatrix} -m \\ n \end{pmatrix} \colon =
 *               \frac{(-m)(-m-1)(-m-2)\hdots(-m-n)}{n!} \quad \left( =: \frac{(-m)!}{(-m-n)!n!} \right) \f]
 *  Recall
 *  \f[ S_m(y-x_*) = \left\{ \begin{matrix} \ln(y-x_*) & \mbox{if} & m = 0 \\
 *                                             \frac{1}{(y-x_*)^m} & \mbox{if} & m \geq 1
 *                              \end{matrix}
 *                      \right.
 *  \f]
 *
 *  With \f$S_m(y-x_*+t)\f$ defined in the same manner (recalling \f$z = y-x_*\f$ and doing the Taylor series for the
 *  case \f$S_0(y-x_*+t)\f$ on the fly)
 *  \f{eqnarray*}{
 *    S_0(y-x_*+t) &=& \ln(y-x_*+t)
 *                  =  \ln(t) + \frac{1}{t} z - \frac{1}{2t^2} z^2 + \frac{1}{3t^3} z^3 - \frac{1}{4t^4} z^4 +
 *                     + \hdots + \frac{(-1)^n}{nt^n}z^n + \hdots \\
 *                 &=& \ln(t) R_0(y-x_*) + \frac{1}{t} R_1(y-x_*) - \frac{1}{2t^2} R_2(y-x_*) +
 *                      + \frac{1}{3t^3} R_3(y-x_*) - \frac{1}{4t^4} R_4(y-x_*) + \hdots + \frac{(-1)^n}{4t^n} R_n(y-x_*) + \hdots \\
 *    S_1(y-x_*+t) &=& \frac{1}{y-x_*+t} = \sum_{n=0}^{\infty} \begin{pmatrix} -1 \\ n \end{pmatrix} (y-x_*)^n t^{-1-n}
 *                  =  \sum_{n=0}^{\infty} \begin{pmatrix} -1 \\ n \end{pmatrix} R_n(y-x_*) t^{-1-n} \\
 *                 & = &
 *                   \begin{pmatrix} -1 \\ 0 \end{pmatrix} R_0(y-x_*) t^{-1-0}
 *                   + \begin{pmatrix} -1 \\ 1 \end{pmatrix} R_1(y-x_*) t^{-1-1}
 *                   + \begin{pmatrix} -1 \\ 2 \end{pmatrix} R_2(y-x_*) t^{-1-2}
 *                   + \hdots
 *                   + \begin{pmatrix} -1 \\ n \end{pmatrix} R_n(y-x_*) t^{-1-n}
 *                   + \hdots \\
 *                 & = &
 *                   \frac{(-1)!}{(-1-0)!0!} R_0(y-x^*) t^{-1-0} + \frac{(-1)!}{(-1-1)!1!} R_1(y-x_*)t^{-1-1}
 *                    + \frac{(-1)!}{(-1-2)!2!} R_2(y-x_*)t^{-1-2} + \hdots
 *                    + \frac{(-1)!}{(-1-n)!n!} R_n(y-x_*)t^{-1-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-1} + \frac{(-1)(-2)!}{(-2)!1!} R_1(y-x_*)t^{-2} +
 *                    + \frac{(-1)(-2)(-3)!}{(-3)!2!} R_2(y-x_*)t^{-3} +
 *                    + \frac{(-1)(-2)(-3)(-4)!}{(-4)!3!} R_3(y-x_*)t^{-4} \\
 *                 && + \frac{(-1)(-2)(-3)(-4)(-5)!}{(-5)!4!} R_4(y-x_*)t^{-5} + \hdots
 *                    + \frac{(-1)(-2)\hdots(-n)(-1-n)!}{(-1-n)!n!} R_n(y-x_*)t^{-1-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-1} - R_1(y-x_*)t^{-2} + R_2(y-x_*)t^{-3}
 *                   - R_3(y-x_*)t^{-4} + R_4(y-x_*)t^{-5} + \hdots
 *                   + (-1)^{n}R_n(y-x_*)t^{-1-n} + \hdots \\
 *    S_2(y-x_*+t) &=& \frac{1}{(y-x_*+t)^2} = \sum_{n=0}^{\infty} \begin{pmatrix} -2 \\ n \end{pmatrix} (y-x_*)^n t^{-2-n}
 *                  =  \sum_{n=0}^{\infty} \begin{pmatrix} -2 \\ n \end{pmatrix} R_n(y-x_*) t^{-2-n} \\
 *                 & = &
 *                   \begin{pmatrix} -2 \\ 0 \end{pmatrix} R_0(y-x_*) t^{-2-0}
 *                   + \begin{pmatrix} -2 \\ 1 \end{pmatrix} R_1(y-x_*) t^{-2-1}
 *                   + \begin{pmatrix} -2 \\ 2 \end{pmatrix} R_2(y-x_*) t^{-2-2}
 *                   + \hdots
 *                   + \begin{pmatrix} -2 \\ n \end{pmatrix} R_n(y-x_*) t^{-2-n}
 *                   + \hdots \\
 *                 & = &
 *                   \frac{(-2)!}{(-2-0)!0!} R_0(y-x^*) t^{-2-0} + \frac{(-2)!}{(-2-1)!1!} R_2(y-x_*)t^{-2-1}
 *                    + \frac{(-2)!}{(-2-2)!2!} R_2(y-x_*)t^{-2-2} + \hdots
 *                    + \frac{(-1)!}{(-1-n)!n!} R_n(y-x_*)t^{-1-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-2} + \frac{(-2)(-3)!}{(-3)!1!} R_1(y-x_*)t^{-3} +
 *                    + \frac{(-2)(-3)(-4)!}{(-4)!2!} R_2(y-x_*)t^{-4} +
 *                    + \frac{(-2)(-3)(-4)(-5)!}{(-5)!3!} R_3(y-x_*)t^{-5} \\
 *                 && + \frac{(-2)(-3)(-4)(-5)(-6)!}{(-6)!4!} R_4(y-x_*)t^{-6} + \hdots
 *                    + \frac{(-2)(-3)\hdots(-1-n)(-2-n)!}{(-2-n)!n!} R_n(y-x_*)t^{-1-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-2} - 2R_1(y-x_*)t^{-3} + 3R_2(y-x_*)t^{-4}
 *                   - 4R_3(y-x_*)t^{-5} + 5R_4(y-x_*)t^{-6} + \hdots
 *                   + (-1)^{n}(n+1)R_n(y-x_*)t^{-2-n} + \hdots
 *  \f}
 *  \f{eqnarray*}{
 *    S_3(y-x_*+t) &=& \frac{1}{(y-x_*+t)^3} = \sum_{n=0}^{\infty} \begin{pmatrix} -3 \\ n \end{pmatrix} (y-x_*)^n t^{-3-n}
 *                  =  \sum_{n=0}^{\infty} \begin{pmatrix} -3 \\ n \end{pmatrix} R_n(y-x_*) t^{-3-n} \\
 *                 & = &
 *                   \begin{pmatrix} -3 \\ 0 \end{pmatrix} R_0(y-x_*) t^{-2-0}
 *                   + \begin{pmatrix} -3 \\ 1 \end{pmatrix} R_1(y-x_*) t^{-2-1}
 *                   + \begin{pmatrix} -3 \\ 2 \end{pmatrix} R_2(y-x_*) t^{-2-2}
 *                   + \hdots
 *                   + \begin{pmatrix} -3 \\ n \end{pmatrix} R_n(y-x_*) t^{-2-n}
 *                   + \hdots \\
 *                 & = &
 *                   \frac{(-3)!}{(-3-0)!0!} R_0(y-x^*) t^{-3-0} + \frac{(-3)!}{(-3-1)!1!} R_2(y-x_*)t^{-3-1}
 *                    + \frac{(-3)!}{(-3-2)!2!} R_2(y-x_*)t^{-3-2} + \hdots
 *                    + \frac{(-3)!}{(-3-n)!n!} R_n(y-x_*)t^{-3-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-3} + \frac{(-3)(-4)!}{(-4)!1!} R_1(y-x_*)t^{-4} +
 *                    + \frac{(-3)(-4)(-5)!}{(-5)!2!} R_2(y-x_*)t^{-5} +
 *                    + \frac{(-3)(-4)(-5)(-6)!}{(-6)!3!} R_3(y-x_*)t^{-6} \\
 *                 && + \frac{(-3)(-4)(-5)(-6)(-7)!}{(-7)!4!} R_4(y-x_*)t^{-7} + \hdots
 *                    + \frac{(-3)(-4)\hdots(-2-n)(-3-n)!}{(-3-n)!n!} R_n(y-x_*)t^{-3-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-3} - 3R_1(y-x_*)t^{-4} + 6R_2(y-x_*)t^{-5}
 *                   - 10R_3(y-x_*)t^{-6} + 15R_4(y-x_*)t^{-7} + \hdots
 *                   + (-1)^{n}\frac{(n+2)(n+1)}{2} R_n(y-x_*)t^{-3-n} + \hdots \\
 *                 & = &
 *                   (-1)^0 \begin{pmatrix} 2 \\ 2 \end{pmatrix} R_0(y-x_*) t^{-3}
 *                   +(-1)^1 \begin{pmatrix} 3 \\ 2 \end{pmatrix} R_1(y-x_*) t^{-4}
 *                   +(-1)^2 \begin{pmatrix} 4 \\ 2 \end{pmatrix} R_2(y-x_*) t^{-5}
 *                   +(-1)^3 \begin{pmatrix} 5 \\ 2 \end{pmatrix} R_3(y-x_*) t^{-6} \\
 *                 &&  +(-1)^4 \begin{pmatrix} 6 \\ 2 \end{pmatrix} R_4(y-x_*) t^{-7}
 *                     + \hdots
 *                     + (-1)^n \begin{pmatrix} n+2 \\ 2 \end{pmatrix} R_n(y-x_*) t^{-3-n} + \hdots \\
 *    S_4(y-x_*+t) &=& \frac{1}{(y-x_*+t)^4} = \sum_{n=0}^{\infty} \begin{pmatrix} -4 \\ n \end{pmatrix} (y-x_*)^n t^{-4-n} \\
 *                 & = &
 *                   \begin{pmatrix} -4 \\ 0 \end{pmatrix} R_0(y-x_*) t^{-4-0}
 *                   + \begin{pmatrix} -4 \\ 1 \end{pmatrix} R_1(y-x_*) t^{-4-1}
 *                   + \begin{pmatrix} -4 \\ 2 \end{pmatrix} R_2(y-x_*) t^{-4-2}
 *                   + \hdots
 *                   + \begin{pmatrix} -4 \\ n \end{pmatrix} R_n(y-x_*) t^{-4-n}
 *                   + \hdots \\
 *                 & = &
 *                   \frac{(-4)!}{(-4-0)!0!} R_0(y-x^*) t^{-4-0} + \frac{(-4)!}{(-4-1)!1!} R_2(y-x_*)t^{-4-1}
 *                    + \frac{(-4)!}{(-4-2)!2!} R_2(y-x_*)t^{-4-2} + \hdots
 *                    + \frac{(-4)!}{(-4-n)!n!} R_n(y-x_*)t^{-4-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-4} + \frac{(-4)(-5)!}{(-5)!1!} R_1(y-x_*)t^{-5} +
 *                    + \frac{(-4)(-5)(-6)!}{(-6)!2!} R_2(y-x_*)t^{-6} +
 *                    + \frac{(-4)(-5)(-6)(-7)!}{(-7)!3!} R_3(y-x_*)t^{-7} \\
 *                 && + \frac{(-4)(-5)(-6)(-7)(-8)!}{(-8)!4!} R_4(y-x_*)t^{-8} + \hdots
 *                    + \frac{(-4)(-5)\hdots(-3-n)(-4-n)!}{(-4-n)!n!} R_n(y-x_*)t^{-4-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-4} - 4R_1(y-x_*)t^{-5} + 10R_2(y-x_*)t^{-6}
 *                   - 20R_3(y-x_*)t^{-7} + 35R_4(y-x_*)t^{-8} + \hdots
 *                   + (-1)^{n}\frac{(n+3)(n+2)(n+1)}{3!} R_n(y-x_*)t^{-4-n} + \hdots \\
 *                 & = &
 *                   (-1)^0 \begin{pmatrix} 3 \\ 3 \end{pmatrix} R_0(y-x_*)t^{-4}
 *                    + (-1)^1 \begin{pmatrix} 4 \\ 3 \end{pmatrix} R_1(y-x_*)t^{-5}
 *                    + (-1)^2 \begin{pmatrix} 5 \\ 3 \end{pmatrix} R_2(y-x_*)t^{-6}
 *                    + (-1)^3 \begin{pmatrix} 6 \\ 3 \end{pmatrix} R_3(y-x_*)t^{-7} \\
 *                 && + (-1)^4 \begin{pmatrix} 7 \\ 3 \end{pmatrix} R_4(y-x_*)t^{-8}
 *                   + \hdots
 *                   + (-1)^n \begin{pmatrix} n+3 \\ 3 \end{pmatrix} R_n(y-x_*)t^{-4-n} + \hdots
 *    \f}
 *    \f{eqnarray*}{
 *    S_5(y-x_*+t) &=& \frac{1}{(y-x_*+t)^5} = \sum_{n=0}^{\infty} \begin{pmatrix} -5 \\ n \end{pmatrix} (y-x_*)^n t^{-5-n} \\
 *                 & = &
 *                   \begin{pmatrix} -5 \\ 0 \end{pmatrix} R_0(y-x_*) t^{-5-0}
 *                   + \begin{pmatrix} -5 \\ 1 \end{pmatrix} R_1(y-x_*) t^{-5-1}
 *                   + \begin{pmatrix} -5 \\ 2 \end{pmatrix} R_2(y-x_*) t^{-5-2}
 *                   + \hdots
 *                   + \begin{pmatrix} -5 \\ n \end{pmatrix} R_n(y-x_*) t^{-5-n}
 *                   + \hdots \\
 *                 & = &
 *                   \frac{(-5)!}{(-5-0)!0!} R_0(y-x^*) t^{-5-0} + \frac{(-5)!}{(-5-1)!1!} R_2(y-x_*)t^{-5-1}
 *                    + \frac{(-5)!}{(-5-2)!2!} R_2(y-x_*)t^{-5-2} + \hdots
 *                    + \frac{(-5)!}{(-5-n)!n!} R_n(y-x_*)t^{-5-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-5} + \frac{(-5)(-6)!}{(-6)!1!} R_1(y-x_*)t^{-6} +
 *                    + \frac{(-5)(-6)(-7)!}{(-7)!2!} R_2(y-x_*)t^{-7} +
 *                    + \frac{(-5)(-6)(-7)(-8)!}{(-8)!3!} R_3(y-x_*)t^{-8} \\
 *                 && + \frac{(-5)(-6)(-7)(-8)(-9)!}{(-9)!4!} R_4(y-x_*)t^{-9} + \hdots
 *                    + \frac{(-5)(-6)\hdots(-4-n)(-5-n)!}{(-5-n)!n!} R_n(y-x_*)t^{-5-n} + \hdots \\
 *                 & = &
 *                   R_0(y-x_*) t^{-5} - 5R_1(y-x_*)t^{-6} + 15R_2(y-x_*)t^{-7}
 *                   - 35R_3(y-x_*)t^{-8} + 70R_4(y-x_*)t^{-9} + \hdots
 *                   + (-1)^{n}\frac{(n+4)(n+3)(n+2)(n+1)}{4!} R_n(y-x_*)t^{-5-n} + \hdots \\
 *                 & = &
 *                   (-1)^0 \begin{pmatrix} 4 \\ 4 \end{pmatrix} R_0(y-x_*)t^{-5}
 *                    + (-1)^1 \begin{pmatrix} 5 \\ 4 \end{pmatrix} R_1(y-x_*)t^{-6}
 *                    + (-1)^2 \begin{pmatrix} 6 \\ 4 \end{pmatrix} R_2(y-x_*)t^{-7}
 *                    + (-1)^3 \begin{pmatrix} 7 \\ 4 \end{pmatrix} R_3(y-x_*)t^{-8} \\
 *                 && + (-1)^4 \begin{pmatrix} 8 \\ 4 \end{pmatrix} R_4(y-x_*)t^{-9}
 *                   + \hdots
 *                   + (-1)^n \begin{pmatrix} n+4 \\ 4 \end{pmatrix} R_n(y-x_*)t^{-5-n} + \hdots \\
 *                 & \vdots &
 *  \f}
 *
 *  The equations define a transformation matrix that is infinite
 *  \f[ \begin{bmatrix} S_0(y-x_*+t) \\ S_1(y-x_*+t) \\ S_2(y-x_*+t) \\ S_3(y-x_*+t) \\ S_4(y-x_*+t) \\ S_5(y-x_*+t) \\
 *                      \vdots \\ S_n(y-x_*+t) \\ \vdots \\ \end{bmatrix}
 *      =
 *      \begin{bmatrix} \ln(t) & \frac{1}{t} & -\frac{1}{2t^2} & \frac{1}{3t^3} & -\frac{1}{4t^4} & \hdots
 *                      & \frac{(-1)^n}{nt^n}                                              & \hdots \\
 *                      t^{-1} & -t^{-2}  & t^{-3}  & -t^{-4}   & t^{-5}   & \hdots
 *                      & (-1)^n t^{-1-n}                                                  & \hdots  \\
 *                      t^{-2} & -2t^{-3} & 3t^{-4} & -4t^{-5}  & 5t^{-6}  & \hdots
 *                      & (-1)^n \begin{pmatrix} n+1 \\ 1 \end{pmatrix} t^{-2-n}        & \hdots  \\
 *                      t^{-3} & -3t^{-4} & 6t^{-5} & -10t^{-6} & 15t^{-7} & \hdots
 *                      & (-1)^n \begin{pmatrix} n+2 \\ 2 \end{pmatrix} t^{-3-n} & \hdots \\
 *                      t^{-4} & -4t^{-5} & 10t^{-6} & -20t^{-7} & 35t^{-8} & \hdots
 *                      & (-1)^n \begin{pmatrix} n+3 \\ 3 \end{pmatrix} t^{-4-n} & \hdots \\
 *                      t^{-5} & -5t^{-6} & 15t^{-7} & -35t^{-8} & 70t^{-9} & \hdots
 *                      & (-1)^n \begin{pmatrix} n+4 \\ 4 \end{pmatrix} t^{-5-n} & \hdots \\
 *                      \vdots & \vdots   & \vdots   & \vdots    & \vdots   & \ddots
 *                      & \vdots                                                           & \hdots \\
 *      \end{bmatrix}
 *      \begin{bmatrix} R_0(y-x_*) \\ R_1(y-x_*) \\ R_2(y-x_*) \\ R_3(y-x_*) \\ R_4(y-x_*) \\ R_5(y-x_*) \\
 *                      \vdots \\ R_n(y-x_*) \\ \vdots \\ \end{bmatrix}
 *  \f]
 *
 *  Note that each expansion of \f$S_m(y-x_*+t)\f$ corresponds to a row of the infinite matrix.
 *  To form
 *  \f[ \ln(y-x_i) = \sum_{m=0}^{\infty} b_m(x_i,x_*-t) S_m(y-x_*+t) \f]
 *  in terms of \f$R_m(y-x_*)\f$ (the \f$S|R\f$ translation) we multiply each equation above
 *  by the corresponding \f$b_m(x_i,x_*-t)\f$ and add the equations.
 *
 *  We write the summation with respect to the \f$R_m(y-x_*)\f$s.
 *  Since each column of the infinite matrix corresponds to a particular \f$R_m(y-x_*)\f$, we seek a formula for the
 *  summation of each column.
 *  Formulas for each column can be found.
 *  However, for this translation the first row is slightly different than its neighbors in the corresponding columns.
 *  In particular, a formula for the first column (column 0) is difficult to determine with the first element being
 *  \f$\ln(t)\f$ (located at (1,1) or (0,0) position).
 *  Due to this problem, the terms involving the first row have been split off from the rest of the summation (as shown below).
 *
 *  \f{eqnarray*}{
 *    \sum_{m=0}^{\infty} b_m(x_i,x_*-t) S_m(y-x_*+t)
 *               & = & b_0(x_i,x_*-t) \ln(t) R_0(y-x_*)
 *                     + \sum_{n=1}^{\infty} b_0(x_i,x_*-t) \frac{(-1)^{n+1}}{n} \frac{1}{t^n} R_m(y-x_*)
 *                   + \sum_{m=0}^{\infty}
 *                      \left( \sum_{n=1}^{\infty} b_n(x_i,x_*-t) (-1)^m \begin{pmatrix} n+m-1 \\ m \end{pmatrix} t^{-m-n}
 *                      \right) R_m(y-x_*) \\
 *               & = &
 *               \sum_{m=0}^{\infty} a_m(x_i,x_*) R_m(y-x_*)
 *  \f}
 *  where
 *  \f{eqnarray*}{
 *    a_0(x_i,x_*) &=& b_0(x_i,x_*-t) \ln(t)
 *     + \sum_{n=1}^{\infty} b_n(x_i,x_*-t) (-1)^0 \begin{pmatrix} n+0 \\ 0 \end{pmatrix} t^{-0-n} \\
 *     & = & b_0(x_i,x_*-t) \ln(t) + \sum_{n=1}^{\infty} b_n(x_i,x_*-t) \frac{1}{t^n} \\
 *    a_1(x_i,x_*) &=& b_0(x_i,x_*-t) \frac{(-1)^2}{1} \frac{1}{t^1}
 *                     + \sum_{n=1}^{\infty} b_n(x_i,x_*-t) (-1)^1 \begin{pmatrix} n+1 \\ 1 \end{pmatrix} t^{-1-n} \\
 *                 &=& - b_0(x_i,x_*-t) \frac{1}{(-t)^1} + \sum_{n=1}^{\infty} - b_n(x_i,x_*-t) (n+1) \frac{1}{t^{n+1}} \\
 *                 &=& \left( \frac{1}{(-t)^1} \sum_{n=1}^{\infty} b_n(x_i,x_*-t) \begin{pmatrix} n+1 \\ 1 \end{pmatrix} \frac{1}{t^n} \right)
 *                           - \frac{b_0(x_i,x_*-t)}{(-t)^1} \\
 *    a_2(x_i,x_*) &=& b_0(x_i,x_*-t) \frac{(-1)^3}{2} \frac{1}{t^2}
 *                     + \sum_{n=1}^{\infty} b_n(x_i,x_*-t)(-1)^2 \begin{pmatrix} n+2-1 \\ 2 \end{pmatrix} t^{-2-n} \\
 *                 &=& \left( \frac{1}{(-t)^2} \sum_{n=1}^{\infty} b_n(x_i,x_*-t) \begin{pmatrix} n+1 \\ 2 \end{pmatrix} \frac{1}{t^n} \right)
 *                     - \frac{b_0(x_i,x_*-t)}{2(-t)^2} \\
 *    a_3(x_i,x_*) &=& b_0(x_i,x_*-t) \frac{(-1)^4}{3} \frac{1}{t^3}
 *                     + \sum_{n=1}^{\infty} b_n(x_i,x_*-t)(-1)^3 \begin{pmatrix} n+3-1 \\ 3 \end{pmatrix} t^{-3-n} \\
 *                 &=& \left( \frac{1}{(-t)^3} \sum_{n=1}^{\infty} b_n(x_i,x_*-t) \begin{pmatrix} n+2 \\ 3 \end{pmatrix} \frac{1}{t^n} \right)
 *                     - \frac{b_0(x_i,x_*-t)}{3(-t)^3} \\
 *                 & \vdots & \\
 *    a_m(x_i,x_*) &=& \left( \frac{1}{(-t)^m} \sum_{n=1}^{\infty} b_n(x_i,x_*-t) \begin{pmatrix} n+m-1 \\ m \end{pmatrix} \frac{1}{t^n} \right)
 *                      - \frac{b_0(x_i,x_*-t)}{m(-t)^m}
 *  \f}
 *
 *  We note that \f$\begin{pmatrix} m+n-1 \\ m \end{pmatrix} = \begin{pmatrix}  m+n-1 \\ n-1 \end{pmatrix} \f$ due to the symmetry in
 *  the formula for combinations
 *  \f[ \begin{pmatrix} m+n-1 \\ m \end{pmatrix} = \frac{(m+n-1)!}{(m+n-1-m)!m!} = \begin{pmatrix} m+n-1 \\ n-1 \end{pmatrix} \f]
 *
 *
 *  <b> Approximation </b>
 *
 *  We will appproximate the SR-translation by truncating the series to \f$p\f$ terms.
 *  This results in a \f$p\f$-by-\f$p\f$ translation matrix.  Note that the indexing
 *  for \f$p\f$ starts on one.
 *
 *  \f{eqnarray*}{
 *      \begin{bmatrix} \ln(t) & \frac{1}{t} & -\frac{1}{2t^2} & \frac{1}{3t^3} & -\frac{1}{4t^4} & \hdots
 *                      & \frac{(-1)^{p-1}}{(p-1)t^{p-1}}                                      \\
 *                      t^{-1} & -t^{-2}     & t^{-3}          & -t^{-4}        & t^{-5}          & \hdots
 *                      & (-1)^{p-1} t^{-p}                                                  \\
 *                      t^{-2} & -2t^{-3}    & 3t^{-4}         & -4t^{-5}       & 5t^{-6}         & \hdots
 *                      & (-1)^{p-1} \begin{pmatrix} p \\ 1 \end{pmatrix} t^{-p-1}             \\
 *                      t^{-3} & -3t^{-4}    & 6t^{-5}         & -10t^{-6}      & 15t^{-7}        & \hdots
 *                      & (-1)^{p-1} \begin{pmatrix} p + 1 \\ 2 \end{pmatrix} t^{-p-2}         \\
 *                      t^{-4} & -4t^{-5}    & 10t^{-6}        & -20t^{-7}      & 35t^{-8}        & \hdots
 *                      & (-1)^{p-1} \begin{pmatrix} p + 2 \\ 3 \end{pmatrix} t^{-p-3}         \\
 *                      t^{-5} & -5t^{-6}    & 15t^{-7}        & -35t^{-8}      & 70t^{-9}        & \hdots
 *                      & (-1)^{p-1} \begin{pmatrix} p + 3 \\ 4 \end{pmatrix} t^{-p-4}         \\
 *                      \vdots & \vdots   & \vdots   & \vdots    & \vdots   & \ddots
 *                      & \vdots                                                           &   \\
 *                      t^{-p+1}  & -(p-1)t^{-p}  &  \begin{pmatrix} p \\ 2 \end{pmatrix} t^{-p-1}  &
 *                      - \begin{pmatrix} p + 1 \\ 3 \end{pmatrix} t^{-p-2}                             &
 *                        \begin{pmatrix} p + 2 \\ 4 \end{pmatrix} t^{-p-3}        & \hdots
 *                      & (-1)^{p-1} \begin{pmatrix} 2p - 3 \\ p - 2 \end{pmatrix} t^{-2p+2}         \\
 *      \end{bmatrix}
 *  \f}
 *
 *
 *  We work with the transpose of this approximation in the code (getSR member function in Potential.cc).
 *  Multiplying the transpose of the matrix by the vector of coefficients for the series having the 'old' center results in the
 *  vector of coefficients for the series that has the 'new' center.  For example, consider the coefficient \f$a_0(x_i,x_*)\f$,
 *  one of the coefficients for the series with the 'new' center \f$x_*\f$, and the set of coefficients \f$b_i(x_i,x_*-t)\f$ for
 *  the series having the 'old' center \f$x_*-t\f$.
 *  The formula for \f$a_0(x_i,x_*)\f$ above is a linear combination of the first column of the translation matrix (above)
 *  with the coefficients of the series having the 'old' center.  Namely,
 *
 *  \f[ \begin{bmatrix} ln(t) & \frac{1}{t} & \frac{1}{t^2} & \hdots   & \frac{1}{t^{p-1}} \end{bmatrix}
 *      \begin{bmatrix} b_0   \\ b_1        \\ b_2          \\ \hdots  \\ b_{p-1}           \end{bmatrix}
 *  \f]
 *
 *  Therefore, it makes sense to work with the transpose matrix to perform the product as a row vector and column vector.
 *  We can see that \f$a_1\f$ is a linear combination of the second column of the translation
 *  matrix and the coefficients of the series having the 'old' center.  In the same manner we can obtain the other
 *  coefficients of the series with the 'new' center.
 *
 *
 *
 *
 *  @subsubsection far_to_far_translation_subsubsection (S|S) Translation (far-to-far translation)
 *
 *  From the previous section, we had the alternate expansion
 *
 *  \f{eqnarray*}{
 *    \frac{1}{\tilde{t}^m} &=& \sum_{n=0}^{\infty} \begin{pmatrix} -m \\ n \end{pmatrix} z^{-m-n} t^n \\
 *                          &=& \sum_{n=0}^{\infty} \begin{pmatrix} -m \\ n \end{pmatrix} \frac{1}{z^{m+n}} t^n \\
 *    \frac{1}{(y-x_*+t)^m} & = & \sum_{n=0}^{\infty} \begin{pmatrix} -m \\ n \end{pmatrix} \frac{1}{(y-x_*)^{m+n}} t^n \\
 *  \f}
 *
 *  Therefore, we have
 *
 *  \f{eqnarray*}{
 *    S_0(y-x+*+t) &=& \ln(y-x_*+t) \\
 *                 &=& \ln(y-x_*) + t \frac{1}{y-x_*} - \frac{t^2}{2} \frac{1}{(y-x_*)^2} + \frac{t^3}{3} \frac{1}{(y-x_*)^3}
 *                    - \frac{t^4}{4} \frac{1}{(y-x_*)^4} + \hdots + \frac{(-1)^nt^n}{n} \frac{1}{(y-x_*)^n} + \hdots + \\
 *                 &=& S_0(y-x_*) + t S_1(y-x_*) - \frac{t^2}{2} S_2(y-x_*) + \frac{t^3}{3} S_3(y-x_*)
 *                    - \frac{t^4}{4} S_4(y-x_*) + \hdots + \frac{(-1)^n t^n}{n} S_n(y-x_*) + \hdots + \\
 *    S_1(y-x_*+t) &=& \frac{1}{(y-x_*+t)}
 *                  = \sum_{n=0}^{\infty} \begin{pmatrix} -1 \\ n \end{pmatrix} \frac{1}{(y-x_*)^{1+n}} t^n \\
 *                 &=& \begin{pmatrix} -1 \\ 0 \end{pmatrix} \frac{1}{(y-x_*)^{1+0}} t^0
 *                   + \begin{pmatrix} -1 \\ 1 \end{pmatrix} \frac{1}{(y-x_*)^{1+1}} t^1
 *                   + \begin{pmatrix} -1 \\ 2 \end{pmatrix} \frac{1}{(y-x_*)^{1+2}} t^2
 *                   + \begin{pmatrix} -1 \\ 3 \end{pmatrix} \frac{1}{(y-x_*)^{1+3}} t^3
 *                   + \hdots + \begin{pmatrix} -1 \\ n \end{pmatrix} \frac{1}{(y-x_*)^{1+n}} t^n + \hdots \\
 *                 &=& \frac{(-1)!}{(-1-0)!0!} S_{1}(y-x_*) + \frac{(-1)!}{(-1-1)!1!} t S_2(y-x_*)
 *                     + \frac{(-1)!}{(-1-2)!2!} t^2 S_3(y-x_*) + \frac{(-1)!}{(-1-3)!3!} t^3 S_3(y-x_*)
 *                     + \hdots + \frac{(-1)!}{(-1-n)!n!} t^n S_{n+1}(y-x_*)+\hdots \\
 *                 &=& S_1(y-x_*) - tS_2(y-x_*) + t^2S_3(y-x_*) - t^3S_4(y-x_*) + \hdots + (-1)^nt^nS_{n+1}(y-x_*)+\hdots \\
 *    S_2(y-x_*+t) &=& \frac{1}{(y-x_*+t)^2}
 *                 = \sum_{n=0}^{\infty} \begin{pmatrix} -2 \\ n \end{pmatrix} \frac{1}{(y-x_*)^{2+n}} t^n \\
 *                 &=& \begin{pmatrix} -2 \\ 0 \end{pmatrix} \frac{1}{(y-x_*)^{2+0}} t^0
 *                   + \begin{pmatrix} -2 \\ 1 \end{pmatrix} \frac{1}{(y-x_*)^{2+1}} t^1
 *                   + \begin{pmatrix} -2 \\ 2 \end{pmatrix} \frac{1}{(y-x_*)^{2+2}} t^2
 *                   + \begin{pmatrix} -2 \\ 3 \end{pmatrix} \frac{1}{(y-x_*)^{2+3}} t^3
 *                   + \hdots + \begin{pmatrix} -2 \\ n \end{pmatrix} \frac{1}{(y-x_*)^{2+n}} t^n + \hdots \\
 *                 &=& \frac{(-2)!}{(-2-0)!0!} S_2(y-x_*) + \frac{(-2)!}{(-2-1)!1!} t S_3(y-x_*)
 *                     + \frac{(-2)!}{(-2-2)!2!} t^2 S_4(y-x_*) + \frac{(-2)!}{(-2-3)!3!}t^3 S_5(y-x_*)
 *                     + \hdots + \frac{(-2)!}{(-2-n)!n!} t^n S_{n+2}(y-x_*) + \hdots \\
 *                 &=& S_2(y-x_*) - 2 t S_3(y-x_*) + 3t^2S_4(y-x_*) - 4t^3S_5(y-x_*)
 *                     + \hdots + (-1)^n (n+1)t^n S_{n+2}(y-x_*)+\hdots\\
 *    S_3(y-x_*+t) &=& \frac{1}{(y-x_*+t)^3}
 *                 = \sum_{n=0}^{\infty} \begin{pmatrix} -3 \\ n \end{pmatrix} \frac{1}{(y-x_*)^{3+n}} t^n \\
 *                 &=& \begin{pmatrix} -3 \\ 0 \end{pmatrix} \frac{1}{(y-x_*)^{3+0}} t^0
 *                   + \begin{pmatrix} -3 \\ 1 \end{pmatrix} \frac{1}{(y-x_*)^{3+1}} t^1
 *                   + \begin{pmatrix} -3 \\ 2 \end{pmatrix} \frac{1}{(y-x_*)^{3+2}} t^2
 *                   + \begin{pmatrix} -3 \\ 3 \end{pmatrix} \frac{1}{(y-x_*)^{3+3}} t^3
 *                   + \hdots + \begin{pmatrix} -3 \\ n \end{pmatrix} \frac{1}{(y-x_*)^{3+n}} t^n + \hdots \\
 *                 &=& \frac{(-3)!}{(-3-0)!0!} S_3(y-x_*) + \frac{(-3)!}{(-3-1)!1!} t S_4(y-x_*)
 *                     + \frac{(-3)!}{(-3-2)!2!} t^2 S_5(y-x_*) + \frac{(-3)!}{(-3-3)!3!}t^3 S_6(y-x_*)
 *                     + \hdots + \frac{(-3)!}{(-3-n)!n!} t^n S_{n+3}(y-x_*) + \hdots \\
 *                 &=& S_3(y-x_*) - 3 t S_4(y-x_*) + 6t^2S_5(y-x_*) - 10t^3S_6(y-x_*)
 *                     + \hdots + (-1)^n \frac{(n+2)(n+1)}{2} t^n S_{n+3}(y-x_*)+\hdots\\
 *    S_4(y-x_*+t) &=& \frac{1}{(y-x_*+t)^4}
 *                 = \sum_{n=0}^{\infty} \begin{pmatrix} -4 \\ n \end{pmatrix} \frac{1}{(y-x_*)^{4+n}} t^n \\
 *                 &=& \begin{pmatrix} -4 \\ 0 \end{pmatrix} \frac{1}{(y-x_*)^{4+0}} t^0
 *                   + \begin{pmatrix} -4 \\ 1 \end{pmatrix} \frac{1}{(y-x_*)^{4+1}} t^1
 *                   + \begin{pmatrix} -4 \\ 2 \end{pmatrix} \frac{1}{(y-x_*)^{4+2}} t^2
 *                   + \begin{pmatrix} -4 \\ 3 \end{pmatrix} \frac{1}{(y-x_*)^{4+3}} t^3
 *                   + \hdots + \begin{pmatrix} -4 \\ n \end{pmatrix} \frac{1}{(y-x_*)^{4+n}} t^n + \hdots \\
 *                 &=& \frac{(-4)!}{(-4-0)!0!} S_4(y-x_*) + \frac{(-4)!}{(-4-1)!1!} t S_5(y-x_*)
 *                     + \frac{(-4)!}{(-4-2)!2!} t^2 S_6(y-x_*) + \frac{(-4)!}{(-4-3)!3!}t^3 S_7(y-x_*)
 *                     + \hdots + \frac{(-4)!}{(-4-n)!n!} t^n S_{n+4}(y-x_*) + \hdots \\
 *                 &=& S_4(y-x_*) - 4 t S_5(y-x_*) + 10t^2S_6(y-x_*) - 20t^3S_7(y-x_*)
 *                     + \hdots + (-1)^n \frac{(n+3)(n+2)(n+1)}{3!} t^n S_{n+4}(y-x_*)+\hdots\\
 *  \f}
 *
 *    The transformation matrix is
 *  \f[ \begin{bmatrix} S_0(y-x_*+t) \\ S_1(y-x_*+t) \\ S_2(y-x_*+t) \\ S_3(y-x_*+t) \\ S_4(y-x_*+t) \\
 *                      S_5(y-x_*+t) \\ S_6(y-x_*+t) \\  \vdots \\ S_n(y-x_*) \\
 *                      S_{n+1}(y-x_*) \\ S_{n+2}(y-x_*) \\ S_{n+3}(y-x_*) \\ S_{n+4}(y-x_*) \\ \vdots
 *      \end{bmatrix}
 *      =
 *      \left[
 *      \begin{array}{cccccccccccc}
 *                      1 & t      & - \frac{t^2}{2}  & \frac{t^3}{3}  & \hdots    &
 *                      & \frac{(-1)^n t^n}{n}  & & \hdots
 *                      & \hdots
 *                      & \hdots & \hdots \\
 *                      0 & 1          & -t     & t^{2}            &                & \hdots    &
 *                      & (-1)^n t^{n}          & \hdots
 *                      & \hdots
 *                      & \hdots & \hdots \\
 *                      0 & 0          & 1      & -2t              & 3t^{2}         &           & \hdots
 *                      &                 &(-1)^n \begin{pmatrix} n+1 \\ 1 \end{pmatrix} t^n
 *                      & \hdots
 *                      & \hdots & \hdots \\
 *                      0 & 0          & 0      & 1                & -3t            & 6t^{2}    &
 *                      & \hdots          &
 *                      & (-1)^n \begin{pmatrix} n+2 \\ 2 \end{pmatrix} t^{n}
 *                      & \hdots & \hdots \\
 *                      0 & 0          & 0      & 0                & 1              & -4t       & 10t^2
 *                      &                 & \hdots
 *                      &
 *                      & (-1)^n \begin{pmatrix} n+3 \\ 3 \end{pmatrix} t^n & \hdots \\
 *                      \vdots & \vdots & \vdots & \ddots   & \ddots   & \ddots    & \ddots
 *                      & \ddots &
 *                      & \hdots & \vdots & \ddots
 *      \end{array}
 *      \right]
 *      \begin{bmatrix} S_0(y-x_*) \\ S_1(y-x_*) \\ S_2(y-x_*) \\ S_3(y-x_*) \\ S_4(y-x_*) \\ S_5(y-x_*) \\ S_6(y-x_*) \\
 *                      \vdots \\ S_n(y-x_*) \\ S_{n+1}(y-x_*) \\ S_{n+2}(y-x_*) \\ S_{n+3}(y-x_*) \\ S_{n+4}(y-x_*) \\ \vdots
 *      \end{bmatrix}
 *  \f]
 *
 *  The formula for the coefficients \f$b_m(x_i,x_*)\f$ in the translation
 *  \f[ \sum_{m=0}^{\infty} b_m(x_i,x_*-t) S_m(y-x_*+t)
 *       =  \sum_{m=0}^{\infty} b_m(x_i,x_*)S_m(y-x_*)
 *       \f]
 *  is obtained by multiplying the equations above by the corresponding
 *  \f$b_m(x_i,x_*-t)\f$ and adding the equations resulting in the coefficients
 *  \f{eqnarray*}{
 *    b_0(x_i,x_*) &=& b_0(x_i,x_*-t) \\
 *    b_1(x_i,x_*) &=& b_0(x_i,x_*-t)t + b_1(x_i,x_*-t) \\
 *    b_2(x_i,x_*) &=& -b_0(x_i,x_*-t)\frac{t^2}{2} - b_1(x_i,x_*-t)t + b_2(x_i,x_*-t) \\
 *    b_3(x_i,x_*) &=& b_0(x_i,x_*-t)\frac{t^3}{3} + b_1(x_i,x_*-t)t^2 - b_2(x_i,x_*-t) 2t + b_3(x_i,x_*-t) \\
 *                 & \vdots & \\
 *    b_n(x_i,x_*) &=& -b_0(x_i,x_*-t) \frac{(-t)^n}{n}
 *                     + \sum_{m=1}^{n} b_m(x_i,x_*-t) \begin{pmatrix} n-1 \\ m -1 \end{pmatrix} (-t)^{n-m}
 *  \f}
 *
 *
 *
 *  <b> Approximation </b>
 *
 *  We will also appproximate the translation.  Truncating the series to \f$p\f$ terms, results
 *  in a \f$p\f$-by-\f$p\f$ matrix.  Note that the indexing for \f$p\f$ starts on one.
 *
 *  \f{eqnarray*}{
 *        \left[
 *      \begin{array}{ccccccccc}
 *                      1 & t      & - \frac{t^2}{2}  & \frac{t^3}{3}  & \frac{t^4}{4} & \frac{t^5}{5} & \frac{t^6}{6} & \hdots
 *                      & \frac{(-1)^{p} t^{p-1}}{p-1}  \\
 *                      0 & 1      & -t               & t^{2}          & -t^{3}        & t^{4}         & t^{5}         & \hdots
 *                      & (-1)^{p} t^{p-2}    \\
 *                      0 & 0      & 1                & -2t            & 3t^{2}        & 4t^{3}        & 5t^{4}        & \hdots
 *                      & \begin{pmatrix} p-2 \\ 1 \end{pmatrix} (-1)^{p-1} t^{p-3} \\
 *                      0 & 0      & 0                & 1              & -3t           & 6t^{2}        & 10t^{3}       & \hdots
 *                      & \begin{pmatrix} p-2 \\ 2 \end{pmatrix} (-1)^{p-2} t^{p-4} \\
 *                      0 & 0      & 0                & 0              & 1             & -4t           & 10t^2         & \hdots
 *                      & \begin{pmatrix} p-2 \\ 3 \end{pmatrix} (-1)^{p-3} t^{p-5} \\
 *                      0 & 0 & 0 & 0 & 0 & 1 & \ddots & \ddots & \vdots \\
 *                      0 & 0 & 0 & 0 & 0 & 0 & 1      & \ddots & \vdots \\
 *                      0 & 0 & 0 & 0 & 0 & 0 & 0      & \ddots & \vdots \\
 *                      0 & \hdots  &   & \hdots  &   &   & \hdots & 0      & 1
 *      \end{array}
 *      \right]
 *  \f}
 *
 *  The transpose (shown below) is used in the source code
 *
 *  \f{eqnarray*}{
 *        \left[
 *      \begin{array}{ccccccccc}
 *                      1                              & 0                  & 0                                                           & 0
 *                      & 0                                                                      & 0             & 0             & \hdots
 *                      &                              0 \\
 *                      t                              & 1                  & 0                                                           & 0
 *                      & 0                                                                      & 0             & 0             & \hdots
 *                      &                              0 \\
 *                      - \frac{t^2}{2}                & -t                 & 1                                                           & 0
 *                      & 0                                                                      & 0             & 0             & \hdots
 *                      &                              0 \\
 *                      \frac{t^3}{3}                  & t^{2}              & -2t                                                         & 1
 *                      & 0                                                                      & 0             & 0             & \hdots
 *                      &                              0 \\
 *                      \frac{t^4}{4}                  & -t^{3}             & 3t^{2}                                                      & -3t
 *                      & 1                                                                      & 0             & 0             & \hdots
 *                      &                              0 \\
 *                      \frac{t^5}{5}                  & t^{4}              & 4t^{3}                                                      & 6t^{2}
 *                      & -4t                                                                    & 1             & 0             & \ddots
 *                      &                              0 \\
 *                      \frac{t^6}{6}                  & t^{5}              & 5t^{4}                                                      & 10t^{3}
 *                      & 10t^2                                                                  & \ddots        & 1             & \ddots
 *                      &                              0 \\
 *                      \vdots                         & \vdots             & \vdots                                                      & \vdots
 *                      & \vdots                                                                 & \ddots        & \ddots         & \ddots
 *                      &                              0 \\
 *                      \frac{(-1)^{p} t^{p-1}}{p-1}   & (-1)^{p} t^{p-2}   & \begin{pmatrix} p-2 \\ 1 \end{pmatrix} (-1)^{p-1} t^{p-3}   & \begin{pmatrix} p-2 \\ 2 \end{pmatrix} (-1)^{p-2} t^{p-4}
 *                      & \begin{pmatrix} p-2 \\ 3 \end{pmatrix} (-1)^{p-3} t^{p-5}              & \hdots        & \hdots        & \hdots
 *                      &                              1 \\
 *      \end{array}
 *      \right]
 *  \f}
 *
 *  In the code, the elements of the transpose are defined recursively with respect to a given row.
 *  We can see the general pattern in the last row of the transpose matrix.  We traverse a row working from the main diagonal
 *  and moving left.  We can see that the power of \f$t\f$ increases by one as we move from element to element.
 *  The sign alternates from negative to positive.  To determine the behavior of the coefficients while progressing left
 *  from the main diagonal element by element, consider the \f$6^{th}\f$ and \f$7^{th}\f$ rows.  Following the progression
 *  from left to right we may recognize if as Pascal's Triangle or the progression of combinations as shown in the formulas for
 *  the coefficients in bottom row of the matrix.
 *
 *
 *  Take for example the \f$6^{th}\f$ row.
 *
 *           The first coefficient from the left of the main diagonal 1 is 4.
 *           This can be obtained from the formula at the bottom of this column with p = 6
 *           \f[ \begin{pmatrix} 6 - 2 \\ 3 \end{pmatrix} = \begin{pmatrix} 4 \\ 3 \end{pmatrix} = \frac{ 4! }{ 3! (4-3)!}
 *                                               = \frac{ 4 \cdot 3 \cdot 2 \cdot 1 }{ 3 \cdot 2 \cdot 1 \cdot 1}
 *                                               = \frac{ 4 }{ 1 }
 *           \f]
 *
 *          We could obtain this from multiplying the main diagonal value 1.0 by j-1 = 5- 1 = 4 and dividing the result by (i-j) = 6-5 = 1
 *
 *          The second coefficient from the left is 6.
 *          From the formula at the bottom of this column with p = 6
 *           \f[ \begin{pmatrix} 6 - 2 \\ 2 \end{pmatrix} = \begin{pmatrix} 4 \\ 2 \end{pmatrix} = \frac{ 4! }{ 2! (4-2)!}
 *                                               = \frac{ 4 \cdot 3 \cdot 2 \cdot 1 }{ 2 \cdot 1 \cdot 2 \cdot 1}
 *                                               = \frac{ 4 \cdot 3 }{ 2 \cdot 1 }
 *                                               = \frac{4}{1} \frac{3}{2}
 *           \f]
 *
 *          We could obtain this from multiplying the previous value 4/1 by j - 1 = 4 - 1 = 3 and dividing the
 *          result by (i-j) = 6 - 4 = 2.
 *
 *          We see this pattern and method in the source code.
 *
 *
 *
 *  @subsection algorithm_subsection MLFMM Algorithm
 *
 *  The Multilevel Fast Multiple Method algorithm can be separated into four steps.
 *
 *  @subsubsection quadtree_partition_subsubsection Quadtree Partition
 *
 *  First we enclose a square around the set of target and source particles.
 *  This code will be based on the assumption that the particles are uniformly distributed.
 *  This assumption allows us to easily divide (or partition) the square up and have the same number of
 *  particles in each subdivision (or level of partition).
 *
 *  We call the original square level zero (\f$l = 0\f$) of the partition and use \f$l\f$ to represent the level.
 *  The subdivision technique is to first take the original square and
 *  connect lines from the opposing edges' midpoints.  This divides the original square into
 *  four smaller squares with the same size (sides have the same length).  This is known as a quadtree.
 *
 *  Due to the uniform distribution
 *  of the particles, each subdivision also has approximately the same number of particles.
 *  We call this first level of subdivision level 1 (\f$l = 1\f$).  We perform the same
 *  technique on each of the four square of \f$l=1\f$.  This results in level \f$l=2\f$
 *  having \f$4 * 4 = 16 \f$ squares.
 *
 *  So far, we see that (i) level zero, \f$l=0\f$, has \f$1 = 4^0\f$ cells (squares) in its partition,
 *  (ii) level one, \f$l=1\f$, has \f$4^1 = 4\f$ square (or cells) making up its partition and \f$l=2\f$, (iii) level two, has
 *  \f$4^2 = 16\f$ cells making up its partition.  Continuing in this manner, level \f$l\f$, \f$l = l \f$ will have
 *  \f$4^l\f$ cells in its partition of the original square.  The highest (or finest) level of refinement will be denoted
 *  level \f$l=L\f$.
 *
 *  We will uniquely identify each square or cell by its ordered pair \f$(n,l\f$ where \f$l\f$ is the level that the cell lives
 *  on and \f$n\f$ is its index number on that level.  We will talk about
 *  the indexing next.
 *
 *  To perform the MLFMM we partition each level into four disjoint sets or domains (see the image below):
 *  (i) the box \f$E_1(n,l) = (n,l)\f$ at level \f$l\f$ having index number \f$n\f$
 *  (ii) the nearest neighbors of the box \f$E_2(n,l)\f$ (which includes the box itself)
 *  (iii) the entire domain (original square) without the nearest neighbors \f$E_3(n,l) = E_1(0,0) \setminus E_2(n,l)\f$ and
 *  (iv) the interaction list \f$E_4(n,l)\f$ made up of the set of parents of the nearest neighbors without the nearest neighbors.
 *
 *  @image html lists.png
 *
 *  Each domain contains particles.  We will consider each particle a source \f$x_i\f$ having a potentials \f$u_i\f$.
 *  Recall that \f$u_i \Phi(y,x_i)\f$ is used to calculate the action of a source \f$x_i\f$ on a target particle \f$y\f$.
 *  The potentials due to sources in each domain are
 *  \f{eqnarray*}{
 *    \Phi_1^{(n,l)} (y) \colon = \sum_{x_i \in E_1(n,l)} u_i \Phi(y,x_i) \\
 *    \Phi_2^{(n,l)} (y) \colon = \sum_{x_i \in E_2(n,l)} u_i \Phi(y,x_i) \\
 *    \Phi_3^{(n,l)} (y) \colon = \sum_{x_i \in E_3(n,l)} u_i \Phi(y,x_i) \\
 *    \Phi_4^{(n,l)} (y) \colon = \sum_{x_i \in E_4(n,l)} u_i \Phi(y,x_i) \\
 *  \f}
 *
 *
 *
 *  @subsubsection far_field_expansion_subsubsection Far-Field Expansions
 *
 *  We start at the highest level of refinement \f$l=L\f$.  Let the particle \f$y\f$ be an arbitrary target
 *  located somewhere in the original square (E_1(0,0) - cell with index \f$n=0\f$ at refinement level \f$l=0\f$).
 *  Consider any cell \f$E_1(n,L)\f$ in the refinement \f$l=L\f$.  For each source \f$x_i \in E_1(n,L)\f$, we write \f$\Phi(y,x_i)\f$
 *  as a Taylor series S-Expansion (far-field expansion).  We do this for every cell \f$E_1(n.L)\f$.  Therefore, we
 *  have two loops: (i) an outside loop over the cells of refinement level \f$L\f$ and (ii) an inside loop over the sources in that cell.
 *
 *  Recall that for the R and S-expansions there is a center, \f$x_*\f$ for example.  For these expansions, we make the center
 *  of the cell \f$x_c^{(n,L)}\f$ the center of the expansion.
 *
 *  Note that for each \f$x_i \in E_1(n,l)\f$, the far-field expansion has the same components \f$S_m(y-x_c) = \frac{1}{y-c_x}^m\f$.
 *  Therefore, a single far-field series can be formed
 *  \f[ \Phi_1^{(n,L)}(y) = \sum_{x_i \in E_1(n,L)} u_i \Phi(y,x_i)
 *                        = \sum_{m=0}^{\infty} C_m^{(n,L)} S_m(y-x_c)
 *  \f]
 *
 *  where
 *  \f[ C^{(n,L)}_m = \sum_{x_i \in E_1(n,L)} u_i b_m(x_i,x_c^{(n,L)}) \f]
 *
 *  @subsubsection upward_pass_subsubsection Upward Pass
 *
 *  We now have far-field expansions for each source in each cell at the highest refinement level \f$l=L\f$.
 *  Our next step is translate these expansions up the refinement levels.
 *
 *  We start at the second highest refinement level \f$l = L-1\f$.  Each cell at this refinement level is
 *  the parent of four children that are cells of refinement level \f$l = L\f$.  We take an arbitrary
 *  cell \f$(n,L-1)\f$ at this refinement level \f$l = L-1\f$ and translate the far-field expansions from the sources
 *  in this parent cell to the center of the cell \f$x_c^{(n,L-1)}\f$.  That is, for each source \f$x_i\f$ in each of
 *  the children cells \f$(child_i(n,L-1),L)\f$ for \f$1 \leq i \leq 4\f$ we translate the far-field expansion
 *  from the center of the child cell \f$x_c^{(child_i(n,L-1),L)}\f$ to the center of the parent cell \f$x_c^{(n,L-1)}\f$.
 *
 *  If, for example, we have 4 source particles in each child cell, then we will be translating \f$4*4=16\f$ expansions
 *  to the center of the parent cell.  However, for each child cell, we have written the four expansions of that cell
 *  as a single series.  Therefore, we are really translating only four expansions, one for each child cell.
 *
 *  We now have all the four far-field expansions corresponding to each of the child cells, centered at the center of the parent cell.
 *  Therefore, each expansion now has the same powers \f$S_m(y-x_c^{(n,L-1)})\f$ and we can write all the expansions as a single series
 *  \f[ \Phi_1^{(n,L-1)}(y)
 *          = \sum_{n=0}^{\infty} C_m^{(n,L-1)} S_m(y-x_c^{(n,L-1)})
 *  \f]
 *  where the coefficient \f$C_m^{(n,L-1)}\f$ depends on the far-field to far-field translation from child centers to the parent center.
 *
 *  Using our work on S|S translations, let's look at the coefficients \f$C_m^{(n,L-1)}\f$ for this first step of translations.
 *  We can write down the exact formula for the \f$C_m^{(n,L-1)}\f$.  From our work on \f$S|S\f$ translations, we have the correspondences
 *  \f[ y-x_* \leftrightarrow y-x_c^{(n',L)} \quad \mbox{and} \quad y-x_*+t = y-(x_*-t) \leftrightarrow y-x_c^{(n,L-1)} \f]
 *  Then \f$y-x_c^{(n',L)} = y-x_c^{(n,L-1)} + t\f$ implies \f$t = x_c^{(n,L-1)} - x_c^{(n',L)}\f$, and
 *
 *  \f{eqnarray*}{
 *    C_0^{(n,L-1)} &=& \sum_{n' \in child(n,L-1)} C_0^{(n',L)} \\
 *    C_m^{(n,L-1)} &=& \sum_{n' \in child(n,L-1)}
 *                      \left(  -C_0^{(n',L)} \frac{(x_c^{(n',L)}-x_c^{(n,L-1)})^m}{m}
 *                                + \sum_{k=1}^{m} C_k^{(n',L)} \begin{pmatrix} m-1 \\ k-1 \end{pmatrix} (x_c^{(n',L)}-x_c^{(n,L-1)})^{m-k}
 *                      \right)
 *  \f}
 *
 *  We repeat this process for each cell (n,L-2) of the next refinement level \f$l=L-2\f$ up by translating the expansion
 *  each of the four children \f$n' \in child(n,L-2)\f$ has to the center \f$x_c^{(n,L-2)}\f$ of the parent cell \f$(n,L-2)\f$.
 *  We continue the process until we finish at level \f$l=2\f$.
 *
 *  @subsubsection downward_pass_subsubsection Downward Pass
 *
 *  In this part of the algorithm we are working from the refinement level two (\f$l = 2\f$) down to refinement level
 *  \f$l = L\f$.
 *
 *  Our first step will be to translate the far-field expansions that we worked on in the upward pass
 *  to R-expansions using an \f$(S|R)\f$ translation.  For each cell \f$(n,2)\f$ of level \f$l=2\f$ we translate the far-field
 *  expansions of the cells in its interactive list to its center \f$x_c^{(n,2)}\f$ using an \f$(S|R)\f$ translation.
 *
 *  Once \f$l=2\f$ is complete, we step down to \f$l=3\f$ and perform the same method for the cells \f$(n,3)\f$ in level \f$l=3\f$.
 *
 *
 *  <b> Local Expansions </b>
 *
 *  Recall that \f$\Phi_4^{(n,l)}(y)\f$ was defined as
 *  \f[ \Phi_4^{(n,l)}(y) = \sum_{x_i \in E_4(n,l)} u_i \Phi(y,x_i) \f]
 *  where \f$E_4(n,l)\f$ was called the interaction list and was made up of the cells that were children
 *  of the nearest neighbors of the parent of the cell of interest \f$(n,l)\f$.  But the list did not include the nearest neighbors
 *  of cell \f$(n,l)\f$ or cell \f$(n,l)\f$ itself (see the figure below).
 *
 *  Note that level \f$l=0\f$ and level \f$l=1\f$ do not have an interaction list (only nearest neighbors), and that
 *  is why we have started our work at level \f$l=2\f$.
 *  Level \f$l=2\f$ is the first level to have an interaction list, where approximations (and not direct
 *  calculations of the logarithmic function) and the fast multipole method can be done.  Further, level
 *  \f$l=2\f$ has the unique property that for any cell \f$(n,2)\f$ the entire domain \f$E_1(0,0)\f$ consists of
 *  that cells nearest neighbors and its interaction list.  This can not be said for higher refinement levels (see the
 *  figure below).
 *
 *  We can say this symbolically by recalling the definition
 *  \f[ \Phi_3^{(n,l)} (y) \colon = \sum_{x_i \in E_3(n,l)} u_i \Phi(y,x_i) \f]
 *  where \f$E_3(n,l)\f$ is the set of all cells in the entire domain \f$E(0,0)\f$ that are not nearest neighbors
 *  of \f$(n,l)\f$.
 *
 *  Therefore, for level \f$l=2\f$
 *  \f[ \Phi_3^{(n,2)} (y) = \Phi_4^{(n,2)}(y). \f]
 *  We want to be able to state \f$\Phi_3^{(n,l)}(y)\f$ for any level \f$l\f$.  We make the argument for the method
 *  to determine \f$\Phi_3^{(n,l)}(y)\f$ for any level \f$l\f$.
 *
 *  To help illustrate the method, shown below are figures of levels \f$l=2\f$ to \f$l=5\f$.  To help establish the
 *  method, the cell of interest
 *  (dark blue) at levels \f$l=3\f$ to \f$l=5\f$ are children, grandchildren or great grandchildren of the cell of
 *  interest at level \f$l=2\f$.
 *
 *  @image html levels2to5.png
 *
 *  For \f$l=2\f$ and each cell \f$(n,2)\f$ of level \f$l=2\f$, we translate the far-field expansion of each cell in the interaction
 *  list of cell \f$(n,2)\f$.  Therefore, each cell \f$(n,2)\f$ of level \f$l=2\f$ has the potential calculations of its interaction list
 *  (each calculation is an series approximation) centered at \f$x_c^{(n,2)}\f$, the center of cell \f$(n,2)\f$.  The interaction
 *  list of a cell consists of at most \f$(3x 3 x 4) - (3 x 3 x 1) = (6 x 6) - (3 x 3) = 36 - 9 = 27\f$ cells.  Therefore, we may
 *  be performing 27 translations to the cell center \f$x_c^{(n,2)}\f$ in this step.  All of the translated series now have the same
 *  center and can be written as a single expansion.  We perform this procedure for every cell \f$(n,2)\f$
 *  at this level \f$l=2\f$; that is, for \f$4^l = 4^2 = 16\f$ cells.
 *
 *   *  \f[ \Phi_3^{(n,2)} (y) = \Phi_4^{(n,2)}(y) \f]
 *
 *  We go down to the next level \f$l=3\f$.  We proceed in the same manner and translate all far-field expansions for each cell in
 *  the interaction list of \f$(n,3)\f$ to its center \f$x_c^{(n,3)}\f$.  Level \f$l=3\f$ is the first level (working down from \f$l=0\f$)
 *  that has cells beyond its interaction list.  The cells are in cyan (light cyan) in the second figure from the left above.
 *  We can see that these are the cells of the interaction list of its parent \f$(Parent(n,3),2)\f$ (seen the first figure, far left).
 *
 *  Now, the parent \f$(Parent(n,3),2)\f$ of cell \f$(n,3)\f$ has already collected the expansions of these cells of its interaction
 *  list in the previous step at \f$l=2\f$.  Therefore, to obtain this information (the single expansion) for cell \f$(n,3)\f$ at its
 *  center \f$x_c^{(n,3)}\f$, we need only perform a near-field to near-field translation from parent cell \f$(Parent(n,3),2)\f$
 *  to cell \f$(n,3)\f$ of the single expansion.
 *
 *  We can now write all of the expansions that we have collected into a single series and the result is \f$\Phi_3^{(n,3)}(y)\f$.
 *
 *  \f[ \Phi_3^{(n,3)}(y) = \Phi_4^{(n,3)}(y) + \Phi_4^{(Parent(n,3),2)}(y) \f]
 *
 *  We continue to work down the levels of refinement to \f$l=4\f$.  We translate the expansions of the cells in the interaction
 *  list of \f$(n,4)\f$ to its center \f$x_c^{(n,4)}\f$.  From the previous step at \f$l=3\f$, the parent \f$(Parent(n,4),3)\f$
 *  of \f$(n,4)\f$ has a series \f$\Phi_3^{(n,3)}(y)\f$ that is a collection of its interaction list expansions as well as its
 *  parent's interaction list expansions.  Therefore, we need only translate this series to the center of cell \f$(n,4)\f$ and we
 *  can form a single series that is \f$\Phi_3^{(n,4)}(y)\f$.
 *
 *   \f[ \Phi_3^{(n,4)}(y) = \Phi_4^{(n,4)}(y) + \Phi_3^{(Parent(n,4),3)}(y) \f]
 *
 *  We continue this process until we reach the last level \f$l=L\f$.  We see that for \f$l \geq 4\f$, the parent has recursively
 *  collected all of the information (series) from the cells outside the interaction list of the child cell (outside the parent's
 *  nearest neighbors - and this is \f$\Phi)_3\f$).  Therefore the following formula holds for \f$l \geq 4\f$
 *
 *   \f[ \Phi_3^{(n,l)}(y) = \Phi_4^{(n,l)}(y) + \Phi_3^{(Parent(n,l),l-1)}(y) \f]
 *
 *  <b> Translation Formulas</b>
 *
 *  Let's construct the formulas for the translations that take place in the downard pass.
 *  The first translation is the \f$S|R\f$ translation of the far-field expansion at each cell in the
 *  interaction list of the cell \f$(n,l)\f$ to its center \f$x_c^{(n,l)}\f$.  Recall that before the
 *  translation, the far-field expansion for each cell in the interaction list is centered at
 *  \f$x_c^{(InteractionList(n,l),l)}\f$.
 *
 *  Recall the far field expansion for each cell \f$(m,l)\f$ at level \f$l\f$
 *  \f[ \Phi_1^{(m,l)}(y) = \sum_{k=0}^{\infty} C_k^{(n,l)} S_k(y-x_c^{(n,l)}) \f]
 *  We start at the first refinement level \f$l=2\f$.  For a given cell \f$(n,l)\f$ at level \f$l=2\f$,
 *  we wish to translate each of these expansions
 *  for the cells \f$(m,l)\f$ that are in the interaction list \f$\left\{ (InteractionList(n,l),l) \right\}\f$
 *  of cell \f$(n,l)\f$ to the center \f$x_c^{(n,l)}\f$ of the cell.
 *
 *  First, let's state the potentials due to sources in the interaction list \f$E_4(n,l)\f$.
 *  \f[ \Phi_4^{(n,l)}(y)
 *       =
 *        \sum_{m \in E_4(n,l)} u_m \Phi_1^{(m,l)}(y) \f]
 *  We wish to translate each expansion \f$\Phi_1^{(m,l)}(y)\f$ from the center \f$x_c^{(m,l)}\f$
 *  to \f$x_c^{(n,l)}\f$ as well as change from a far-field \f$S\f$ to a near-field \f$R\f$ expansion.
 *  Using the \f$(S|R)\f$-translation with \f$y-x_*+t=y-x_c^{(m,l)}\f$ and \f$y-x_*=y-x_c^{(n,l)}\f$
 *  implying \f$y-x_c^{(m,l)} = y-x_c^{(n,l)}+t\f$ and \f$t = x_c^{(n,l)} - x_c^{(m,l)}\f$ we have
 *  \f{eqnarray*}{
 *    \Phi_1^{(m,l)}(y) &=& \sum_{k=0}^{\infty} C_k^{(m,l)} S_k(y-x_c^{(m,l)}) \\
 *      & = &
 *        \sum_{k=0}^{\infty} \hat{D}_k^{(m,l)} R_k(y-x_c^{(n,l)})
 *  \f}
 *
 *  where
 *
 *  \f{eqnarray*}{
 *    \hat{D}_0^{(m,l)} & = & C_0^{(m,l)} ln(x_c^{(n,l)} - x_c^{(m,l)})
 *                      + \sum_{j=1}^{\infty} C_j^{(m,l)} \frac{1}{(x_c^{(n,l)} - x_c^{(m,l)})} \\
 *    \hat{D}_k^{(m,l)} & = & - \frac{C_0^{(m,l)}}{k  (x_c^{(m,l)} - x_c^{(n,l)})^k}
 *                      + \frac{1}{(x_c^{(m,l)} - x_c^{(n,l)})^k} \sum_{j=1}^{\infty} C_j^{(m,l)} \begin{pmatrix} j+k-1 \\ k \end{pmatrix}
 *                                                           \frac{1}{(x_c^{(n,l)} - x_c^{(m,l)})^j}
 *  \f}
 *
 *  We perform the translation for the far field expansion of every cell in the interaction list.  The resulting expansions can then
 *  be written as a single series
 *  \f[  \Phi_4^{(n,l)}(y) = \sum_{k=0}^{\infty} D_k^{(n,l)} R_k(y-x_c^{(n,l)}) \f]
 *
 *  where
 *  \f{eqnarray*}{
 *    D_0^{(n,l)} &=& \sum_{m \in E_4(n,l)} u_m \hat{D}_0^{(m,l)} \\
 *    D_k^{(n,l)} &=& \sum_{m \in E_4(n,l)} u_m \hat{D}_k^{(m,l)}
 *  \f}
 *
 *  The last translation in the Downward Pass is the \f$(R|R)\f$-translation from the parent to the child
 *  of the far-field expansions that are outside the parent's nearest neighbors in \f$E_3(Parent(n,l),l-1)\f$
 *  (the parents interaction list as well as everything in the domain beyond the interaction list).  Note that
 *  the nearest neighbors of the parent contain the cell of interest, its nearest neighbors and its interaction list.
 *  Therefore, the last translation in the Downward Pass takes into account the far field expansions that are in
 *  cells that are beyond cell \f$(n,l)\f$'s interaction list.
 *
 *  Start with level \f$l=2\f$.  The interaction list of cell \f$(n,l)\f$ is all that is in the domain \f$E_1(0,0)\f$
 *  beyond the nearest neighbors of the cell.  We have just taken into account the interaction list in the previous
 *  step.  We are done for \f$l=2\f$ and do not have to consider the parent's interaction list or beyond.
 *
 *  We move on to level \f$l=3\f$.  We have taken into account the interaction list of cell \f$(n,l)\f$ at this level.
 *  There is a parent interaction list for each cell at \f$l=3\f$ as well.  The far-field expansions for the cells in this interaction
 *  list have also been translated already to the parent cell and are written as a single expansion (as explained above).
 *  It is at this point that we do our first \f$(R|R)\f$-translation.  We translate  to the center \f$x_c^{(n,l)}\f$
 *  \f{eqnarray*}{
 *    \Phi_4^{(Parent(n,l),l-1))}(y) &=& \sum_{k=0}^{\infty} D_k^{(Parent(n,l),l-1)} R_k(y-x_c^{(Parent(n,l),l-1)}) \\
 *      \sum_{k=0}^{\infty} E_k^{(n,l)} R_k(y-x_c^{(n,l)})
 *  \f}
 *  Here, \f$y-x_*+t = y - x_c^{(Parent(n,l),l-1)}\f$ and \f$y-x_* = y - x_c^{(n,l)}\f$ implies \f$t = x_c^{(n,l)} - x_c^{(Parent(n,l),l-1)}\f$.
 *  Therefore
 *  \f{eqnarray*}{
 *    E_k^{(n,l)} &=& \sum_{j=k}^{\infty} D_j^{(Parent(n,l),l-1)} \begin{pmatrix} j \\ k \end{pmatrix} (x_c^{(n,l)} - x_c^{(Parent(n,l),l-1)})^{j-k}
 *  \f}
 *
 *  For each cell \f$(n,l)\f$ at level \f$l=3\f$, we now have a way to calculate \f$\Phi_3^{(n,l)}(y)\f$
 *  \f{eqnarray*}{
 *   \Phi_3^{(n,l)}(y) &=& \Phi_4^{(n,l)}(y) + \Phi_3^{(Parent(n,l),l-1)}(y) \\
 *                     &=& \Phi_4^{(n,l)}(y) + \Phi_4^{(Parent(n,l),l-1)}(y) \\
 *                     &=& \sum_{k=0}^{\infty} \left( D_k^{(n,l)} + E_k^{(n,l)} \right) R_k(y-x_c^{(n,l)})
 *  \f}
 *  where the second equation above is due to \f$\Phi_3^{(n,l)}(y) = \Phi_4^{(n,l)}(y)\f$ for \f$l=2\f$.  Note that \f$E_k^{(n,l)}\f$
 *  is dependent on the \f$D_k^{(Parent(n,l),l-1)}\f$'s determined in the previous step of the downward pass for the refinement \f$l-1\f$.
 *
 *  Now that we have \f$\Phi_3^{(n,l)}(y)\f$ for level \f$l=3\f$, we can move on to level \f$l=4\f$.  First we carry out the
 *  \f$(S|R)\f$-translation for the interaction list and then the translation for \f$\Phi_3^{(Parent(n,l),l-1)}(y)\f$.  In the
 *  cases \f$l \geq 4\f$, we have a different situation since \f$\Phi_3^{(Parent(n,l),l-1)}(y) \neq \Phi_4^{(Parent(n,l),l-1)}(y)\f$
 *  (as was the case for \f$l=3\f$).
 *  \f{eqnarray*}{
 *    \Phi_3^{(n,l)}(y) &=& \Phi_4^{(n,l)}(y) + \Phi_3^{(Parent(n,l),l-1)}(y) \\
 *                      &=& \sum_{k=0}^{\infty} D_k^{(n,l)} R_k(y-x_c^{(n,l)})
 *                          + \sum_{k=0}^{\infty} \left( D_k^{(Parent(n,l),l-1)} + E_k^{(Parent(n,l),l-1)} \right)
 *                                                   R_k(y-x_c^{(Parent(n,l),l-1)}) \\
 *                      &=& \sum_{k=0}^{\infty} D_k^{(n,l)} R_k(y-x_c^{(n,l)})
 *                          + \sum_{k=0}^{\infty} \hat{F}_k^{(Parent(n,l),l-1)}
 *                                                   R_k(y-x_c^{(Parent(n,l),l-1)}) \\
 *  \f}
 *  where \f$\hat{F}_k^{(Parent(n,l),l-1)} = D_k^{(Parent(n,l),l-1)} + E_k^{(Parent(n,l),l-1)} \f$.
 *
 *  We translate the second series to \f$x_c^{(n,l)}\f$.  As with the \f$(R|R)\f$-translation for
 *  \f$\Phi_4^{(Parent(n,l),l-1)}(y)\f$ for level \f$l=3\f$
 *  we have
 *  \f{eqnarray*}{
 *  \Phi_3^{(Parent(n,l),l-1)}(y) &=& \sum_{k=0}^{\infty} \hat{F}_k^{(Parent(n,l),l-1)}
 *                                                   R_k(y-x_c^{(Parent(n,l),l-1)}) \\
 *                                &=& \sum_{k=0}^{\infty} F_k^{(n,l)}
 *                                                   R_k(y-x_c^{(n,l)})
 *  \f}
 *  where
 *  \f[ F_k^{(n,l)} = \sum_{j=k}^{\infty} \hat{F}_j^{(Parent(n,l),l-1)} \begin{pmatrix} j \\ k \end{pmatrix}
 *                                               (x_c^{(n,l)} - x_c^{(Parent(n,l),l-1)})^{j-k}
 *  \f]
 *
 *  We repeat the process for \f$l \geq 4\f$ until we reach \f$l=L\f$.
 *
 *  <b> Potential at Particle \f$y\f$ </b>
 *
 *  We are not in position to determine the potential \f$\Phi(y)\f$ of particle \f$y\f$
 *  \f{eqnarray*}{
 *    \Phi(y) &=& \sum_{i=1}^N u_i \Phi (y,x_i) \\
 *            &=& \sum_{x_i \in E_2(n,l) \bigcup E_3(n,l)} u_i \Phi(y,x_i) \\
 *            &=& \Phi_2^{(n,l)} (y) + \Phi_3^{(n,l)}(y) \\
 *  \f}
 *
 *  where for each \f$x_i\f$ in \f$E_2(n,l)\f$ we calculate \f$\Phi (y,x_i) = \ln(\|y-x_i\|)\f$
 *  directly (not using a series approximation).
 */
#include <stdlib.h>
#include <vector>
#include <complex>
#include <iostream>
#include <string>
#include <limits>
#include <cmath>


#include "Point.h"
#include "Potential.h"
#include "FmmTree.h"
#include "Example1.h"

using namespace std;

int main()
{
  int p = 5;            // p may be upper index of summation in series approximation
  int MAX_NUM_LEVEL = 8;                       // code's maximum refinement level possible
  int DEFAULT_NUM_LEVEL = 3;                   // default refinement level
  int maxClusterThreshold = 5;                 // threshold for particles per cell

  // number of refinement levels is 4 (1, 2, 3, 4)
  // if count starts on 0, then number of refinement
  // levels is 3 (0, 1, 2, 3)
  Example1 example1(4);
  std::vector<Point>  x = example1.getX();
  std::vector<Point>  y = example1.getY();
  std::vector<double> u = example1.getU();


  Potential potential(p);

  // This is the 'build' member function of FmmTree.
  // The call was recursive so took it out of class FmmTree.
  // Build determines the refinement level for cells to stay
  // below the maximum cluster
  // threshold (maximum number of particles per cell)
  int lowest_level_L = MAX_NUM_LEVEL;
  for (int i=DEFAULT_NUM_LEVEL; i<MAX_NUM_LEVEL; ++i)
  {
    FmmTree trial_tree(i, x, y, potential);
    int trial_tree_num_of_levels = trial_tree.getNumOfLevels();
    int trial_tree_cluster_threshold = trial_tree.getClusterThreshold();
    if (trial_tree_cluster_threshold <= maxClusterThreshold)
    {
      lowest_level_L = trial_tree_num_of_levels;
      i = MAX_NUM_LEVEL; // kick out of for loop
    }
    else
    { // do nothing - lowest_level_L will be MAX_NUM_LEVEL
    }
  }

  // for our constructed problem above, the refinement level that
  // allows 4 points per box (cell) is L = 4.  This value is for
  // when counting the levels with level l = 1 as the first level.
  // However, in C++ the index for counts starts with zero.
  // So, you will see the counts start with l = 0.  This means
  // that the refinement level that allows 4 points per box is
  // L = 4 - 1 = 3 when counting the levels of refinement starts
  // with l = 0.
  std::cout << "lowest_level_L = " << lowest_level_L << "\n";

  FmmTree fmmtree(lowest_level_L, x, y, potential);

  std::vector<double> indirect = fmmtree.solve(u);
  std::vector<double> direct = fmmtree.solveDirect(u);


  for (unsigned int i=0; i<direct.size(); ++i)
  {
    std::cout << "direct[" << i << "] = " << direct[i] << " versus "
    		  << "indirect[" << i << "] = " << indirect[i] << "\n";
  }


  double error = 0.0;
  for (unsigned int i = 0; i<direct.size(); ++i)
  {
	  double tmp = std::abs(direct[i]-indirect[i]);
	  if (tmp>error)
        error = tmp;
  }

  std::cout << "Error = " << error << "\n";

  std::cout << "Finished" << "\n";


};




// check to make sure fmmtree works, coords correct
//  std::cout << fmmtree.getNumOfLevels() << '\n';

//  fmmtree.printX();
//  fmmtree.printTreeStructure();

//  int level_for_cluster_threshold = fmmtree.build(threshold);
//  std::cout << "level for cluster threshold = " << threshold
//		    << " is " << level_for_cluster_threshold << '\n';



/**
 * Source and target points explanation
 *
 * In the code below, coordinates for source (x vector) and target (y vector) particles are
 * set at the quarter and 3/4 lengths of each cell in the x and y direction
 * There are four particles per cell with coordinates
 * -  ll_corner_point = P1(x,y) = P1(ll_corner_cell + quarter_length, ll_corner + quarter_length)
 *    lr_corner_point = P2(x,y) = P2(ll_corner_cell + three_quarter_length, ll_corner + quarter_length)
 *    ul_corner_point = P3(x,y) = P3(ll_corner_cell + quarter_length, ll_corner + three_quarter_length)
 *    ur_corner_point = P4(x,y) = P4(ll_corner_cell + three_quarter_length, ll_corner + three_quarter_length)
 * where ll, lr, ul, and ur means lower-left, lower-right, upper-left, and upper-right respectively
 *
 * If the lowest refinement level is L = 3
 * Then the number of cells at L = 3 is 4^L = 4^3 = 64
 * And the number of cells along the length of the unit square domain is 2^L = 2^3 = 8 (8 * 8 = 64)
 * The length of the side of a cell at refinement level l = 3 is therefore (1.0-0.0)/8 = 0.125
 * And if each cell has cell_length = 0.125, then dividing the side length (1.0-0.0) of the unit square domain
 * by the cell length results in (1.0-0.0)/0.125 = 8

 * We work through each cell, incrementing from lower-left cell corner to lower-left cell corner
 * starting at the lower left corner (0.0,0.0) of the unit square domain.
 * The inner for loop increments along the x-axis from lower-left cell corner x-coordinate to lower-left cell
 * corner x-coordinate until reaching the last lower-left cell corner of the 8th cell along the bottom of the
 * unit square domain.  The outer for loop then increments up along the y-axis to the y-coordinate of the lower-
 * left corner for the next set of 8 cells that the inner loop works across.
 *
 * Four target and four source points are collect at each stop of the inner loop.  The points are as described
 * above: ll_corner_point, lr_corner_point, ul_corner_point, and ur_corner_point.  For our example l = 3, there
 * are 4 points per cell * 64 cells = 256 source and target points.  Therefore, the each loop does 8 iterations
 * x_coord and y_coord increment from 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, to 0.875
 *
 */


//  int N = 100;        // set size of vector of points (particles) for x and y
//  int threshold = 5;    // looks like threshold for points in a cell
//  int p = 5;            // p may be upper index of summation in series approximation

//  int level = 3;                               // lowest level of refinement L
//  int nParticlesPerCell = 4;                   // particles per cell at the lowest refinement level L
//  int nTotalParticles                          // total particles in unit square domain
//          = nParticlesPerCell * pow(4,level);
//  int MAX_NUM_LEVEL = 8;                       // code's maximum refinement level possible
//  int DEFAULT_NUM_LEVEL = 3;                   // default refinement level
//  int maxClusterThreshold = 5;                 // threshold for particles per cell



//  std::vector<Point> x;                        // x points - source particles that operate on y particles
//  x.resize((nTotalParticles), Point());
//  std::vector<Point> y(nTotalParticles);       // y points - target particles that are operated on by x particles
//  std::vector<double> u(nTotalParticles);      // u - x points (particles) charge potentials


//  double cells_per_side = pow(2.0,level);         // number of cells along each side of unit square domain
//  double cell_length = (1.0-0.0)/cells_per_side;  // dividing up length of unit square domain into intervals
                                                  // having length of cells of lowest refinement level L
//  double quarter_length = 0.25*cell_length;
//  double three_quarter_length = 0.75*cell_length;

/*
  unsigned int n_points = 0;
  for (double y_coord=0.0; y_coord<1.0; y_coord+=cell_length)
    for (double  x_coord=0.0; x_coord<1.0; x_coord+=cell_length)
    {
      x[n_points].setX(x_coord + quarter_length);       // source particle (point)
      x[n_points].setY(y_coord + quarter_length);       // lower left corner of cell
      y[n_points].setX(x_coord + quarter_length);       // target particle (point)
      y[n_points].setY(y_coord + quarter_length);       // lower left corner of cell

      ++n_points;
      x[n_points].setX(x_coord + three_quarter_length); // source particle (point)
      x[n_points].setY(y_coord + quarter_length);       // lower right corner of cell
      y[n_points].setX(x_coord + three_quarter_length); // target particle (point)
      y[n_points].setY(y_coord + quarter_length);       // lower right corner of cell

      ++n_points;
      x[n_points].setX(x_coord + quarter_length);       // source particle (point)
      x[n_points].setY(y_coord + three_quarter_length); // upper left corner of cell
      y[n_points].setX(x_coord + quarter_length);       // target particle (point)
      y[n_points].setY(y_coord + three_quarter_length); // upper left corner cell

      ++n_points;
      x[n_points].setX(x_coord + three_quarter_length); // source particle (point)
      x[n_points].setY(y_coord + three_quarter_length); // upper right corner of cell
      y[n_points].setX(x_coord + three_quarter_length); // target particle (point)
      y[n_points].setY(y_coord + three_quarter_length); // upper right corner of cell

      ++n_points;
    }

  for (unsigned int i=0; i<u.size(); ++i)
    u[i] = 1.0;

  bool flag_different_x = false;
  bool flag_different_y = false;
  bool flag_different_u = false;

  std::cout << "flag_different_x = " << flag_different_x << "\n";
  std::cout << "flag_different_y = " << flag_different_y << "\n";
  std::cout << "flag_different_u = " << flag_different_u << "\n";

  for (unsigned int i=0; i<x.size(); ++i)
  {
      double maxX = std::max(std::abs(x[i].getCoord()),
    	                         std::abs(x2[i].getCoord()));
      double maxXOne = std::max(1.0,maxX);
      if (std::abs(x[i].getCoord()-x2[i].getCoord())
              <= std::numeric_limits<double>::epsilon()*maxXOne)
      {
        // do nothing - values are same
      }
      else
      {
        flag_different_x = true;
      }

      double maxY = std::max(std::abs(y[i].getCoord()),
    	                         std::abs(y2[i].getCoord()));
      double maxYOne = std::max(1.0,maxY);
      if (std::abs(y[i].getCoord()-y2[i].getCoord())
              <= std::numeric_limits<double>::epsilon()*maxYOne)
      {
        // do nothing - values are same
      }
      else
      {
        flag_different_y = true;
      }

      double maxU = std::max(std::abs(u[i]),
    	                         std::abs(u2[i]));
      double maxUOne = std::max(1.0,maxU);
      if (std::abs(u[i]-u2[i])
              <= std::numeric_limits<double>::epsilon()*maxUOne)
      {
        // do nothing - values are same
      }
      else
      {
        flag_different_u = true;
      }

  }

  std::cout << "Are the x-values different? " << flag_different_x << "\n";
  std::cout << "Are the y-values different? " << flag_different_y << "\n";
  std::cout << "Are the u-values different? " << flag_different_u << "\n";
*/


  /** generate random numbers between 0 and 1 */
//  for (unsigned int i=0; i<N; ++i)
//  {
//    x[i].x = ((double) rand() / (RAND_MAX));
//    x[i].y = ((double) rand() / (RAND_MAX));
//    y[i].x = ((double) rand() / (RAND_MAX));
//    y[i].y = ((double) rand() / (RAND_MAX));
//    u[i] = ((double) rand() / (RAND_MAX));
//  }


//    int level = 3;                               // lowest level of refinement L
//    int nParticlesPerCell = 4;                   // particles per cell at the lowest refinement level L
//    int nTotalParticles                          // total particles in unit square domain
//            = nParticlesPerCell * pow(4,level);
  //  int MAX_NUM_LEVEL = 8;                       // code's maximum refinement level possible
  //  int DEFAULT_NUM_LEVEL = 3;                   // default refinement level
  //  int maxClusterThreshold = 5;                 // threshold for particles per cell



//    std::vector<Point> x;                        // x points - source particles that operate on y particles
//    x.resize((nTotalParticles), Point());
//    std::vector<Point> y(nTotalParticles);       // y points - target particles that are operated on by x particles
//    std::vector<double> u(nTotalParticles);      // u - x points (particles) charge potentials

//    for (unsigned int i=0; i<u.size(); ++i)
//      u[i] = 1.0;
/*
    double cells_per_side = pow(2.0,level);         // number of cells along each side of unit square domain
    double cell_length = (1.0-0.0)/cells_per_side;  // dividing up length of unit square domain into intervals
                                                    // having length of cells of lowest refinement level L
    double quarter_length = 0.25*cell_length;
    double three_quarter_length = 0.75*cell_length;


    unsigned int n_points = 0;
    for (double y_coord=0.0; y_coord<1.0; y_coord+=cell_length)
      for (double  x_coord=0.0; x_coord<1.0; x_coord+=cell_length)
      {
        x[n_points].setX(x_coord + quarter_length);       // source particle (point)
        x[n_points].setY(y_coord + quarter_length);       // lower left corner of cell
        y[n_points].setX(x_coord + quarter_length);       // target particle (point)
        y[n_points].setY(y_coord + quarter_length);       // lower left corner of cell

        ++n_points;
        x[n_points].setX(x_coord + three_quarter_length); // source particle (point)
        x[n_points].setY(y_coord + quarter_length);       // lower right corner of cell
        y[n_points].setX(x_coord + three_quarter_length); // target particle (point)
        y[n_points].setY(y_coord + quarter_length);       // lower right corner of cell

        ++n_points;
        x[n_points].setX(x_coord + quarter_length);       // source particle (point)
        x[n_points].setY(y_coord + three_quarter_length); // upper left corner of cell
        y[n_points].setX(x_coord + quarter_length);       // target particle (point)
        y[n_points].setY(y_coord + three_quarter_length); // upper left corner cell

        ++n_points;
        x[n_points].setX(x_coord + three_quarter_length); // source particle (point)
        x[n_points].setY(y_coord + three_quarter_length); // upper right corner of cell
        y[n_points].setX(x_coord + three_quarter_length); // target particle (point)
        y[n_points].setY(y_coord + three_quarter_length); // upper right corner of cell

        ++n_points;
      }
*/













