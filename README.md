# multipole-conv: A Multipole Moment Converter

## Purpose

In electrodynamics (or gravitational theory) the electric field $\mathbf{E}$ (or
the gravitational field $\mathbf{g}$) is described with the help of a scalar
function, called potential and commonly denoted with $\phi$:
$\mathbf{E} = -\nabla \phi(\mathbf{r})$. The potential is the solution to
Poisson's equation and it is given by the following expression:

$$
4\pi \epsilon\_0 \phi(\mathbf{r}) = \int \frac{\rho(\mathbf{r}')}{|\mathbf{r} -\mathbf{r}'|} \mathrm{d}r'^3 .
$$

Where $\rho(\mathbf{r})$ is a charge density.

The potential $\phi$ is very often approximated using the multipole expansion.
There are two forms of the multipole expansion:

1. The Cartesian multipole expansion

$$
4\pi \epsilon\_0 \phi(\mathbf{r}) = \sum\_{l=0}^{\infty} \frac{1}{l!} \frac{1}{r^{2l + 1}} Q\_{i\_{1} \cdots i\_{l}}^{(l)} x^{i\_{1}} \cdots x^{i\_{l}}
$$

2. The spherical (harmonic) multipole expansion

$$
4\pi \epsilon\_0 \phi(\mathbf{r}) = \sum\_{l=0}^{\infty} \sum\_{m=-l}^{l} \frac{4\pi}{2l + 1} \frac{1}{r^{l + 1}} q\_{l}^{m} Y\_{l}^{m}(\theta, \varphi)
$$

The symmetric and traceless tensors $Q\_{i\_{1} \cdots i\_{l}}^{(l)}$ are called
*(Cartesian) multipole moments*. The coefficients $q\_{m}^{l}$ are called
*spherical multipole moments*.

The purpose of this command-line tool is to convert between these two kind of
multipole moments. For example, if you have computed the Cartesian multipole
moments for your application, but now you need the spherical multipole moments,
then this command-line tool computes formulae in which you can plug in your
Cartesian multipole moments to get the $q\_{l}^{m}$ you are looking for.

## Compilation

The `multipole-conv` tool needs an additional library to be compiled, namely
the `boost` library. This library is contained in the packages repositories of
must Linux distributions. In case, you have not installed it yet, you can
install it with the following commands:

Ubuntu/Debian
```shell
sudo apt-get install libboost-all-dev
```

Arch Linux/Manjaro
```shell
sudo pacman -Syu boost
```

In a next step clone the repository and call CMake and make, i.e.
```shell
git clone https://github.com/nils-schween/multipole-conv.git
cd multipole-conv
cmake -S . -B build
cd build
make
```

Now there should be a `build` directory, which contains the `multipole-conv`
binary.

## Usage

`multipole-conv` is a command-line tool. The call `./multipole-conv -h` provides
you with an overview of all possible options, i.e.
```
Options:
  -h [ --help ]              Help screen
  -v [ --version ]           Displays the version number.
  -d [ --degree ] arg        Degree of spherical or Cartesian multipole
                             moments.
  -c [ --convention ] arg    Conventions are predefined set of options.
                             Possible values are solid_harmonics, jackson,
                             johnston, real_solid_harmonics.
  --complex arg              Use complex solid harmonics.
  --complex-conjugate arg    Use the complex conjugate of the solid harmonics.
  --normalisation arg        Normalise the (real or complex) solid harmonics.
  --remove-csp arg           Remove the Condon-Shortley phase.
  --include-addition-thm arg Include the addition theorem factor in the
                             definition of the spherical multipole moments.
  --split-addition-thm arg   Include the square root of the addition theorem
                             factor in the definition of the spherical
                             multipole moments. (Requires that the option
                             "inlcude-addition-thm" is set.)
  --cartesian arg            Compute the Cartesian multipole moments.
  --dependent-components arg Compute the dependent components of the Cartesian
                             multipole moments.
  --include-l-factorial arg  Include l! from the Taylor expansion in the
                             definition of the Cartesian multipole moment
```
To get started the most important options are `--degree, -d`, `--convention,c`
and `--cartesian`.

The `-d` options expects an integer as argument and corresponds to the order of
the two multipole expansions. For example, if `multipole-conv` is called with
the option `-d 2`, then it will provide you with the formulae which relate the
quadrupole moment $Q\_{ij}^{(2)}$ with the corresponding spherical multipole
moments of degree two, namely $q\_{2}^{m}$ where $m \in \\{-2,-1,0,1,2\\}$.

The `-c` options expects a string as argument. All possible values are listed in
the help text. This option is meant to ease the handling of the available
options. Every "convention" represents a predefined collection of options. For
example, the convention `jackson` represents a combination of options such that
the definition of the (Cartesian) multipole moments and of the spherical
multipole moments are in agreement with the definitions in Jackson's textbook
"Classical Electrodynamics".

The `--cartesian` option expects either 0 or 1. It is a switch. If it is one,
then the (Cartesian) multipole moments are given as a sum over the spherical
multipole moments. And if it is 0 (which is the default), the `multipole-conv`
tool provides the user with expressions of the spherical multipole moments in
terms of the (Cartesian) multipole moments.

The other options are explained in detail in the section [Options](#options).

An example call of `multipole-conv` would be
``` shell
./multipole-conv -d 2 -c jackson
```

## Interpreting the Output

The output of the above command is

```
Spherical multipole moments

q^+2 _2 = -0.257516-0i Q_0,2,0 -0.128758-0i Q_0,0,2 +0-0.257516i Q_1,1,0
q^+1 _2 = +0+0.257516i Q_0,1,1 -0.257516-0i Q_1,0,1
q^+0 _2 = +0.315392-0i Q_0,0,2
q^-1 _2 = +0+0.257516i Q_0,1,1 +0.257516-0i Q_1,0,1
q^-2 _2 = -0.257516-0i Q_0,2,0 -0.128758-0i Q_0,0,2 +0+0.257516i Q_1,1,0
```
On the left-hand side the spherical multipole moments $q\_{l}^{m}$ show up and
the right-hand side the (Cartesian) multipole moments. If you had computed the
quadrupole moment of your application, namely $Q\_{ij}^{(2)}$, you could plug in
your results into the above formulae and would be able to convert the
(Cartesian) multipole moments to the spherical multipole moments.

Note that the $Q$'s on the right-hand side have three indices (instead of the
two indices, which you would expect the quadrupole moment to have). These three
indices have the following background: The tensors
$Q\_{i\_{1} \cdots i\_{l}}^{(l)}$ are symmetric.
This means that you can interchange two arbitrary
indices and you get the same component. For example,

$$
Q\_{3\mathbf{1}21\mathbf{3}}^{(5)} = Q\_{3\mathbf{3}21\mathbf{1}}^{(5)}
$$

Hence, it is always possible to sort the indices such that the component becomes
$Q\_{1 \cdots 2 \cdots 3 \cdots}^{(l)}$. This implies that only components with
different numbers of indices whose value is 1, 2 or 3 are different and we can
represent all the components with $p$ indices whose value is 1, $q$ indices
whose value is 2 and $r$ indices whose value is 3 with the symbol
$Q\_{pqr}^{(l)}$. Note that $p + q + r = l$, since there are $l$ indices in
total. We take a look at some examples from the first line of the above output,
namely

$$
\begin{align}
Q\_{020}^{(2)} &= Q\_{22}^{(2)} \\
Q\_{002}^{(2)} &= Q\_{33}^{(2)} \\
Q\_{110}^{(2)} &= Q\_{12}^{(2)}
\end{align}
$$

Another important point to note is that on the right-hand side of the output
there are no components of the (Cartesian) multipole moments with $p > 1$. This
means that it is not necessary to compute all the components of the (Cartesian)
multipole moments. It is enough to compute the ones with $p = 0$ or $p=1$. The
reason for this is that the tensors $Q\_{i\_{1} \cdots i\_{l}}$ are traceless,
i.e. you can contract any two of their indices and you will get zero. This
statement is equivalent to

$$
\sum\_{k=1}^{3} Q\_{i_{1} \cdots k \cdots k \cdots i\_{l}}^{(l)} = 0
$$

and it results in a set of equations which allows to express all the components
of the tensors $Q^{(l)} with $p > 1$ in terms of the components with $p=0$ or
$p=1$.

This becomes clear when an example is considered: In Jackson's textbook
"Classical Mechanics" you find the following equation for $q\_{2}^{2}$
(corresponding to the first line of the above output)

$$
q\_{2}^{2} = \frac{1}{12} \sqrt{\frac{15}{2\pi}} \left(Q\_{11}^{(2)} - 2i Q\_{12}^{(2)} - Q\_{22}^{(2)}\right) .
$$

Note that the component $Q\_{11}^{(2)}$ is $Q\_{200}^{(2)}$ in the $p,q,r$
notation and hence a component where $p > 1$. The fact that tensors $Q^{(l)}$
are traceless implies that

$$
Q\_{11}^{(2)} = -Q\_{22}^{(2)} -Q\_{33}^{(2)}
$$

and plugging this into Jackson's expression for $q\_{2}^{2}$, we recover the
first line of the output of the `multipole-conv` tool.

Up to now, we looked at the output of the `multipole-conv` tool for the
spherical multipole moments. Its output for the (Cartesian) multipole moments is
obtained with the call

``` shell
./multipole-conv -d 2 -c jackson --cartesian 1
```
and it is
```
Cartesian multipole moments (independent components = multipole basis functions)

Q_0,2,0 = -1.94163-0i q^+2_2 -1.58533-0i q^+0_2 -1.94163-0i q^-2_2
Q_0,1,1 = +0-1.94163i q^+1_2 +0-1.94163i q^-1_2
Q_0,0,2 = +3.17066-0i q^+0_2
Q_1,0,1 = -1.94163-0i q^+1_2 +1.94163-0i q^-1_2
Q_1,1,0 = +0+1.94163i q^+2_2 +0-1.94163i q^-2_2

Cartesian multipole moments (dependent components)

Q_2,0,0 = +1.94163+0i q^+2_2 -1.58533+0i q^+0_2 +1.94163+0i q^-2_2
```
On the left-hand side we see the components of the quadrupole moment expressed
in terms of the spherical multipole moments of degree two. For the components of
the quadrupole moment the $p,q,r$ notation is used. The components with
$p \leq 1$ are called independent components and the ones with $p > 1$ are
called dependent components.

**Note:** Currently, the `multipole-conv` command-line tool works without
problems up to degree 15. Afterwards numerical problems appear and its results
must be looked upon with care. In case, you need more precision, just file an
issue and we will try to make it work.


## Options

The `multipole-conv` command-line tool provides options to deal with the
different forms of the multipole expansions used in different areas of physics.
The available options can be split into two parts: options changing the
spherical multipole expansion and options changing the (Cartesian) multipole
expansion.

### Spherical multipole expansion

The options for the spherical multipole expansion concern the spherical harmonic
functions, which are part of the definition of the spherical multipole moments.

```
--complex and --complex-conjugate
```
If **complex spherical harmonics** are used the spherical multipole expansion is

$$
4\pi \epsilon\_0 \phi(\mathbf{r}) = \sum\_{l=0}^{\infty}  \sum\_{m=-l}^{l} \frac{4\pi}{2l + 1} \frac{1}{r^{l + 1}} q\_{l}^{m} Y\_{l}^{m}(\theta, \varphi)
$$

and the spherical multipole moments are defined as

$$
q\_{l}^{m} \equiv \int \rho(\mathbf{r}') r'^{l} Y\_{l}^{m*}(\theta', \varphi') \mathrm{d}r'^3
$$

Instead of using complex spherical harmonics $Y_{l}^{m}$ it is possible to use
**real spherical harmonics** $Y_{lms}$. Then the expansion looks like

$$
4\pi \epsilon\_0 \phi(\mathbf{r}) = \sum\_{l=0}^{\infty} \sum\_{m=0}^{l} \sum\_{s=0}^{l} \frac{4\pi}{2l + 1} \frac{1}{r^{l + 1}} q\_{lms} Y\_{lms}(\theta, \varphi)
$$

and the definition of the real spherical multipole moments is

$$
q\_{lms} \equiv \int \rho(\mathbf{r}') r'^{l} Y\_{lms}(\theta', \varphi') \mathrm{d}r'^3
$$

If complex spherical harmonics are used in the multipole expansion, then the
spherical multipole moment's definition contains complex-conjugate spherical
harmonics and it is necessary to add the option `--complex 1` and
`--complex-conjugate 1` . If this options are not set, then the real spherical
harmonics are used, because they are contained in definition of the real
spherical multipole moments $q\_{lms}$.

```
--normalisation
```

Real and complex spherical harmonics can be normalised or not. It they are
normalised, as in the case as in the definitions of $q\_{m}^{l}$ and $q\_{lms}$,
then the option `--normalisation 1` must be included in the call of the
`multipole-conv` tool.

We denote the **complex spherical harmonics without normalisation** with
$\tilde{Y}\_{l}^{m}$. If these are used in the multipole expansion, then we get

$$
4\pi \epsilon\_0 \phi(\mathbf{r}) = \sum\_{l=0}^{\infty} \sum\_{m=-l}^{l} \frac{(l-m)!}{(l+m)!} \frac{1}{r^{l + 1}} \tilde{q}\_{l}^{m} \tilde{Y}\_{l}^{m}(\theta, \varphi)
$$

and, accordingly,

$$
\tilde{q}\_{l}^{m} \equiv \int \rho(\mathbf{r}') r'^{l} \tilde{Y}\_{l}^{m*}(\theta', \varphi') \mathrm{d}r'^3
$$

The **real spherical harmonics without normalisation** are denoted with
$\tilde{Y}\_{lms}$ and are part of the following multipole expansion

$$
4\pi \epsilon\_0 \phi(\mathbf{r}) = \sum\_{l=0}^{\infty} \sum\_{m=0}^{l} \sum\_{s=0}^{l} \frac{2}{1 + \delta_0m} \frac{(l-m)!}{(l+m)!} \frac{1}{r^{l + 1}} \tilde{q}\_{lms} \tilde{Y}\_{lms}(\theta, \varphi)
$$

$$
\tilde{q}\_{lms} \equiv \int \rho(\mathbf{r}') r'^{l} \tilde{Y}\_{lms}(\theta', \varphi') \mathrm{d}r'^3
$$

The default value of `--normalisation` is `0`. Moreover, when `./multipole-conv`
is called with only `--degree` (or, alternatively `-d`), then the above
multipole expansion with its definition of $\tilde{q}\_{lms}$ is used.

```
--include-addition-theorem
```

In all four versions of the spherical multipole expansion there is an additional
numerical factor (e.g. $4\pi/(2l + 1)) right after the summation symbols. This
factor is a consequence of the [addition theorem for spherical
harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics#Addition_theorem).
Some authors include it in their definition of the spherical multipole moments.
If this is the case, then `multipole-conv` should be called with
`--include-addition-theorem 1`. Its default value is `0`.

``` 
--split-addition-theorem
```
This *addition theorem factor* is sometimes split, namely

$$
\frac{4\pi}{2l + 1} = \sqrt{\frac{4\pi}{2l + 1}} \sqrt{\frac{4\pi}{2l + 1}}
$$

and its square root is then included in the definition of the spherical
multipole moment. If this is the case, `multipole-conv` should be called with
`--split-addition-theorem 1`. Its default value is `0`.

```
--remove_condon_shortley_phase
```
 
 Very often the definition of the associated Legendre Polynomials contains a
 factor $(-1)^{m}$, which is called the Condon-Shortley phase. Sometimes this
 factor is not included in the definition of the Legendre Polynomials, but in
 the definition of the (real) spherical harmonics.
 
 If this factor is not included in the associated Legendre Polynomials **nor**
 in the spherical harmonics, it is not included in the definition of the
 spherical multipole moments and needs to be removed. In this case
 `multipole-conv` should be called with `--remove_condon_shortley_phase 1`.
 
### Cartesian multipole moments

```
--cartesian
```

If set to `1`, `multipole-conv` computes the (Cartesian) multipole moments in
terms of the spherical multipole moments as described at the end of the
[Usage](#usage) section. Its default value is `0`.

```
--dependent_components
```

As explained in the section [Interpreting the Output](#interpreting-the-output),
we distinguish between independent components and dependent components of the
(Cartesian) multipole moments. If `multipole-conv` is called with
`--dependent_componets 1`, then it also outputs formulae for the dependent
components. Its default value is `0`.

```
--include_l_factorial
```

The (Cartesian) multipole expansion is

$$
4\pi \epsilon\_0 \phi(\mathbf{r}) = \sum\_{l=0}^{\infty} \frac{1}{l!} \frac{1}{r^{2l + 1}} Q\_{i\_{1} \cdots i\_{l}}^{(l)} x^{i\_{1}} \cdots x^{i\_{l}}
$$

The factor $1/l!$ can be included in the definition of the (Cartesian) multipole
moment $Q^{(l)}$. If this is the case, `multipole-conv` should be called with
`--include_l_factorial 1`. Its default value is `0`.

## Further applications

### Plasma physics

The dynamics of a plasma are modelled with the Boltzmann equation. In some cases
it can be useful to expand the distribution function $f$ in terms of spherical
harmonics or in terms of Cartesian "tensors", namely

$$
\begin{align}
	f(\mathbf{r}, \mathbf{p}, t) &= \sum\_{l=0}^{\infty} \sum\_{m=0}^{l} \sum\_{s=0}^{1} (-1)^m \tilde{f}\_{lms}(\mathbf{r},p,t) \tilde{Y}\_{lms}(\theta,\varphi) \\
	f(\mathbf{r}, \mathbf{p}, t) &= \sum\_{l=0}^{\infty} F\_{i\_{1} \cdots i\_{l}}^{(l)}(\mathbf{r},p,t) \frac{p^{i\_{1}} \cdots p^{i\_{l}}}{p^{l}}
\end{align}
$$

as done by Johnston in his 1960 paper "Cartesian Tensor Scalar Product and
Spherical Harmonic Expansions in Boltzmann's Equation". Note that Johnston
defines $f\_{lms} \equiv (-1)^m \tilde{f}\_{lms}$. To express the components of
the Cartesian tensors $F^{(l)}$ in terms of the $f\_{lms}$ you can use the
convention `johnston`. For example, the call `./multipole-conv -d 2 -c johnston`
produces

``` 
Cartesian multipole moments (independent components = multipole basis functions)

Q_0,2,0 = -3 q_2,2,0 -0.5 q_2,0,0
Q_0,1,1 = 1.5 q_2,1,1
Q_0,0,2 = 1 q_2,0,0
Q_1,0,1 = 1.5 q_2,1,0
Q_1,1,0 = 3 q_2,2,1

Cartesian multipole moments (dependent components)

Q_2,0,0 = 3 q_2,2,0 -0.5 q_2,0,0
```

It is, of course, possible to work with complex and normalised spherical
harmonics instead of using $\tilde{Y}\_{lms}$. If you have troubles finding the
right set options, contact us.

### Basis transformation in the space of homogeneous and harmonic polynomials

The space of homogeneous and harmonic polynomials of degree $(l)$ is denoted
with $\mathcal{H}^{l}(\mathbb{R}^{3})$. A basis of this space are the solid
harmonics of degree $l$, i.e.

$$
\mathcal{H}^l(\mathbb{R}^{3}) = \text{span}\\{r^l Y\_{l}^{m} \mid -l \leq m \leq l\\}
$$

Another basis of this space are the multipole basis functions $M\_{pqr}^(l)$ as
defined in our paper (see section [Cite our paper](#cite-our-paper)), i.e.

$$
\mathcal{H}^l(\mathbb{R}^{3}) = \text{span}\\{ M\_{pqr}^(l) \mid p \leq 1 \text{ and } p+q+r=l\\}
$$

The `multipole-conv` tool can compute the basis transformation between these two
bases. If you call it with the convention `solid_harmonics` ( or
`real_solid_harmonics`) it outputs expressions for the solid harmonics (real
solid harmonics) in terms of the multipole basis functions.

## 
