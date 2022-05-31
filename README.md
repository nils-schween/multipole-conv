# multipole-conv: A Multipole Moment Converter

## Purpose

In electrodynamics (or gravitational theory) the electric field $\mathbf{E}$ (or
the gravitational field $\mathbf{g}$) is described with the help of a scalar
function, called potential and commonly denoted with $\phi$:
$\mathbf{E} = -\nabla \phi(\mathbf{r})$. The potential is the solution to
Poisson's equation and it is given by the following expression:

$$
4\pi \epsilon_0 \phi(\mathbf{r}) = \int \frac{\rho(\mathbf{r}')}{|\mathbf{r} -\mathbf{r}'|} \mathrm{d}r'^3 .
$$

Where $\rho(\mathbf{r})$ is a charge density.

The potential $\phi$ is very often approximated using the multipole expansion.
There are two forms of the multipole expansion:

1. The Cartesian multipole expansion

$$
	4\pi \epsilon_0 \phi(\mathbf{r}) = \sum_{l=0}^{\infty} \frac{1}{l!} \frac{1}{r^{2l + 1}} Q_{i_{1} \cdots i_{l}}^{(l)} x^{i_{1}} \cdots x^{i_{l}}
$$

2. The spherical (harmonic) multipole expansion

$$
	4\pi \epsilon_0 \phi(\mathbf{r}) = \sum_{l=0}^{\infty} \frac{4\pi}{2l + 1} \frac{1}{r^{l + 1}} q_{l}^{m} Y_{l}^{m}(\theta, \varphi)
$$

The symmetric and traceless tensors $Q_{i_{1} \cdots i_{l}}^{(l)}$ are called
_(Cartesian) multipole moments_. The coefficients $q_{m}^{l}$ are called
_spherical multipole moments_.

The purpose of this command-line tool is to convert between these two kind of
multipole moments. For example, if you have computed the Cartesian multipole
moments for your application, but now you need the spherical multipole moments,
then this command-line tool computes formulae in which you can plug in your
Cartesian multipole moments to get the $q_{l}^{m}$ you are looking for.

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
quadrupole moment $Q_{ij}^{(2)}$ with the corresponding spherical multipole
moments of degree two, namely $q_{2}^{m}$ where $m \in \\{-2,-1,0,1,2\\}$.

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

The other options are explained in detail in the Section [Options](#options).

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
On the left-hand side the spherical multipole moments $q_{l}^{m}$ show up and
the right-hand side the (Cartesian) multipole moments. If you had computed the
quadrupole moment of your application, namely $Q_{ij}^{(2)}$, you could plug in
your results into the above formulae and would be able to convert the
(Cartesian) multipole moments to the spherical multipole moments.

Note that the $Q$'s on the right-hand side have three indices (instead of the
two indices, which you would expect the quadrupole moment to have). These three
indices have the following background: The tensors
$Q_{i_{1} \cdots i_{l}}^{(l)}$ are symmetric.
This means that you can interchange two arbitrary
indices and you get the same component. For example,

$$
Q_{3\mathbf{1}21\mathbf{3}}^{(5)} = Q_{3\mathbf{3}21\mathbf{1}}^{(5)}
$$

Hence, it is always possible to sort the indices such that the component becomes
$Q_{1 \cdots 2 \cdots 3 \cdots}^{(l)}$. This implies that only components with
different numbers of indices whose value is 1, 2 or 3 are different and we can
represent all the components with $p$ indices whose value is 1, $q$ indices
whose value is 2 and $r$ indices whose value is 3 with the symbol
$Q_{pqr}^{(l)}$. Note that $p + q + r = l$, since there are $l$ indices in
total. We take a look at some examples from the first line of the above output,
namely

$$
\begin{align}
Q_{020}^{(2)} &= Q_{22}^{(2)} \\
Q_{002}^{(2)} &= Q_{33}^{(2)} \\
Q_{110}^{(2)} &= Q_{12}^{(2)}
\end{align}
$$

Another important point to note is that on the right-hand side of the output
there are no components of the (Cartesian) multipole moments with $p > 1$. This
means that it is not necessary to compute all the components of the (Cartesian)
multipole moments. It is enough to compute the ones with $p = 0$ or $p=1$. The
reason for this is that the tensors $Q_{i_{1} \cdots i_{l}}$ are traceless, i.e.
you can contract any two of their indices and you will get zero. This statement
is equivalent to

$$
\sum_{k=1}^{3} Q_{i_{1} \cdots k \cdots k \cdots i_{l}}^{(l)} = 0
$$

and it results in a set of equations which allows to express all the components
of the tensors $Q^{(l)} with $p > 1$ in terms of the components with $p=0$ or
$p=1$.

This becomes clear when an example is considered: In Jackson's textbook
"Classical Mechanics" you find the following equation for $q_{2}^{2}$
(corresponding to the first line of the above output)

$$
q_{2}^{2} = \frac{1}{12} \sqrt{\frac{15}{2\pi}} \left(Q_{11}^{(2)} - 2i Q_{12}^{(2)} - Q_{22}^{(2)}\right) .
$$

Note that the component $Q_{11}^{(2)}$ is $Q_{200}^{(2)}$ in the $p,q,r$
notation and hence a component where $p > 1$. The fact that tensors $Q^{(l)}$
are traceless implies that

$$
Q_{11}^{(2)} = -Q_{22}^{(2)} -Q_{33}^{(2)}
$$

and plugging this into Jackson's expression for $q_{2}^{2}$, we recover the first
line of the output of the `multipole-conv` tool.

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

