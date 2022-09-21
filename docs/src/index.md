```@meta
CurrentModule = BatchReactor
```

# BatchReactor
Batch reactor is a package for simulating plug flow reactor with detailed surface  or gasphase chemistry. Additionally you may use user defined function for the calculation of reaction rates. 

Documentation for [BatchReactor](https://github.com/vinodjanardhanan/BatchReactor.jl).

## Installation
To install the package, use the following commands in the julia REPL
```julia
julia> using Pkg
julia> Pkg.add("BatchReactor")
```

## General interfaces
```@index
```

```@autodocs
Modules = [BatchReactor]
```

## Governing equations
The batch model simulates a constant volume batch reactor.  The code integrates the following governing equation 
```math
\frac{d(\rho_k)}{dt} = \dot{s_k} M_k A_{sv} + \dot{\omega_k}M_k
```

```math
\rho = \sum \rho_k
```

```math
p = \frac{\rho RT}{\bar{M}}
```
Here $\rho_k$ is the mass density of species $k$ defined as $\rho Y_k$, where $Y_k$ is the mass fraction of species $k$ and $\rho$ is the density (kg/m$^3$) of the mixture. $M_k$ is the molecular weight of species $k$ (Kg/mol), $\dot{s}_k$ is the molar production rate (mol/m$^2$-s) of species $k$ due to surface reactions, $\dot{\omega}_k$ is the molar production rate (mol/m$^3$-s) of species $k$ due to gas phase reactions, and $A_{sv}$ is the surface area to volume (1/m). 

# Executing the code
## Surface chemistry
For solving a surface chemistry problem: On the Julia REPL 
```julia
julia>using Plug
julia>plug("batch.xml","lib/", surfchem=true)
```
## Gasphase chemistry
For solving a gasphase chemistry problem: On the Julia REPL 
```julia
julia>using Plug
julia>plug("batch.xml", "lib/", gaschem=true)
```

In the above calls, it is assumed that the input file *plug.xml* is present in the working directory and *../lib/* is the path to the *lib* directory relative to the current working directory. The function can be called from any directory and in that case the first argument must point to the *plug.xml* file relative to the current working directory. The output files will be generated in the directory where *plug.xml* is present. In general the function takes three optional keyword arguments *gaschem*, *surfchem*, and *sens*. *gaschem* must be true to simulate gasphase chemistry, *surfchem* must be true for surface chemistry, and *sens* must be true whenever sensitivity analysis is performed. 

## User defined chemistry
For solving the model with user defined chemistry: On the Julia REPL 
```julia
julia>using Plug, ReactionCommons
julia>plug("batch.xml", "lib/", udf)
```
*udf* is a function having the following signature
```
function udf(state::UserDefinedState)
```
where state is a structure defined as follows
```
struct UserDefinedState
    T::Float64
    p::Float64
    molefracs::Array{Float64,1}
    molwt::Array{Float64,1}
    species::Array{String,1}
    source::Array{Float64,1}
end
```
The program expects the species source terms in *source* mols/m3-s depending on whether you are solving surface chemistry problem or gasphase chemistry problem. The example call provided in the *runtests.jl* returns zero molar production and consumption rates. Within the code the source terms are multiplied with the molecular weight. The order ing of species source terms must be same as the order in wich the species appears in UserState.species.


## Input file for surface chemistry
The method takes *file\_name* as the argument. The file_name points to the input XML file. The structure of the XML input file is shown below.

```
?xml version="1.0" encoding="ISO-8859-1"?>
<batch>
    <gasphase>CH4 H2O H2 CO CO2 O2 N2</gasphase>
    <molefractions>CH4=0.25,H2O=0.25,N2=0.5</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <Asv>10</Asv>
    <time>10</time>
    <surface_mech>ch4ni.xml</surface_mech>
</batch>
```
## Input file for gasphase chemistry
The method takes *file\_name* as the argument. The file_name points to the input XML file. The structure of the XML input file is shown below.

```
?xml version="1.0" encoding="ISO-8859-1"?>
<batch>    
    <molefractions>CH4=0.25,H2O=0.25,N2=0.5</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <Asv>10</Asv>
    <time>10</time>
    <gas_mech>h2o2.dat</gas_mech>
</batch>
```

The major difference between the input file for surface chemistry problem and gasphase chemistry problem is the <gasphase> tag of xml input. In the case of gasphase chemistry problem, the participating species are read from the mechanism input file, which is specified using the <gas_mech> tag

The meaning of different tags is specified below.

- <batch> : The root XML tag for batch
- <gasphase> : list of gas-phase species. The species names must be separated by white spaces or tab
- <molefractions> : inlet mole fraction of the gas-phase species. Instead of mole fractions, mass fractions may also be specified. In that case, the tag must be <massfractions>. You must ensure that the sum of mass or mole fractions specified is unity. There are no internal checks to ascertain this. 
- <T>: operating temperature in K
- <p>: initial pressure in Pa
- <Asv> : Surface area to volume ratio (1/m)
- <time> : integration time in s
- <surface_mech> : name of the surface reaction mechanism
- <gas_mech> : name of the gasphase reaction mechanism

## Output
The code generates two output files in the same directory where the input file **`batch.xml`** is present. 
The file **`gas_profile.dat`** contains the mole fraction of the gas phase species as a function of time.
The file **`surf_covg.dat`** contains the surface coverages of adsorbed species as a function of time. 
In addition to these two files, the code also generates terminal output, which shows integration progress.
The terminal output is nothing by the integration time. 

An example terminal output for surface chemistry problem is shown below.

```
julia> batch("batch_surf/batch.xml","lib/", surfchem=true)
0.000000e+00
2.116852e-16
4.233704e-16
9.294022e-16
1.435434e-15
2.617619e-15
..
..
..
9.925172e+00
9.953705e+00
9.973533e+00
9.993361e+00
1.000000e+01
:Success

```

A sample output of **`gas_profile.dat`** is shown below. 
```
         t           T           p         rho         CH4         H2O          H2          CO         CO2          O2          N2
0.0000e+00  1.0732e+03  1.0000e+05  2.5240e-01  2.5000e-01  2.5000e-01  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  5.0000e-01
5.3800e-12  1.0732e+03  1.0000e+05  2.5240e-01  2.5000e-01  2.5000e-01  2.8224e-17  2.8064e-31  3.6828e-43  3.5083e-34  5.0000e-01
7.3222e-12  1.0732e+03  1.0000e+05  2.5241e-01  2.5000e-01  2.5000e-01  7.2345e-17  5.7521e-31  5.7355e-42  1.1037e-33  5.0000e-01
1.3893e-11  1.0732e+03  1.0000e+05  2.5241e-01  2.5000e-01  2.5000e-01  5.7775e-16  5.2239e-29  1.5426e-39  2.3067e-32  5.0000e-01
1.7496e-11  1.0732e+03  1.0000e+05  2.5241e-01  2.5000e-01  2.5000e-01  1.2389e-15  2.4341e-28  1.3782e-38  6.7062e-32  5.0000e-01
...
...
9.8894e+00  1.0732e+03  1.0205e+05  2.5241e-01  2.3492e-01  2.2984e-01  3.5247e-02  4.9651e-03  5.0865e-03  3.8647e-17  4.8994e-01
9.9840e+00  1.0732e+03  1.0206e+05  2.5241e-01  2.3483e-01  2.2971e-01  3.5461e-02  4.9881e-03  5.1229e-03  1.5330e-17  4.8988e-01
1.0000e+01  1.0732e+03  1.0207e+05  2.5241e-01  2.3481e-01  2.2969e-01  3.5498e-02  4.9920e-03  5.1290e-03  1.8337e-17  4.8987e-01
```

A sample output of **`surf_profile.dat`** is shown below. 
```
         t           T        (NI)       H(NI)       O(NI)     CH4(NI)     H2O(NI)     CO2(NI)      CO(NI)      OH(NI)       C(NI)     HCO(NI)      CH(NI)     CH3(NI)     CH2(NI)
0.0000e+00  1.0732e+03  6.0000e-01  0.0000e+00  0.0000e+00  0.0000e+00  4.0000e-01  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00
5.3800e-12  1.0732e+03  6.0842e-01  2.9832e-04  3.2879e-05  1.1712e-09  3.9101e-01  2.7168e-33  5.9257e-21  2.3256e-04  1.0792e-12  2.3196e-23  8.8535e-13  4.6995e-11  8.9031e-12
7.3222e-12  1.0732e+03  6.1140e-01  4.1988e-04  5.9207e-05  1.1838e-09  3.8782e-01  5.7960e-32  1.8386e-20  3.0147e-04  2.8680e-12  8.2649e-23  1.6840e-12  5.9925e-11  1.4892e-11
1.0920e-11  1.0732e+03  6.1683e-01  6.6111e-04  1.2488e-04  1.1946e-09  3.8197e-01  1.7656e-30  2.0089e-19  4.1135e-04  9.6114e-12  3.6587e-22  3.4522e-12  7.8856e-11  2.6942e-11
...
...
...
9.8894e+00  1.0732e+03  7.7820e-01  1.0109e-01  3.5083e-02  1.4451e-09  5.2625e-04  5.6835e-06  8.4964e-02  1.2547e-04  9.1985e-06  6.5892e-13  3.1341e-11  4.1784e-10  2.2512e-10
9.9840e+00  1.0732e+03  7.7791e-01  1.0137e-01  3.4840e-02  1.4442e-09  5.2584e-04  5.6713e-06  8.5216e-02  1.2498e-04  9.2286e-06  6.6385e-13  3.1254e-11  4.1700e-10  2.2450e-10
1.0000e+01  1.0732e+03  7.7787e-01  1.0141e-01  3.4799e-02  1.4441e-09  5.2576e-04  5.6693e-06  8.5259e-02  1.2490e-04  9.2337e-06  6.6469e-13  3.1239e-11  4.1686e-10  2.2439e-10
```