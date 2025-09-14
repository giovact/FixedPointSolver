# FixedPointSolver

Iterative fixed-point solver for integral self-consistency systems of equations


#### General description
It's a basic library to solve systems of self-consistent equations, the usual problem in computing equilibrium phase diagrams for Mean-Field theories in the thermodynamic limit. In statistical physics, the equilibrium behavior of a system in the thermodynamic limit is determined by an extremization condition (through a saddle-point integration) of a free energy function depending on a finite set of order parameters. Thermodynamics is determined by (local) minima of this object. Typical situation in spin-glass theories there are integrals involved, reminiscent of the quenched noise. 

In the end it's just a library to try to help make things easier and faster for a small community of people. 

- Is it easy to customize? YES

### Features
1. modular structure: a new model and its self-consistent equations can be added as a single julia file in the folder `src/models/` and all the built-in functions can be used
2. parallelized implementation (multi-thread) of custom scripts to 
   - explore a 2-dimensional space in the control parameters with different initial conditions 
   - parallel scan solution by following previous fixed points
   - computation of critical lines in 2-d planes 

## Installation

The package is not yet registered. To use it, clone the repository locally:
> git clone https://github.com/giovact/FixedPointSolver


## Basic Usage
The bulk abstract type of this module is `FPModel`. 

Any model must be defined as a struct `MyModel<:FPModel`, with the following 3 required fields:
```
struct MyModel<:FPModel
    params  :: NamedTuple{ (:param_1,:param_2,...,:param_n),
                           Tuple{type_1,type_2,...,type_n}}
    O       :: NamedVec                                                                 
    aux     :: Dict{String,Float64}                        
```

NamedVec is an alias of NamedArray with default type Float64 (see `src/types.jl`).

* `params <: NamedTuple`: NamedTuple (immutable) containing name/value of each control parameter of the model. All names are parsed as Symbols, value types can be any scalar, but must be consistent with constructors. 
* `O::NamedVec`: NamedArray containing name/value of the order paramers, to be updated until convergence.  
* `aux::Dict{String,Float64}`: auxiliary dictionary (so with mutable scalar entries) used to store additional optional observables. The free energy computed at convergence will be stored in this dictionary.

Additional fields can be declared if needed. For instance, the typical case is to use a conjugate set of order parameters, this can be added as an additional field `Oconj::NamedVec` in the structure. However, currently convergence criterion is defined only on the field `O`.

Suppose now that the `src/models` folder you created a file `MyModel.jl` with the above structure declaration. There are some other ingredients the file must contain. A first rule of thumb, all functions written in `MyModel.jl` will be dispatched on the type `MyModel` to avoid inconsistencies. 

A detailed guide on how to define your custom model is given in `src/models/examples/model_template.jl`

The file `src/models/HopfieldExample.jl` contains an example file containing the self-consistent equations defining the MF theory for the Hopfield model of associative memory.


## Parallelized scripts
In the folder `scripts/` there are 3 files that perform different parallel operations on multi-threads

Remember to launch any one of the above scripts with option `julia -t $nthreads ` or `julia -t auto ` when using either one of these for faster computation.  

## TODO
1. Fix logarithms for numerical stability (e.g. `log(cosh(x))` to be changed in `log2cosh(x)` or `log(Herf(x))` to be changed in `logHerf(x)`)