# Lossy DistFlow
This codebase implements the DistFlow Equations in both single and multiphase with a loss parameterization.

## Single Phase Basic Use
The main function is `distflow_lossy` which takes a MATPOWER case and possible options structure as input.
To run the lossless DistFlow simply use:
```
[v, pf, qf] = distflow_lossy(mpc);
```

Paramerization options are available via the `opt` function.
For example,
```
opt = struct('alpha_method', 7, 'alpha', [0.483, 0.499]);
[v, pf, qf] = distflow_lossy(mpc, opt);
```

## Multiphase Basic Use
The main function is `distflow_multi`, which takes structure arrays `bus` and `branch` as arguments as well as a possible `opt` structure.
Description of the data format is available in the functions.
To run DistFlow using the lossless approximation simply use:
```
[rbus, rbranch] = distflow_multi(bus, branch);
```

Parametrization options are available via the `opt` structure.
For example,
```
opt = struct('alpha_method', 11, 'alpha', [0.492, 0.501]);
[rbus, rbranch] = distflow_multi(bus, branch, opt);
```

Some example data structures are available:
 - `IEEE_13.mat`
 - `IEEE_34.mat`
 - `IEEE_37.mat`
 - `IEEE_123.mat`
