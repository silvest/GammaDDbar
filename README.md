
---

# A tool to determine simultaneously the CKM angle $\gamma$ and the charm mixing and CP-violating parameters 

This repository contains code for determining the CKM angle $\gamma$, charm mixing, and CP-violating parameters in the framework of approximate universality.

## Features

We combine charm and beauty observables in a Bayesian framework. We use the product of Gaussian distributions for each set of correlated data, while the priors of the parameters are assumed to follow uniform distributions. The posterior p.d.f.s of the parameters are obtained through a Markov Chain Monte Carlo (MCMC) algorithm implemented in the Bayesian Analysis Toolkit (BAT) software package. This code can be used to obtain:

- **posterior of the CKM Angle $\gamma$:** four different estimates can be obtained by choosing to combine measurements involving different types of $B$ mesons separately.
- **posteriors of charm mixing and CP-violating parameters in both the dispersive-absorptive and familiar formalisms:** $x_{12}$, $y_{12}$, $\phi_2^M$ and $\phi_2^{\Gamma}$ or, equivalently, $x$, $y$, $\phi_2$ and $\vert q/p \vert -1$.
- **posteriors of beauty and charm decay parameters:** ratios of magnitudes of decay amplitudes and strong phases for the most precise modes available to date.
- **extensible:** new inputs and parameters can be added comfortably by modifying the model class.

Results are stored in BAT output files and ROOT files as one- and two-dimensional histograms.

## Usage 

After having installed ROOT and BAT, the code can be executed by running the bash scripts in the home directory "GammaDDbar". 

1. The classes can be compiled by running the code

   ```
    sh compile_classes.sh
    ```

2. Run the main script using
   ```
    sh compile_main.sh Nchains Nevents_pre Nevents Output_name CombType Var_file
    ```
Here:
- **Nchains** is the number of Markov Chains used.
- **Nevents_pre** is the number of events used to thermalize the MCMC algorithm.
- **Nevents** is the number of configurations of the parameters employed to reconstruct the posteriors.
- **Output_name** is the name of the folder containing the results.
- **CombType** is a number identifying the kind of beauty observables you want to include in the combination. Put '0' for only charged $B$ modes, '1' for only neutral $B$, '2' for only neutral $B_s$ modes and '3' for all the observables.
- **Var_file** is the name of the file containing the parameters for which you want to store the 1D and 2D Histograms in the ROOT file. Examples are already stored in the "Variables" folder.

New inputs and parameters can be added to the combination by editing the class ```MixingModel```. 
The data are stored using the classes ```dato``` and  ```CorrelatedGaussianObservables```.

## Dependencies

- ROOT v6.26/04
- BAT v1.0

## Contributing

Contributions are welcome! If you have suggestions, bug reports, or would like to contribute new features, please open an issue or submit a pull request.

---

