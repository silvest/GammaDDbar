
---

# A tool to determine simultaneously the CKM angle $\gamma$ and the charm mixing and CP-violating parameters 

This repository contains code for determining the CKM angle $\gamma$, charm mixing and CP-violating parameters in the framework of approximate universality.

## Features

We combine charm and beauty observables in a Bayesian framework. We use the product of Gaussian distributions for each set of correlated data, while the priors of the parameters are assumed to follow uniform distributions. The posterior p.d.f.s of the parameters are obtained through a Markov Chain Monte Carlo (MCMC) algorithm implemented in the Bayesian Analysis Toolkit (BAT) software package. This code can be used to obtain:

- **posterior of the CKM Angle $\gamma$:** four different estimates can be obtained by choosing to combine measurements involving different initial $B$ mesons separately.
- **posteriors of charm mixing and CP-violating parameters in both the dispersive-absorptive and familiar formalisms:** $x_{12}$, $y_{12}$, $\phi_2^M$ and $\phi_2^{\Gamma}$ or, equivalently, $x$, $y$, $\phi_2$ and $\vert q/p \vert -1$.
- **posteriors of beauty and charm decay parameters:** ratios of magnitudes of decay amplitudes and strong phases for the most precise modes available to date.
- **extensible:** new inputs and parameters can be added comfortably by modifying the model class.

Results are stored in BAT output files and in ROOT files as one- and two-dimensional contours for the variables contained in the folder "GammaDDbar/Variables".

## Usage 

After having installed ROOT and BAT, the code can be executed by running the bash scripts in the home directory "GammaDDbar".

1. The classes can be compiled by running the code

   ```
    cd GammaDDbar
    sh compile_classes.sh
    ```

3. Run the main script using "compile_main.sh" by sepcifying the Number of Markov chains, Numer of events to thermalize the chains, Number of events to sample to reconstruct the posterior of the parameters, Name of the output folder, combination type and name of the .txt file containing the variables to populate the histograms in the ROOT file.

## Dependencies

- ROOT v 6.26/04
- BAT v 1.0

## Contributing

Contributions are welcome! If you have suggestions, bug reports, or would like to contribute new features, please open an issue or submit a pull request.

---

