# KDE-Based Ensemble Divergence Estimators (EnDive)

[![DOI:10.3390/e20080560](https://zenodo.org/badge/DOI/10.3390/e20080560.svg)](https://doi.org/10.3390/e20080560)

Source code for the EnDive estimator of $f$-divergence functional integrals. Based on the paper [here](https://doi.org/10.3390/e20080560).

## $f$-Divergence Functionals
Divergence functionals are integral functionals of two probability distributions. Divergence functionals play a big role in applications in machine learning, information theory, statistics, and signal processing. In particular, divergence functionals [can be related](https://doi.org/10.1109/DSP-SPE.2015.7369520) to the best probability of error rate of a classification problem. 

While the paper referenced above covers general divergence functionals, this code is written for $f$-divergence functionals which include the [Kullback-Leibler (KL)](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence) divergence, the [Renyi-$\alpha$](https://en.wikipedia.org/wiki/R%C3%A9nyi_entropy) divergence integral, the [Hellinger distance](https://en.wikipedia.org/wiki/Hellinger_distance), the [total variation distance](https://en.wikipedia.org/wiki/Total_variation_distance_of_probability_measures), and the Henze-Penrose or [Dp](https://doi.org/10.1109/TSP.2015.2477805) divergence. $f$-divergence functionals have the form of: 

$D_f(q,p)=\int f\left(\frac{q(x)}{p(x)}\right) p(x)dx,$

where $q$ and $p$ are probability densities and $f$ is some function. For $D_f$ to be a true divergence, certain properties need to be [satisfied](https://en.wikipedia.org/wiki/F-divergence).

## Ensemble Divergence Estimation
EnDive computes an ensemble of kernel density estimators ([KDE](https://en.wikipedia.org/wiki/Kernel_density_estimation)) of the densities $p$ and $q$ with different bandwith values for each estimator. Plugging in the KDEs gives an ensemble of divergence estimators. The EnDive estimator then takes a weighted sum of the ensemble, where the weights are chosen to minimize the mean squared error (MSE). This results in an estimator that achieves an MSE of $O(1/N)$ where $N$ is the number of samples from each density.

The bandwidths can be provided by the user. Otherwise, the default is to compute the set of bandwidths based on the $k$-nearest neighbor distances. 

## References
If you find this work useful, please cite:
```
@article{moon2018endive,
  title={Ensemble estimation of information divergence},
  author={Moon, Kevin R and Sricharan, Kumar and Greenewald, Kristjan and Hero, Alfred O},
  journal={Entropy},
  year={2018},
  volume={20},
  number={8},
  pages={560}
}
```
Other related papers that may be of interest:

[[1](https://doi.org/10.1109/ISIT.2016.7541476)] K.R Moon, K. Sricharan, K. Greenewald, A.O. Hero III, "Improving convergence of divergence functional ensemble estimators," _IEEE International Symposium on Information Theory (ISIT)_, pp. 1133-1137, July 2016.

[[2](https://doi.org/10.1109/DSP-SPE.2015.7369520)] K.R. Moon, V. Delouille, and A.O. Hero III, "Meta learning of bounds on the Bayes classifier error," _IEEE Signal Processing and SP Education Workshop_, pp. 13-18, Aug. 2015.
