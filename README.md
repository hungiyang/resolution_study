# resolution_study
## abstract
In this work, we characterize the accuracy with which simulations of a given resolution can model galaxy clustering statistics.  We focus on models in which galaxies are directly matched to resolved halos and subhalos in simulations, using the subhalo abundance matching approach.
We investigate different abundance matching proxies and a range of galaxy number densities.  We present a model which characterizes the missing pairs in the 3D two-point correlation function due to inadequate resolution. We use a suite of simulations to test these effects explicitly for mass resolutions ranging from $\sim 3 \times 10^8$ to $\sim 1 \times 10^{10} \Msunh$, in all cases by using a paired simulation with at least eight times higher resolution. We suggest that a rough rule of thumb for the necessary resolution may be achieved by assuring that the impact of these missing pairs does not dominate over the sample variance of the relevant sample.
We use this approach to calculate the required resolution for several current and future galaxy samples, including those from the SDSS, GAMA, BOSS, and DESI surveys.

## ipython notebooks
`filepath.py` contains the directory in which all the pair counts resoults are saved. `notebook_load_functions.py` contains all the `load_*()` functions that load the pair count results, and also the luminosity-number density relation. We import the functions in `notebook_load_functions.py` at the beginning of each notebook, and then proceed to make the plots we want.

### `paper_plots.ipynb`
This notebook contains the code to produce all the plots in the paper.

### `simulate_low_resolution.ipynb`
Simulate the low resolution mass function with the high resolution halo catalog by randomly selecting samples in each mass bin. The plot compares the pair counts of simulate low resolution with high resolution results, and we don't seem as many missing pairs as in the actual low resolution case.

### `Multidark.ipynb`
Shows sample variance with different box sizes, redshifts, and simulations (high, low resolution Multidark, and Darksky). We also investigated the variation of proxy cuts for the subboxes.

### `projected_correlation_function.ipynb`
Investigate the the sample variance of the projected correlation function.(instead of the 3D correlation function) The peak around 1 Mpc still exists, but is slighly smoothed out.

### `resolution_vs_sample_variance.ipynb`
The notebook I used to figure out how to model missing pairs and sample variance. The final results are applied in `paper_plots.ipynb`

### `requirement_plot.ipynb`
Contains all the plots used in the paper, but also some additional code/plots.


## scripts (in `scipts/`)
`scripts/` contains the code to produce the pair counts in the analysis. Note that they are copied from `/nfs/slac/g/ki/ki22/cosmo/iameric/resolution_study`, so it is necessary to change the save directory in `save_\*.py` files before running them. `save_\*.py` contains the code, and `submit_\*.py` submits the batch jobs.
#### missing pairs: 
pair counts of all halos (Chinchilla boxes)  

* `save_pairs.py`  
in `save_pairs.py`, there is the option `--JK` so that the code uses `count_pairs_JKindex.py` and save the Jackknife index of the pairs. We then use the 
* `submit_all.py`  

#### sample variance: 
pair counts of all halos for subboxes (Multidark, Darksky)
* `save_multidark.py`  
* `submit_all_multidark.py`

#### HOD scheme: 
pair counts for central halos only (Chinchilla)
* `save_pair_counts_HOD.py`
* `submit_all_HOD_pair_counts.py`

#### projected correlation function: 
projected pair counts of all halos for subboxes (Multidark)
* `save_projected_pair_counts_Multidark.py`
* `submit_all_projected_correlation.py`

#### proxy cut for sample variance:
save the proxy cuts for subboxes of different box sizes in the sample variance study
* `save_proxycut.py`
* `submit_all_proxycut.py`







