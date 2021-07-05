# Pulsar Glitch Analyze with Stride Boxcar Fitting

## Recipe for fitting

Turn on tempo2 fit for RAJ, DECJ, F0, F1, F2, JUMP. If exist PMRA, PMDEC, turn on tempo2 fit for DM, DM1, …, PX, PMRA, PMDEC in dbr_psrn.par. (Remove GLPH)

```vim dbr_psrn.par```

Fit dbr_psrn.par chp_psrn.tim with ENTERPRISE

```/local/scratch/yangliu/glitch/run_enterprise/run_enterprise.py --gl-all --auto-add --dynesty --nlive 500 -t 2 dbr_psrn.par chp_psrn.tim --truth-file dbr_psrn_truth.txt --dynesty-plots --plot-derived -j --red-prior-log --measured-prior -A -8 --tspan-mult 1.1 --glitch-alt-f0 --glitch-alt-f0t 200 --alt-f0t-gltd --glitch-epoch-range 100 --measured-sigma 50 --glitch-td-min 1 --glitch-td-max 2.5 --glitch-f0d-range 3.0 --glitch-f0-range 0.8 --glitch-f1-range 0.8 --glitch-td-split 2```

Open fix_dbr_psrn_d.par.post with fix_chp_psrn.tim in tempo2 and save a new par file fnl_psrn_d.par, \_d stands for dynesty fit

```tempo2 -gr plk -f dbr_psrn_d.par.post chp_psrn.tim```

Copy fnl_psrn_d.par to tst_psrn_d.par

```cp fnl_psrn_d.par tst_psrn_d.par```

Trun off the fit for all parameters other than F0, F1, F2 in tst_psrn_d.par. Remove all glitch parameters except GLEP, remove red/white noise parameters

```vim tst_psrn_d.par```

Make combined plots and nudot.asc etc.

```python new_make_pulsar_plots.py fnl_psrn_d.par chp_psrn.tim tst_psrn_d.par```

Generate stride epochs and save it in text files, c stands for cadence, the expression in brackets are reference values for the choice of boxcar width and step size. Do tempo2 fitting for F0, F1, F2 in all the boxcars and save results in psrn_g_w_s_data.txt, g: numbber of glitches, w: boxcar width, s: step size

```python stride_fit.py -p tst_psrn_d.par -t chp_psrn.tim -w (5*c) -s (2.5*c) -g (GLEP_1, …)```

Calculate the results for panels and save them in text files. Plot the evolution of nu and nudot

```python plot_panels.py -p fnl_psrn_d.par -s psrn_g?_w?_s?_data.txt```


## Reference
* Ellis J. A., Vallisneri M., Taylor S. R., Baker P. T., 2019, Astrophysics Source Code Library,p. ascl:1912.015
* Foreman-Mackey D., Hogg D. W., Lang D., Goodman J., 2013, Publications of the Astro-nomical Society of the Pacific, 125, 306
* Hobbs G. B., Edwards R. T., Manchester R. N., 2006, Monthly Notices of the Royal Astro-nomical Society, 369, 655
* Shaw B., et al., 2018b, Monthly Notices of the Royal Astronomical Society, 478, 3832
* Shaw B., Keith M. J., Lyne A. G., Mickaliger M. B., Stappers B. W., Turner J. D., WeltevredeP., 2021, Monthly Notices of the Royal Astronomical Society: Letters, 505, L6
* Speagle J. S., 2019, arXiv:1904.02180 [astro-ph, stat 10/ghjwgz
