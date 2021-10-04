.. leavitt documentation master file, created by
   sphinx-quickstart on Wed Sep 22 13:03:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**********
Leavitt
**********

Introduction
============
|Leavitt| [#f1]_ is a general-purpose variable star fitting software.

.. toctree::
   :maxdepth: 1

   install
   gettingstarted
   examples
   

Description
===========
|Leavitt| fits variable star light-curves using a multi-step approach.

1. Period estimation using the Saha and Vivas (2017) hybrid method.
2. Run "Upsilon" variable star random forest classification software (only for 80+ measurements).
3. Template fitting using "self template".
4. Generic template fitting (each band separately) using the Layden templates.
5. RR Lyrae-specifc template using RRL templates and holding band-to-band amplitude ratios to the RRL values.
6. Cepheid-specific template fitting using period-dependent light curves and olding band-to-band amplitude ratios to the Cepheid values.
7. Eclipsing Binary fitting using EBAI.
8. Classification using all above information.
9. Estimate distance using period-luminosity relationship for the classified type.

|Leavitt| can be called from python directly or the command-line script ``leavitt`` can be used.


Examples
========

.. toctree::
    :maxdepth: 1

    examples
    gettingstarted


Leavitt
=======
Here are the various input arguments for command-line script ``leavitt``::

  usage: doppler [-h] [--outfile OUTFILE] [--payne] [--fitpars FITPARS]
                 [--fixpars FIXPARS] [--figfile FIGFILE] [-d OUTDIR] [-l] [-j]
                 [--snrcut SNRCUT] [-p] [-c] [-m] [-r READER] [-v]
                 [-nth NTHREADS] [--notweak] [--tpoly] [--tpolyorder TPOLYORDER]
                 files [files ...]

  Run Doppler fitting on spectra

  positional arguments:
    files                 Spectrum FITS files or list

  optional arguments:
    -h, --help            show this help message and exit
    --outfile OUTFILE     Output filename
    --payne               Fit a Payne model
    --fitpars FITPARS     Payne labels to fit (e.g. TEFF,LOGG,FE_H
    --fixpars FIXPARS     Payne labels to hold fixed (e.g. TEFF:5500,LOGG:2.3
    --figfile FIGFILE     Figure filename
    -d OUTDIR, --outdir OUTDIR
                          Output directory
    -l, --list            Input is a list of FITS files
    -j, --joint           Joint fit all the spectra
    --snrcut SNRCUT       S/N threshold to fit spectrum separately
    -p, --plot            Save the plots
    -c, --corner          Make corner plot with MCMC results
    -m, --mcmc            Run MCMC when fitting spectra individually
    -r READER, --reader READER
                          The spectral reader to use
    -v, --verbose         Verbose output
    -nth NTHREADS, --nthreads NTHREADS
                          Number of threads to use
    --notweak             Do not tweak the continuum using the model
    --tpoly               Use low-order polynomial for tweaking
    --tpolyorder TPOLYORDER
                          Polynomial order to use for tweaking

			  
.. rubric:: Footnotes

.. [#f1] For `Henrietta Leavitt <https://en.wikipedia.org/wiki/Henrietta_Swan_Leavitt>`_ who was american astronomer that discovered the important period-luminosity relationship of variable stars.

