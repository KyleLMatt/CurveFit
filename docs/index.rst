.. leavitt documentation master file, created by
   sphinx-quickstart on Wed Sep 22 13:03:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**********
Leavitt
**********

Introduction
============
|Leavitt| [#f1]_ is a general-purpose variable star fitting software.  The goal is to provide easy-to-use variale star software includig robust classification.

.. toctree::
   :maxdepth: 1

   install
   gettingstarted
   modules   

Description
===========
|Leavitt| has a variety of variable star tools.

- Download NSC time series data
- Generate multi-band Lomb-Scargle periodogram
- Period estimation using the Saha and Vivas (2017) hybrid method.
- Generic template fitting (each band separately) using the Layden templates.
- Template fitting using "self template".
- RR Lyrae-specifc template fitting using RRL templates and holding band-to-band amplitude ratios to the RRL values.
- Generate Cepheid light curve using the Yoachim et al. (2009) templates 

.. rubric:: Footnotes

.. [#f1] For `Henrietta Leavitt <https://en.wikipedia.org/wiki/Henrietta_Swan_Leavitt>`_ who was american astronomer that discovered the important period-luminosity relationship of variable stars.

