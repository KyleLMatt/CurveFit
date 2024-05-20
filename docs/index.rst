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


.. rubric:: Footnotes

.. [#f1] For `Henrietta Leavitt <https://en.wikipedia.org/wiki/Henrietta_Swan_Leavitt>`_ who was american astronomer that discovered the important period-luminosity relationship of variable stars.

