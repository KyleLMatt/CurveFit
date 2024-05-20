*****************************
Getting Started with Variable
*****************************



Imports
=======

.. code-block:: python

	from leavitt.leavitt.timeseries import Variable
	from leavitt.leavitt.utils import phase_fold, plot_phased_lightcurve


There are two main ways to initialize a ```Variable``` object. 
1. Using the NSC ID (objectid) of a star. We'll do that here. This performs a synchronous query to the DataLab catalog.
2. Giving it a TimeSeries object.

.. code-block:: python
   
	star = Variable('150537_4644')

This object now has an attribute called ```timeseries``` which contains the timeseries data.

.. code-block:: python

	star.timeseries

Now, we can perform a Lomb-Scargle multiband to find the period of the star (it may take a while). This is a RRc star with a period of 0.3367 days.

.. code-block:: python

	frequency, power = star.ls_mb_periodogram()
	period, error = star.get_period(frequency, power)

The period and errors are stored in the new variables, but also a new attribute called ```period``` that stores the period only.

.. code-block:: python

	print('{:.5f} +/- {:.5f} days'.format(period.value,error.value))
	print('{:.5f} days'.format(star.period.value))

We can now get the data to construct a phased lightcurve.

.. code-block:: python

	phase = phase_fold(star.timeseries['time'],period)

	plot_phased_lightcurve(phase, star.timeseries['mag_auto'],mags_errs=star.timeseries['magerr_auto'],filters=star.timeseries['filter'])

