********
Examples
********


Running leavitt
===============
The simplest way to run |Leavitt| is with the command-line script ``leavitt``.  The only required argument is the name of a spectrum fits file.

.. code-block:: bash

    leavitt catalog.fits

This will create an output file called ``catalog_leavitt.fits``.

By default, |Leavitt| doesn't print anything to the screen.  So let's set the ``--verbose`` or ``-v`` parameter.

.. code-block:: bash
		
    leavitt catalog.fits -v

You can also have it save a figure of the best-fitting model using the ``--plot`` or ``-p`` parameter.

.. code-block:: bash
		
    leavitt catalog.fits -v -p
    ...
    chisq =  1.96
    Figure saved to catalog_leavitt.png
    dt =  4.08 sec.
    Writing output to catalog_leavitt.fits

By default, the figure filename is ``catalog_leavitt.png``, but the figure filename can be set with the
``--figfile`` parameter.  Note that you can set the figure type with the extension.  For example,
to get a PDF just use ``--figfile catalog_leavitt.pdf``.

    
Running |Leavitt| from python
=============================

All of the |Leavitt| functionality is also available directory from python.

First, import |Leavitt| and read in the example catalog.

.. code-block:: python

	import leavitt
	cat = leavitt.read('catalog.fits')

Now fit the time-series data.

.. code-block:: python
		
	out,model = leavitt.fit(cat)

The output will be a table with the final results, the best-fitting model
