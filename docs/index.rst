PySME documentation
===================

More than two decades ago `Valenti & Piskunov (1996) <https://ui.adsabs.harvard.edu/abs/1996A&AS..118..595V>`_ developed SME - Spectroscopy Made Easy, a high-precision stellar-spectra synthesis/analysis engine that has powered hundreds of studies.
PySME is its modern Python front-end: a wrapper around the original C++/Fortran core that lets you (1) compute accurate, high-resolution synthetic spectra from a linelist + model atmosphere, (2) invert observed spectra to derive stellar parameters, and (3) explore NLTE corrections — all from an interactive notebook or scripted pipeline. The same capabilities make PySME invaluable for exoplanet work, where characterising the host star is essential for understanding its planets.

🎉 PySME (version 0.6.14+) is now fully compatible with the arm64 MacOS! Feel free to try it out. 🎉

.. admonition:: Key features

   * Plane-parallel and spherical radiative-transfer engine  
   * LTE & 1-D NLTE line formation with pre-computed grids  
   * Automatic :math:`\chi^2` fitting for :math:`T_\mathrm{eff}`, :math:`\log{g}`, :math:`v_\mathrm{mic}`, [X/Fe] …  
   * Seamless use of ATLAS and MARCS model atmospheres and VALD line lists

.. raw:: html

   <section class="home-cards">
     <a class="home-card" href="getting_started/installation.html">
       <h3>Installation</h3>
       <p>Set up PySME and verify your code before synthesizing.</p>
       <span class="home-card-cta">Go to installation</span>
     </a>
     <a class="home-card" href="getting_started/index.html">
       <h3>Getting started</h3>
       <p>Run the first spectrum and first fit with a minimal, beginner-friendly path.</p>
       <span class="home-card-cta">Open getting started</span>
     </a>
     <a class="home-card" href="advance/index.html">
       <h3>Advanced usage</h3>
       <p>Learn core usage and structure details for deeper and more controlled use.</p>
       <span class="home-card-cta">To advanced usage</span>
     </a>
   </section>

.. toctree::
   :maxdepth: 1

   getting_started/index
   fundamentals/index
   advance/index
   concepts/index
   dev/index

Citation
~~~~~~~~

- Jian et al. (2026; in prep.)
- `Wehrhahn et al. (2023) <https://ui.adsabs.harvard.edu/abs/2023A&A...671A.171W>`_

Indices and tables
~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. rubric:: Quick links
:GitHub repository: https://github.com/MingjieJian/SME
:Issue tracker:     https://github.com/MingjieJian/SME/issues
