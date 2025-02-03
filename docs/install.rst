Installation
============

First clone the repository and ``cd`` into the repository root folder:

.. parsed-literal::

   git clone https://github.com/omsf/OpenPharmMDFlow.git
   cd OpenPharmMDFlow

.. warning::
   This software is still in early developlment

Install dependencies using Micromamba_:

.. _Micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

.. parsed-literal::
   micromamba create -n openpharmmdflow -c conda-forge --file conda-envs/macos-latest.yaml

or

.. parsed-literal::
   micromamba create -n openpharmmdflow -c conda-forge --file conda-envs/ubuntu-latest.yaml

depending on your operating system.

.. note::
   If you run into network errors, you may need to configure a proxy.

Then activate the environment and install the package:

.. parsed-literal::
   micromamba activate openpharmmdflow
   python -m pip install . --no-deps

Tests
-----

First install test dependencies:

.. parsed-literal::
   micromamba install -c conda-forge pytest pytest-xdist

Now run the tests:

.. parsed-literal::
   pytest -n auto -v
