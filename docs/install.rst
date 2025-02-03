Installation
============

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


Then install the package:

.. parsed-literal::
   python -m pip install . --no-deps

Tests
-----

First install test dependencies:

.. parsed-literal::
   micromamba install -c conda-forge pytest pytest-xdist

Now run the tests:

.. parsed-literal::
   pytest -n auto -v
