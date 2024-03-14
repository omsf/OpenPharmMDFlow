Installation
============

.. warning::
   This software is still in early developlment

Install dependencies using micromamba:

.. parsed-literal::
   micromamba create -n openpharmmdflow --file examples/mini_mvp/env.yaml

Then install the package:

.. parsed-literal::
   python -m pip . --no-deps

Tests
-----

First install test dependencies:

.. parsed-literal::
   micromamba install pytest pytest-xdist

Now run the tests:

.. parsed-literal::
   pytest -n auto -v
