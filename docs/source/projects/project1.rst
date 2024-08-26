.. _project1-label:

A study of argon
================

Using Monte Carlo simulations, the equation of state of 3D fluid of argon is simulated
for a large range of density.

.. figure:: project1/avatar-dm.webp
    :alt: The fluid made of argon atoms and simulated using monte carlo and python.
    :height: 200
    :align: right
    :class: only-dark

.. figure:: project1/avatar.webp
    :alt: The fluid made of argon atoms and simulated using monte carlo and python
    :height: 200
    :align: right
    :class: only-light

In 1957, Wood and Parker reported the first numerical study of a 3D fluid that was
using an attractive Lennard-Jones potential (previous studies by Metropolis et al.
only used repultive hard-sphere potentials :cite:`metropolis1953equation, metropolis1953simulated`).
In their study, Wood and Parker used a Monte Carlo algorithm to predit the 
Equation of State (EoS) of neutral particles interacting with Lennard-Jones potential,
whose parameters were chosen to match those of argon gas. Their results show
good agreement with experimental measurements of Michels on argon :cite:`michels1949isotherms`.

Here, we take advantage of our code and try to reproduce the results of Wood and Parker
for varying density. To follow this project, only Monte Carlo move are needed,
and all the chapters up to :ref:`chapter7-label` must have been followed.

Parameters choice
-----------------

Some parameters are taken *exactly* the same as those from Wood and Parker:

- :math:`\sigma = 3.405~\text{Å}` (calculated as :math:`\sigma = r^*  2^{-1/6}`
  where :math:`r^*  = 3.822~\text{Å}`),
- :math:`\epsilon = 0.238~\text{kcal/mol}`,
- :math:`m = 39.948~\text{g/mol}`,
- :math:`T = 328.15~\text{K}` (or :math:`T = 55~^\circ\text{C}`).

In addition, some parameters that are not specified by Wood and Parker are
freely chosen, but it should not impact the result too much, except may be
for the cut-off:

- :math:`d_\text{mc} = \sigma/5`,
- :math:`r_\text{c} = 2.5 \sigma`.

Finally, since an average modern laptop is much faster than Wood and Parker's
IBM 701 calculators (if you are curious, here is the |IBM-wiki| page of the IBM 701),
let us choose a number of atoms that is slightly larger
than theirs, i.e. :math:`N_\text{atom} = 200`.

.. |IBM-wiki| raw:: html

    <a href="https://en.wikipedia.org/wiki/IBM_701" target="_blank">wikipedia</a>

In order to cover the same density range as Wood and Parker, let us vary the volume
of the box, :math:`V`, so that :math:`V/V^* \in [0.75, 7.6]`, where
:math:`V^* = 2^{-1/2} N_\text{A} r^{*3} = 23.77 \text{centimeter}^3/\text{mole}`
is the molar volume.

Equation of State
-----------------

The Equation of State (EoS) is a fundamental relationship in thermodynamics and
statistical mechanics that describes how the state variables of a system, such as
pressure :math:`p`, volume :math:`V`, and temperature :math:`T`, are interrelated.

Here, let us extract the pressure of a simple fluid for different density values.