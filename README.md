[![DOI](https://zenodo.org/badge/443812727.svg)](https://zenodo.org/doi/10.5281/zenodo.13354423)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Workflow Status](https://github.com/mdcourse/mdcourse.github.io/actions/workflows/tests.yml/badge.svg)
![Workflow Status](https://github.com/mdcourse/mdcourse.github.io/actions/workflows/gh-pages.yml/badge.svg)

# Learn Molecular Simulations with Python

<img src="docs/source/_static/logo/logo-py.png" width="30%" align="right"/></a>

The goal of [Learn Molecular Simulations with Python](https://mdcourse.github.io/)
is to write a simple code containing most of the basic functionalities of molecular
simulations, such as:
- Energy minimization,
- Molecular dynamics,
- Monte Carlo move.

The Python code that is written here is used to realize molecular
scientific projects. Note that the code is slow and that efficiency is not the
primary objective here.

### Contributing

We welcome contributions from the community. Before you start, please read our
[Contributing Guidelines](./CONTRIBUTING.md). These guidelines will help you
understand the process and expectations for contributing to the project.

### Prerequisite

The target audience includes people who are either completely new to molecular simulations
or users of open-source codes such as LAMMPS and GROMACS who want to better understand
what is behind those codes. Although some basic knowledge of coding, thermodynamics, and
statistical physics is recommended for a full understanding of molecular simulations,
[Learn Molecular Simulations with Python](https://mdcourse.github.io/) can be followed
even without deep expertise in these fields. Annexes with key concepts and suggested
readings are provided as needed.

### What is not (yet) in the code

- molecules/shake,
- electrostatic
- Monte Carlo
- thermostats and barostats other than Berendsen,
- energy minimization methods other than the steepest descent.

### Authors

#### Project creator

- [Simon Gravelle](https://simongravelle.github.io/): computer physicist in soft matter
  and fluids at interfaces.

#### Contributors

-  [Jake Anderson](https://github.com/jaketanderson) 
-  [Kenneth Ngo](https://github.com/KennethQNgo)

### Acknowledgments and license

This project has received funding from the European
Union's Horizon 2020 research and innovation programme
under the Marie Skłodowska-Curie grant agreement No 101065060.
All inputs, scripts, and data files are released under the GNU
General Public License v3.0. The released files have been uploaded
to [Zenodo](https://zenodo.org/), under its [DOI](https://zenodo.org/records/13624976).

### TODO

- Choose a different title? (Build Your Own Molecular Simulations with Python,
  Create Molecular Simulations with Python, Construct Molecular Simulations with Python)
- Find a cool acronym. MS-PyLearn, LMS-py ?
