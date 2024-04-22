# README for the HEMCO source code repository

[![Latest Pre-Release](https://img.shields.io/github/v/release/geoschem/HEMCO?include_prereleases&label=Latest%20Pre-Release)](https://github.com/geoschem/HEMCO/releases) [![Latest Stable Release](https://img.shields.io/github/v/release/geoschem/HEMCO?label=Latest%20Stable%20Release)](https://github.com/geoschem/HEMCO/releases) [![Release Date](https://img.shields.io/github/release-date/geoschem/HEMCO)](https://github.com/geoschem/HEMCO/releases) [![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/geoschem/HEMCO/blob/main/LICENSE.txt) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4618253.svg)](https://doi.org/10.5281/zenodo.4618253) [![RTD](https://img.shields.io/readthedocs/hemco?label=ReadTheDocs)](https://hemco.readthedocs.io/en/latest/) [![Azure](https://dev.azure.com/geoschem/hemco/_apis/build/status/Quick%20Build?branchName=dev)](https://dev.azure.com/geoschem/hemco/_apis/build/status/Quick%20Build?branchName=dev)

This repository (https://github.com/geoschem/HEMCO) contains the Harmonized Emissions Component
(HEMCO) source code. HEMCO is a software component for computing (atmospheric) emissions from
different sources, regions, and species on a user-defined grid. It can combine, overlay, and
update a set of data inventories ('base emissions') and scale factors, as specified by the user
through the HEMCO configuration file. Emissions that depend on environmental variables and
non-linear parameterizations are calculated in separate HEMCO extensions. HEMCO can be run
in standalone mode or coupled to an atmospheric model. A more detailed description of HEMCO
is given in Keller et al. (2014) and Lin et al (2021).

HEMCO has been coupled to several atmospheric and Earth System Models, and can be coupled with
or without using the Earth System Modeling Framework (ESMF). A detailed description of HEMCO
coupled with other models is given in Lin et al. (2021).

## Documentation

### Reference

C. A. Keller, M. S. Long, R. M. Yantosca, A. M. Da Silva, S. Pawson, D. J. Jacob, *HEMCO v1.0: a versatile,
ESMF-compliant component for calculation emissions in atmospheric models*, <u>Geosci. Model Dev.</u>, **7**, 1409-1417, 2014.

Lin, H., Jacob, D. J., Lundgren, E. W., Sulprizio, M. P., Keller, C. A., Fritz, T. M., Eastham, S. D., Emmons, L. K., Campbell, P. C., Baker, B., Saylor, R. D., and Montuoro, R.: *Harmonized Emissions Component (HEMCO) 3.0 as a versatile emissions component for atmospheric models: application in the GEOS-Chem, NASA GEOS, WRF-GC, CESM2, NOAA GEFS-Aerosol, and NOAA UFS models*, <u>Geosci. Model Dev.</u>, **14**, 5487–5506, 2021.

### Online user's manual

Installation and usage instructions are posted online at [hemco.readthedocs.io](http://hemco.readthedocs.io)

## Support
We encourage GEOS-Chem users to use [the Github issue tracker attached to this repository](https://github.com/geoschem/HEMCO/issues/new/choose) to report bugs or technical issues with the HEMCO code.

## License

HEMCO is distributed under the MIT license. Please see the license documents LICENSE.txt and
AUTHORS.txt in the root folder.