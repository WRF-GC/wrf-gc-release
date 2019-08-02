# WRF-GC | A two-way coupled meteo/chem model

This repository (`wrf-gigc`) holds the `chem` directory for the WRF-GC model, which contains all the custom code necessary to run WRF-GC under a WRF installation
and a built-in version of GEOS-Chem (currently `12.2.1` official version)

**This version of WRF-GC is version:** `1.0`

(c) 2017-2019 Haipeng Lin <linhaipeng@pku.edu.cn>, Xu Feng <fengx7@pku.edu.cn>, Tzung-May Fu <tmfu@pku.edu.cn>

(c) 2017-2019 Atmospheric Chemistry & Climate Group, Peking University

GEOS-Chem, GEOS-Chem High Performance, HEMCO, ESMF/MAPL Frameworks are (c) their original authors.

## Maintainers
Correspondence To: Tzung-May Fu (`fuzm at sustech.edu.cn`)

Atmospheric Chemistry & Climate Group, Southern University of Science and Technology (SUSTECH)

* Haipeng Lin (`hplin at g.harvard.edu`) - Code Architectural Design & Maintenance, Lead Developer
* Xu Feng (`fengx7 at pku.edu.cn`) - Scientific Programming

## Primary Components for Coupling

* WRF-to-Chemistry Abstraction Layer (Codename "Pumpkin") - `wrf-gchp-pumpkin` Project
* WRF-Grid-Independent-GEOS-Chem ("GIGC") Chemistry Driver - `chemics_init` & `chem_driver`
* Stateful Conversion Module - `GIGC_Convert_State_Mod`

## Support Components
The below support components are based off GEOS-Chem High Performance ("GCHP") technology.

* GIGC Column Code - `GIGC_Chunk_Mod`, based off GCHP's original column code
* GIGC Stateful Variables Container - `GIGC_Stateful_Mod`, redesigned to support multiple domains within the same CPU in GEOS-Chem

## License
```
(c) 2017-2019 Haipeng Lin <hplin@g.harvard.edu>, Xu Feng <fengx7@pku.edu.cn>, Tzung-May Fu <fuzm@sustech.edu.cn>
(c) 2017-2019 Atmospheric Chemistry & Climate Group, Southern University of Science and Technology

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to 
use, copy, modify the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

- The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

- The Software, modified in part or in full may not be redistributed without
express permission from the copyright holder(s).

Except as contained in this notice or in attribution, the name of the WRF-GC model
shall not be used as an endorsement for distributing modified copies of the
Software without prior written permission from the copyright holder(s).

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```