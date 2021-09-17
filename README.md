# WRF-GC | A two-way coupled meteo/chem model

This repository (`wrf-gc`) holds the `chem` directory for the WRF-GC model, which contains all the custom code necessary to run WRF-GC under a WRF installation
and a built-in version of GEOS-Chem (currently `12.8.2` official version)

**This version of WRF-GC is version:** `2.0.1-release` (2021/09/17)

(c) 2017-2021 Haipeng Lin <hplin@seas.harvard.edu>, Xu Feng <fengx7@pku.edu.cn>, Tzung-May Fu <fuzm@sustech.edu.cn>

(c) 2017-2021 Atmospheric Chemistry & Climate Group, SUSTech

GEOS-Chem, GEOS-Chem High Performance, HEMCO, ESMF/MAPL Frameworks are (c) their original authors.

## Maintainers
Correspondence To: Tzung-May Fu (`fuzm at sustech.edu.cn`)

Atmospheric Chemistry & Climate Group, SUSTech

### Two-way coupled version (`2.0`)
* Xu Feng (`fengx7 at pku.edu.cn`) - Two-way architecture design and lead development
* Haipeng Lin (`hplin at seas.harvard.edu`) - Nested-grid model development, code maintenance

### One-way coupled version (`0.1` - `1.2`)
* Haipeng Lin (`hplin at seas.harvard.edu`) - Code Architectural Design & Maintenance, lead development
* Xu Feng (`fengx7 at pku.edu.cn`) - Scientific programming

## Primary Components for Coupling

### One-way coupled version, two-way coupled version
* WRF-to-Chemistry Abstraction Layer (Codename "Pumpkin") - `wrf-gchp-pumpkin` Project
* WRF-Grid-Independent-GEOS-Chem ("GIGC") Chemistry Driver - `chemics_init` & `chem_driver`
* Stateful Conversion Module - `WRFGC_Convert_State_Mod` (formerly `GIGC_Convert_State_Mod`)

### Nested model
* State management module - `GC_Stateful_Mod`

## Support Components
The below support components are based off GEOS-Chem High Performance ("GCHP") technology.

* GEOS-Chem Column Code - `GIGC_Chunk_Mod`, based off GCHP's original column code
* GEOS-Chem Stateful Variables Container - `GC_Stateful_Mod`, redesigned to support multiple domains within the same CPU in GEOS-Chem

## License
```
(c) 2017-2020 Haipeng Lin <hplin@seas.harvard.edu>, Xu Feng <fengx7@pku.edu.cn>, Tzung-May Fu <tmfu@pku.edu.cn>
(c) 2017-2018 Atmospheric Chemistry & Climate Group, Peking University
(c) 2018-2021 Atmospheric Chemistry & Climate Group, SUSTech

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
