# Licensing

## Original project contributions: BSD-3-Clause

Copyright (c) 2026 Zhong-Xuan Mao and Jing-Feng Tian. The original SE-RR kernels
(`src/python/marcum/serr.py`, `src/cpp/serr.cpp`, `src/fortran/serr.f90`) carry
BSD-3-Clause SPDX notices. The grant also covers the authors' independently
authored contributions. See the root LICENSE.md and LICENSES/BSD-3-Clause.txt.
Please cite https://arxiv.org/abs/2611.29209 when using this work in research;
this citation request is not an additional BSD license condition.

## Algorithm 939: ACM terms

The unchanged upstream archive/extracted source, Python/C++ component ports,
RR/LAE bridges and fixed-work adaptations retain the ACM terms. Project headers
do not transfer upstream ownership or license these components under BSD.
The adapted files are modified works, not the original Algorithm 939 software.
See third_party/algorithm939/LICENSE-ACM.txt and LICENSE-SOURCE.md.
The ACM agreement permits academic/research/noncommercial use subject to its
conditions, including notice retention, no charge for modified distributions,
identifying modifications and notification to ACM and the original authors.
Commercial use requires a separate ACM license.

## GammaCHI 2021: GPLv2

The bundled GammaCHI source identifies the 2021 version. Its official software
record lists GPLv2: https://data.mendeley.com/datasets/d2kwwvsyny/1 .
The full GPLv2 text is included in third_party/gammachi/COPYING. Original author
comments and README are retained. GammaCHI and works derived from or linked
with it must satisfy the applicable GPLv2 terms, including source availability
when distributing binaries. BSD-3-Clause contributions can be included in a
GPLv2 work without changing the upstream license. Separate source components
in this repository are not all relicensed as one blanket BSD/GPL package.

The GammaCHI application is built separately from the Algorithm 939 baseline;
the repository distributes source rather than precompiled combined binaries.
Any new integration must respect both applicable upstream terms.

## References

- BSD-3-Clause: https://opensource.org/license/bsd-3-clause
- CALGO: https://calgo.acm.org/
- ACM agreement mirror: https://mipav.cit.nih.gov/license/CopyrightInfo.html
- GammaCHI 2021: https://doi.org/10.17632/d2kwwvsyny.1
- GPLv2: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
