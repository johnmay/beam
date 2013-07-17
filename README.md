# Grins

A free and open source Java library for parsing and generating [Simplified molecular-input line-entry system (SMILES)](http://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) line notations. The primary focus of the library is to provide an implementation of the [OpenSMILES](http://www.opensmiles.org) specification for the [Chemistry Development Kit (CDK)](http://sourceforge.net/projects/cdk/). The code is separate to avoid cluttering the CDK API with data structures specific to Grins.

## Goals
 - Handle tetrahedral and double bond stereo chemistry and maintain configuration through different orderings
 - Canonization via a permutation of vertices in the chemical graph
 - Unification of double bond stereo using atom-based configuration - ['@' for Cis/Trans around Double Bonds
](http://www.opensmiles.org/opensmiles.html#_tt_tt_for_cis_trans_around_double_bonds)

## Planned 
 - Generation of Universal SMILES ([O’Boyle, 2012](http://www.jcheminf.com/content/4/1/22))
 - Allene, Square Planar, Trigongal Bipyramidal and Octahedral topologies

## License BSD 2-Clause

Copyright (c) 2013, European Bioinformatics Institute (EMBL-EBI)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the authors and should not be interpreted as representing official policies, either expressed or implied, of the FreeBSD Project.

## Disclaimer

This software is not associated or endorsed by [Daylight Chemical Information Systems](http://www.daylight.com) and is in no way related to their discontinued GRaphical INput for Smiles (GRINS) chemical editors (xvgrins, cgi-grins and JavaGrins).