# grins

A free and open source Java library for parsing and generating [Simplified molecular-input line-entry system (SMILES)](http://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) line notations.

### Objectives

The primary goal of the library is to provide an implementation of the OpenSMILES specification for the [Chemistry Development Kit (CDK)](http://sourceforge.net/projects/cdk/). The library is separate to allow an object model that precisely represents SMILES structures. In addition there are several design objectives the library aims to fulfil.

#### Round tripping

Parsing and generating a _standard form_ structure without modifying it should produce the exact same SMILES.

#### External Canonicalisation

Canonization will be achieved by providing a permutation on the vertices of the chemical graph. This allows external procedures such as the InChI numberings of Universal SMILES ([O’Boyle, 2012](http://www.jcheminf.com/content/4/1/22)) to be easily applied.

#### Correct handling of stereo chemistry

Persistent stereo configuration of tetrahedral and double bond configurations (likely to extend to other forms at a later date).

#### Validation

Report errors and inconsistencies with SMILES input.

