# Grins

A free and open source Java toolkit dedicated to parsing and generating [Simplified molecular-input line-entry system (SMILES)](http://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) line notations. The primary focus of the library is to handle the SMILES syntax and to provide an implementation of the [OpenSMILES](http://www.opensmiles.org) specification for the [Chemistry Development Kit (CDK)](http://sourceforge.net/projects/cdk/). The code is separate to avoid cluttering the CDK API with data structures specific to Grins.

## Goals
 - Handle tetrahedral and double bond stereo chemistry and maintain configuration through different orderings
 - Canonization via a permutation of vertices in the chemical graph
 - Unification of double bond stereo using atom-based configuration - ['@' for Cis/Trans around Double Bonds
](http://www.opensmiles.org/opensmiles.html#_tt_tt_for_cis_trans_around_double_bonds)

## Planned 
 - Generation of Universal SMILES ([O’Boyle, 2012](http://www.jcheminf.com/content/4/1/22))
 - Allene, Square Planar, Trigongal Bipyramidal and Octahedral topologies

## Examples

_grins is still in a development and some APIs will likely change until a release is made._

The main 'molecule' class in _grins_ is the 'ChemicalGraph' it provides convenience methods for reading SMILES directly.

```java
ChemicalGraph g = ChemicalGraph.fromSmiles("CCO");
```

and for writing it back to SMILES notation.

```java
String smi = g.toSmiles();
```

Grins provides excellent round tripping, preserving exactly how the input was specified. Disregarding inputs with redundant brackets and erroneous/repeated ring numbers - the actually input will generally be identical to the output.

```java
// bond labels
ChemicalGraph.fromSmiles("C1=CC=CC=C1").toSmiles();    // kekule      (implicit single bonds)
ChemicalGraph.fromSmiles("C-1=C-C=C-C=C1").toSmiles(); // kekule      (explicit single bonds)
ChemicalGraph.fromSmiles("c1ccccc1").toSmiles();       // delocalised (implicit aromatic bonds)
ChemicalGraph.fromSmiles("c:1:c:c:c:c:c1").toSmiles(); // delocalised (explicit aromatic bonds)

// bracket atoms stay as bracket atoms
ChemicalGraph.fromSmiles("[CH]1=[CH][CH]=[CH][CH]=[CH]1").toSmiles();
ChemicalGraph.fromSmiles("[CH]1=[CH]C=C[CH]=[CH]1").toSmiles();       // mix bracket and subset atoms
```

Although preserving the representation was one of the design goals for grins it is common practise to normalise.
As such there are several utilities for standardising outputs.

One can easily convert a molecule with implicit hydrogens `[CH3][CH2][OH]` to
one with inferred hydrogens `CCO`.

```java
ChemicalGraph g = ChemicalGraph.fromSmiles("[CH3][CH2][OH]");
ChemicalGraph h = Functions.collapse(g);
h.toSmiles().equals("CCO");
```

Likewise, conversion of a molecule with inferred hydrogens `CCO` to
one with implicit hydrogens `[CH3][CH2][OH]` is also easy.

```java
ChemicalGraph g = ChemicalGraph.fromSmiles("CCO");
ChemicalGraph h = Functions.expand(g);
h.toSmiles().equals("[CH3][CH2][OH]");
```

Randomly generate different SMILES notations preserving stereo-configuration.

```java
ChemicalGraph g  = ChemicalGraph.fromSmiles("CCC[C@@](C)(O)[C@H](C)N");
StringBuilder sb = new StringBuilder(g.toSmiles());
for (int i = 0; i < 25; i++)
    sb.append('.').append(Functions.randomise(g).toSmiles());
System.out.println(sb);
```

Using atom-based double-bond configuration.

```java
ChemicalGraph g   = ChemicalGraph.fromSmiles("F/C=C/F");
ChemicalGraph h   = Functions.atomBasedDBStereo(g);
String        smi = h.toSmiles(); // F[C@H]=[C@@H]F
```

Normalise directional labels. There are two possible ways to write each
bond-based configuration. Grins allows you to normalise the labels such
that the first symbol is always a forward slash (`/`). Some examples are
shown below.

```java
ChemicalGraph g   = ChemicalGraph.fromSmiles("F\\C=C/F");
ChemicalGraph h   = Functions.normaliseDirectionalLabels(g);
String        smi = h.toSmiles();
```

<table>
<tr><th>Original</th>               <th>Normalised</th>
<tr><td>`F/C=C/F`</td>              <td>`F/C=C/F`</td>
<tr><td>`F\C=C\F`</td>              <td>`F/C=C/F`</td>
<tr><td>`F/C=C\F`</td>              <td>`F/C=C\F`</td>
<tr><td>`F\C=C/F`</td>              <td>`F/=C\F`</td>
<tr><td>`C(\F)(/C)=C\F`</td>        <td>`C(/F)(\C)=C/F`</td>
<tr><td>`C(\\F)=C\\F`</td>          <td>`C(/F)=C(\F)`</td>
<tr><td>`FC=C(F)C=C(F)\\C=C\\F`</td><td>`FC=C(F)C=C(F)/C=C/F`</td>
</table>

## How to Grin

Grins is still in development but you can obtain the latest build from the [EBI snapshots repository](http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo-snapshots/). An example configuration for maven is shown below.

```xml
<project>
...
<repositories>
   <repository>
      <id>ebi-repo</id>
      <url>http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo/</url>
   </repository>
   <repository>
      <id>ebi-repo-snapshots</id>
      <url>http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo-snapshots/</url>
   </repository>
</repositories>
...
<dependencies>
    <dependency>
        <groupId>uk.ac.ebi.grins</groupId>
        <artifactId>grins-core</artifactId>
        <version>LATEST</version>
    </dependency>
    <dependency>
        <groupId>uk.ac.ebi.grins</groupId>
        <artifactId>grins-func</artifactId>
        <version>LATEST</version>
    </dependency>
</dependencies>
...
</project>
```

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