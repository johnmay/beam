# Beam

Beam - _to express by means of a radiant smile_ 

Beam is a free toolkit dedicated to parsing and generating Simplified
molecular-input line-entry system - [SMILES&trade;](http://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
line notations. The primary focus of the library is to elegantly handle the
SMILES&trade; syntax.

## Beaming

*Note: Beam is still in a development and some APIs will likely change until a release is made.*

One of the primary types in Beam is the `Graph` it provides convenience
methods for reading SMILES&trade; notation directly.

```java
Graph g = Graph.fromSmiles("CCO");
```

and for writing it back to SMILES&trade; notation.

```java
String smi = g.toSmiles();
```

Beam provides excellent round tripping, preserving exactly how the input was
specified. Disregarding inputs with redundant brackets and erroneous/repeated
ring numbers - the actually input will generally be identical to the output.

```java
// bond labels
Graph.fromSmiles("C1=CC=CC=C1").toSmiles();    // kekule      (implicit single bonds)
Graph.fromSmiles("C-1=C-C=C-C=C1").toSmiles(); // kekule      (explicit single bonds)
Graph.fromSmiles("c1ccccc1").toSmiles();       // delocalised (implicit aromatic bonds)
Graph.fromSmiles("c:1:c:c:c:c:c1").toSmiles(); // delocalised (explicit aromatic bonds)

// bracket atoms stay as bracket atoms
Graph.fromSmiles("[CH]1=[CH][CH]=[CH][CH]=[CH]1").toSmiles();
Graph.fromSmiles("[CH]1=[CH]C=C[CH]=[CH]1").toSmiles();       // mix bracket and subset atoms
```

Although preserving the representation was one of the design goals for beam it
is common to normalise output SMILES&trade;.

_Collapse_ a graph with labelled hydrogens `[CH3][CH2][OH]` to one with implicit
hydrogens `CCO`.

```java
Graph g = Graph.fromSmiles("[CH3][CH2][OH]");
Graph h = Functions.collapse(g);
h.toSmiles().equals("CCO");
```

_Expand_ a graph where the hydrogens are implicit `CCO` to one with labelled
hydrogens `[CH3][CH2][OH]`.

```java
Graph g = Graph.fromSmiles("CCO");
Graph h = Functions.expand(g);
h.toSmiles().equals("[CH3][CH2][OH]");
```

Stereo specification is persevered through rearrangements. The example below 
randomly generates arbitrary SMILES&trade; preserving correct stereo-configuration.

```java
Graph g  = Graph.fromSmiles("CCC[C@@](C)(O)[C@H](C)N");
StringBuilder sb = new StringBuilder(g.toSmiles());
for (int i = 0; i < 25; i++)
    sb.append('.').append(Functions.randomise(g).toSmiles());
System.out.println(sb);
```

Bond based double-bond configuration is normal is SMILES but can be problematic.
The issue is that a single symbol may be specifying two adjacent configurations.
A proposed extension was to use atom-based double-bond configuration.

Beam will input, output and convert atom and bond-based double-bond stereo 
specification. 

```java
Graph  g   = Graph.fromSmiles("F/C=C/F");
Graph  h   = Functions.atomBasedDBStereo(g);
String smi = h.toSmiles();
smi.equals("F[C@H]=[C@@H]F");
```

```java
Graph  g   = Graph.fromSmiles("F[C@H]=[C@@H]F");
Graph  h   = Functions.bondBasedDBStereo(g);
String smi = h.toSmiles();
smi.equals("F/C=C/F");
```

Convert a graph with delocalised bonds to kekul&eacute; representation.

```java
Graph  furan        = Graph.fromSmiles("o1cccc1");
Graph  furan_kekule = furan.kekule();
String smi          = furan_kekule.toSmiles();
smi.equals("O1C=CC=C1");
```

With bond-based double-bond stereo specification there are two possible ways to
write each bond-based configuration. beam allows you to normalise the labels such
that the first symbol is always a forward slash (`/`). Some examples are shown
below.

```java
Graph   g   = Graph.fromSmiles("F\\C=C/F");
Graph   h   = Functions.normaliseDirectionalLabels(g);
String  smi = h.toSmiles();
smi.equals("F/C=C\\F");
```

```
F/C=C/C              is normalised to F/C=C/C
F\C=C\C              is normalised to F/C=C/C
F/C=C\C              is normalised to F/C=C\C
F\C=C/C              is normalised to F/=C\C
C(\F)(/C)=C\C        is normalised to C(/F)(\C)=C/C
FC=C(F)C=C(F)\C=C\C  is normalised to FC=C(F)C=C(F)/C=C/C
```

## Beam me up

beam is still in development but you can obtain the latest build from the [EBI snapshots repository](http://www.ebi.ac.uk/intact/maven/nexus/content/repositories/ebi-repo-snapshots/). An example configuration for maven is shown below.

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
        <groupId>uk.ac.ebi.beam</groupId>
        <artifactId>beam-core</artifactId>
        <version>LATEST</version>
    </dependency>
    <dependency>
        <groupId>uk.ac.ebi.beam</groupId>
        <artifactId>beam-func</artifactId>
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

---------------------------------------

&trade;: SMILES is a trademark of [Daylight Chemical Information Systems](http://daylight.com/)
