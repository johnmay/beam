/*
 * Copyright (c) 2013, European Bioinformatics Institute (EMBL-EBI)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the FreeBSD Project.
 */

package uk.ac.ebi.grins;

import org.junit.Test;

import static org.hamcrest.CoreMatchers.hasItems;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.not;
import static org.hamcrest.CoreMatchers.sameInstance;
import static org.junit.Assert.assertThat;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;
import static uk.ac.ebi.grins.Bond.AROMATIC;
import static uk.ac.ebi.grins.Bond.DOUBLE;
import static uk.ac.ebi.grins.Bond.IMPLICIT;
import static uk.ac.ebi.grins.Bond.SINGLE;

/** @author John May */
public class ImplicitToExplicitTest {

    @Test public void cycloHexane() throws Exception {
        ChemicalGraph g = new ChemicalGraph(6);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addEdge(new Edge(0, 1, IMPLICIT));
        g.addEdge(new Edge(1, 2, IMPLICIT));
        g.addEdge(new Edge(2, 3, IMPLICIT));
        g.addEdge(new Edge(3, 4, IMPLICIT));
        g.addEdge(new Edge(4, 5, IMPLICIT));
        g.addEdge(new Edge(5, 0, IMPLICIT));

        ChemicalGraph h = new ImplicitToExplicit().transform(g);

        assertThat(g, is(not(sameInstance(h))));

        for (int u = 0; u < h.order(); u++) {
            for (Edge e : h.edges(u)) {
                assertThat(e.bond(), is(SINGLE));
            }
        }
    }

    @Test public void aromaticBenzene() throws Exception {
        ChemicalGraph g = new ChemicalGraph(6);
        g.addAtom(Atom.AromaticSubset.Carbon);
        g.addAtom(Atom.AromaticSubset.Carbon);
        g.addAtom(Atom.AromaticSubset.Carbon);
        g.addAtom(Atom.AromaticSubset.Carbon);
        g.addAtom(Atom.AromaticSubset.Carbon);
        g.addAtom(Atom.AromaticSubset.Carbon);
        g.addEdge(new Edge(0, 1, IMPLICIT));
        g.addEdge(new Edge(1, 2, IMPLICIT));
        g.addEdge(new Edge(2, 3, IMPLICIT));
        g.addEdge(new Edge(3, 4, IMPLICIT));
        g.addEdge(new Edge(4, 5, IMPLICIT));
        g.addEdge(new Edge(5, 0, IMPLICIT));

        ChemicalGraph h = new ImplicitToExplicit().transform(g);

        assertThat(g, is(not(sameInstance(h))));

        for (int u = 0; u < h.order(); u++) {
            for (Edge e : h.edges(u)) {
                assertThat(e.bond(), is(AROMATIC));
            }
        }
    }

    @Test public void kekuleBenzene() throws Exception {
        ChemicalGraph g = new ChemicalGraph(6);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addAtom(Atom.AliphaticSubset.Carbon);
        g.addEdge(new Edge(0, 1, IMPLICIT));
        g.addEdge(new Edge(1, 2, DOUBLE));
        g.addEdge(new Edge(2, 3, IMPLICIT));
        g.addEdge(new Edge(3, 4, DOUBLE));
        g.addEdge(new Edge(4, 5, IMPLICIT));
        g.addEdge(new Edge(5, 0, DOUBLE));

        ChemicalGraph h = new ImplicitToExplicit().transform(g);

        assertThat(g, is(not(sameInstance(h))));


        assertThat(h.edges(0), hasItems(new Edge(0, 1, SINGLE),
                                        new Edge(0, 5, DOUBLE)));
        assertThat(h.edges(1), hasItems(new Edge(1, 0, SINGLE),
                                        new Edge(1, 2, DOUBLE)));
        assertThat(h.edges(2), hasItems(new Edge(2, 1, DOUBLE),
                                        new Edge(2, 3, SINGLE)));
        assertThat(h.edges(3), hasItems(new Edge(3, 2, SINGLE),
                                        new Edge(3, 4, DOUBLE)));
        assertThat(h.edges(4), hasItems(new Edge(4, 3, DOUBLE),
                                        new Edge(4, 5, SINGLE)));
        assertThat(h.edges(5), hasItems(new Edge(5, 0, DOUBLE),
                                        new Edge(5, 4, SINGLE)));
    }

    @Test public void aromaticType() {
        Atom a = mock(Atom.class);
        Atom b = mock(Atom.class);
        when(a.aromatic()).thenReturn(true);
        when(b.aromatic()).thenReturn(true);
        assertThat(ImplicitToExplicit.type(a, b), is(Bond.AROMATIC));
    }

    @Test public void singleType() {
        Atom a = mock(Atom.class);
        Atom b = mock(Atom.class);

        when(a.aromatic()).thenReturn(true);
        when(b.aromatic()).thenReturn(false);
        assertThat(ImplicitToExplicit.type(a, b), is(Bond.SINGLE));

        when(a.aromatic()).thenReturn(false);
        when(b.aromatic()).thenReturn(true);
        assertThat(ImplicitToExplicit.type(a, b), is(Bond.SINGLE));

        when(a.aromatic()).thenReturn(false);
        when(b.aromatic()).thenReturn(false);
        assertThat(ImplicitToExplicit.type(a, b), is(Bond.SINGLE));
    }

    @Test public void toExplicitEdge_NonImplicitIdentity() {
        ChemicalGraph g = new ChemicalGraph(0);
        for (Bond b : Bond.values()) {
            if (b != IMPLICIT) {
                Edge e = new Edge(0, 1, SINGLE);
                assertThat(ImplicitToExplicit
                                   .toExplicitEdge(g, e), is(sameInstance(e)));
            }
        }
    }

    @Test public void toExplicitEdge() {
        ChemicalGraph g = new ChemicalGraph(2);

        Atom u = mock(Atom.class);
        Atom v = mock(Atom.class);

        when(u.aromatic()).thenReturn(false);
        when(v.aromatic()).thenReturn(false);

        g.addAtom(u);
        g.addAtom(v);

        Edge e = new Edge(0, 1, IMPLICIT);
        assertThat(ImplicitToExplicit.toExplicitEdge(g, e),
                   is(not(sameInstance(e))));
        assertThat(ImplicitToExplicit.toExplicitEdge(g, e),
                   is(new Edge(0, 1, SINGLE)));


    }
}
