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

import static org.hamcrest.CoreMatchers.hasItem;
import static org.hamcrest.CoreMatchers.hasItems;
import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

/** @author John May */
public class ChemicalGraphTest {

    @Test public void addAtoms() {
        ChemicalGraph g = new ChemicalGraph(5);
        assertThat(g.addAtom(mock(Atom.class)), is(0));
        assertThat(g.addAtom(mock(Atom.class)), is(1));
        assertThat(g.addAtom(mock(Atom.class)), is(2));
        assertThat(g.addAtom(mock(Atom.class)), is(3));
        assertThat(g.addAtom(mock(Atom.class)), is(4));
    }

    @Test public void addAtomsResize() {
        ChemicalGraph g = new ChemicalGraph(2);
        assertThat(g.addAtom(mock(Atom.class)), is(0));
        assertThat(g.addAtom(mock(Atom.class)), is(1));
        assertThat(g.addAtom(mock(Atom.class)), is(2));
        assertThat(g.addAtom(mock(Atom.class)), is(3));
        assertThat(g.addAtom(mock(Atom.class)), is(4));
    }

    @Test public void atomAccess() {
        Atom[] atoms = new Atom[]{
                mock(Atom.class),
                mock(Atom.class),
                mock(Atom.class),
                mock(Atom.class)
        };
        ChemicalGraph g = new ChemicalGraph(5);
        for (Atom a : atoms)
            g.addAtom(a);
        assertThat(g.atom(0), is(atoms[0]));
        assertThat(g.atom(1), is(atoms[1]));
        assertThat(g.atom(2), is(atoms[2]));
        assertThat(g.atom(3), is(atoms[3]));
    }

    @Test public void testOrder() {
        ChemicalGraph g = new ChemicalGraph(5);
        assertThat(g.order(), is(0));
        g.addAtom(mock(Atom.class));
        assertThat(g.order(), is(1));
        g.addAtom(mock(Atom.class));
        assertThat(g.order(), is(2));
        g.addAtom(mock(Atom.class));
        assertThat(g.order(), is(3));
        g.addAtom(mock(Atom.class));
        assertThat(g.order(), is(4));
        g.addAtom(mock(Atom.class));
        assertThat(g.order(), is(5));
    }

    @Test public void testSize() {
        ChemicalGraph g = new ChemicalGraph(5);
        g.addAtom(mock(Atom.class));
        g.addAtom(mock(Atom.class));
        g.addAtom(mock(Atom.class));

        Edge e1 = new Edge(0, 1, Bond.IMPLICIT);
        Edge e2 = new Edge(0, 1, Bond.IMPLICIT);

        assertThat(g.size(), is(0));
        g.addEdge(e1);
        assertThat(g.size(), is(1));
        g.addEdge(e2);
        assertThat(g.size(), is(2));
    }

    @Test public void testEdges() {
        ChemicalGraph g = new ChemicalGraph(5);
        g.addAtom(mock(Atom.class));
        g.addAtom(mock(Atom.class));
        g.addAtom(mock(Atom.class));
        g.addEdge(new Edge(0, 1, Bond.IMPLICIT));
        g.addEdge(new Edge(1, 2, Bond.IMPLICIT));
        assertThat(g.edges(0).size(), is(1));
        assertThat(g.edges(0), hasItem(new Edge(0, 1, Bond.IMPLICIT)));
        assertThat(g.edges(1).size(), is(2));
        assertThat(g.edges(1), hasItems(new Edge(0, 1, Bond.IMPLICIT),
                                        new Edge(1, 0, Bond.IMPLICIT)));
    }

    @Test public void testDegree() {
        ChemicalGraph g = new ChemicalGraph(5);
        g.addAtom(mock(Atom.class));
        g.addAtom(mock(Atom.class));
        g.addAtom(mock(Atom.class));
        g.addEdge(new Edge(0, 1, Bond.IMPLICIT));
        g.addEdge(new Edge(1, 2, Bond.IMPLICIT));
        assertThat(g.degree(0), is(1));
        assertThat(g.degree(1), is(2));
    }

    @Test(expected = IllegalArgumentException.class)
    public void testInvalidEdges() {
        ChemicalGraph g = new ChemicalGraph(5);
        g.edges(4);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testInvalidDegree() {
        ChemicalGraph g = new ChemicalGraph(5);
        g.degree(4);
    }

    @Test(expected = IllegalArgumentException.class)
    public void addingEdgeNonVertex1() {
        ChemicalGraph g = new ChemicalGraph(5);
        g.addAtom(mock(Atom.class));
        g.addEdge(new Edge(0, 1, Bond.IMPLICIT));
    }

    @Test(expected = IllegalArgumentException.class)
    public void addingEdgeNonVertex2() {
        ChemicalGraph g = new ChemicalGraph(5);
        g.addAtom(mock(Atom.class));
        g.addEdge(new Edge(1, 0, Bond.IMPLICIT));
    }

    @Test(expected = IllegalArgumentException.class)
    public void invalidAtom() {
        new ChemicalGraph(5).atom(-1);
    }

    @Test(expected = IllegalArgumentException.class)
    public void invalidAtom2() {
        new ChemicalGraph(5).atom(2);
    }

    @Test public void addTopology() {
        Topology t = mock(Topology.class);
        when(t.atom()).thenReturn(5);
        ChemicalGraph g = new ChemicalGraph(5);
        g.addTopology(t);
        assertThat(g.topologyOf(5), is(t));
    }

    @Test public void defaultTopology() {
        ChemicalGraph g = new ChemicalGraph(5);
        assertThat(g.topologyOf(5), is(Topology.unknown()));
    }
}
