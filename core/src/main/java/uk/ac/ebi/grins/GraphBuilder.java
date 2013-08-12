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

/**
 * Provides the ability to incrementally build up a chemical graph from atoms
 * and their connections.
 *
 * <blockquote><pre>
 * ChemicalGraph g = GraphBuilder.create(3)
 *                               .add(Carbon, 3)
 *                               .add(AtomBuilder.aliphatic(Carbon)
 *                                               .hydrogens(2)
 *                                               .build())
 *                               .add(Oxygen, 1)
 *                               .add(0, 1)
 *                               .add(1, 2)
 *                               .add(2, 3)
 *                               .build();
 * </pre></blockquote>
 *
 * @author John May
 */
public final class GraphBuilder {

    /** Current we just use the non-public methods of the actual graph object. */
    private final ChemicalGraph g;

    /**
     * Internal constructor.
     *
     * @param nAtoms expected number of atoms
     */
    private GraphBuilder(int nAtoms) {
        this.g = new ChemicalGraph(nAtoms);
    }

    public static GraphBuilder create(int n) {
        return new GraphBuilder(n);
    }

    /**
     * Add an aliphatic element with the specified number of carbons.
     *
     * @param e      element
     * @param hCount number of hydrogens
     * @return graph builder for adding more atoms/connections
     */
    public GraphBuilder add(Element e, int hCount) {
        return add(AtomBuilder.aliphatic(e)
                              .hydrogens(hCount)
                              .build());
    }

    /**
     * Add an atom to the graph.
     *
     * @param a the atom to add
     * @return graph builder for adding more atoms/connections
     */
    public GraphBuilder add(Atom a) {
        g.addAtom(a);
        return this;
    }

    /**
     * Add an edge to the graph.
     *
     * @param e the edge to add
     * @return graph builder for adding more atoms/connections
     */
    public GraphBuilder add(Edge e) {
        g.addEdge(e);
        return this;
    }

    /**
     * Connect the vertices u and v with an {@link Bond#IMPLICIT} bond label.
     *
     * @param u a vertex
     * @param v another vertex
     * @return graph builder for adding more atoms/connections
     */
    public GraphBuilder add(int u, int v) {
        add(u, v, Bond.IMPLICIT);
        return this;
    }

    /**
     * Connect the vertices u and v with the specified bond label.
     *
     * @param u a vertex
     * @param v another vertex
     * @return graph builder for adding more atoms/connections
     */
    public GraphBuilder add(int u, int v, Bond b) {
        add(b.edge(u, v));
        return this;
    }

    /**
     * Connect the vertices u and v with a single bond.
     *
     * @param u a vertex
     * @param v another vertex
     * @return graph builder for adding more atoms/connections
     */
    public GraphBuilder singleBond(int u, int v) {
        return add(u, v, Bond.SINGLE);
    }

    /**
     * Connect the vertices u and v with an aromatic bond.
     *
     * @param u a vertex
     * @param v another vertex
     * @return graph builder for adding more atoms/connections
     */
    public GraphBuilder aromaticBond(int u, int v) {
        return add(u, v, Bond.AROMATIC);
    }

    /**
     * Connect the vertices u and v with a double bond.
     *
     * @param u a vertex
     * @param v another vertex
     * @return graph builder for adding more atoms/connections
     */
    public GraphBuilder doubleBond(int u, int v) {
        return add(u, v, Bond.DOUBLE);
    }

    /**
     * Finalise and build the chemical graph.
     *
     * @return chemical graph instance
     */
    public ChemicalGraph build() {
        return g;
    }
}
