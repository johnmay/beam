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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Defines a labelled graph with atoms as vertex labels and bonds as edge
 * labels. Topological information around atoms can also be stored.
 *
 * @author John May
 */
final class ChemicalGraph {

    /** The vertex labels, atoms. */
    private Atom[] atoms;

    /** Incidence list storage of edges with attached bond labels. * */
    private List<Edge>[] edges;

    /** Topologies indexed by the atom which they describe. */
    private Map<Integer, Topology> topologies;

    /** Vertex and edge counts. */
    private int order, size;

    /**
     * Create a new chemical graph with expected size.
     *
     * @param expSize expected size
     */
    @SuppressWarnings("unchecked") ChemicalGraph(int expSize) {
        this.order = 0;
        this.size = 0;
        this.edges = new List[expSize];
        for (int i = 0; i < expSize; i++)
            edges[i] = new ArrayList<Edge>(4);
        this.atoms = new Atom[expSize];
        this.topologies = new HashMap<Integer, Topology>(10);
    }

    /**
     * Add an atom to the graph and return the index to which the atom was
     * added.
     *
     * @param a add an atom
     * @return index of the atom in the graph (vertex)
     */
    int addAtom(Atom a) {
        int index = order;
        if (order == atoms.length) {
            atoms = Arrays.copyOf(atoms, order * 2);
            edges = Arrays.copyOf(edges, order * 2);
            for(int i = order; i < edges.length; i++)
                edges[i] = new ArrayList<Edge>(4);
        }
        atoms[order++] = a;
        return index;
    }

    /**
     * Access the atom at the specified index.
     *
     * @param i index of the atom to access
     * @return the atom at that index
     * @throws IllegalArgumentException no atom exists
     */
    Atom atom(int i) {
        if (i < 0 || i >= order)
            throw new IllegalArgumentException("no atom at index " + i);
        return atoms[i];
    }

    /**
     * Add an labelled edge to the graph.
     *
     * @param e new edge
     * @throws IllegalArgumentException attempt to create an edge between atoms
     *                                  which do not exist
     */
    void addEdge(Edge e) {
        int u = e.either(), v = e.other(u);
        if (u < 0 || u >= order)
            throw new IllegalArgumentException("cannot add edge, vertex "
                                                       + u + " does not exist");
        if (v < 0 || v >= order)
            throw new IllegalArgumentException("cannot add edge, vertex "
                                                       + v + " does not exist");
        edges[u].add(e);
        edges[v].add(e);
        size++;
    }

    /**
     * Access the degree of vertex 'u'.
     *
     * @param u a vertex
     * @return the degree of the specified vertex
     * @throws IllegalArgumentException attempting to access the degree of an
     *                                  atom which does not exist
     */
    int degree(int u) {
        if (u < 0 || u >= order)
            throw new IllegalArgumentException("no atom at index " + u);
        return edges[u].size();
    }

    /**
     * Access the edges of which vertex 'u' is an endpoint.
     *
     * @param u a vertex
     * @return edges incident to 'u'
     * @throws IllegalArgumentException attempting to access the edges of an
     *                                  atom which does not exist
     */
    List<Edge> edges(int u) {
        if (u < 0 || u >= order)
            throw new IllegalArgumentException("no atom at index " + u);
        return Collections.unmodifiableList(edges[u]);
    }

    /**
     * Add a topology description to the graph. The topology describes the
     * configuration around a given atom.
     *
     * @param t topology
     * @return whether the topology replaced an existing configuration
     */
    boolean addTopology(Topology t) {
        return topologies.put(t.atom(), t) != null;
    }

    /**
     * Access the topology of the vertex 'u'. If no topology is defined then
     * {@link Topology#unknown()} is returned.
     *
     * @param u a vertex to access the topology of
     * @return the topology of vertex 'u'
     */
    Topology topologyOf(int u) {
        Topology t = topologies.get(u);
        return t != null ? t : Topology.unknown();
    }

    /**
     * The order is the number vertices in the graph, |V|.
     *
     * @return number of vertices
     */
    int order() {
        return order;
    }

    /**
     * The size is the number edges in the graph, |E|.
     *
     * @return number of edges
     */
    int size() {
        return size;
    }

    void clear() {
        topologies.clear();
        for (int i = 0; i < order; i++) {
            atoms[i] = null;
            edges[i].clear();
        }
        order = 0;
        size = 0;
    }
}
