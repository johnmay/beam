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
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Defines a labelled graph with atoms as vertex labels and bonds as edge
 * labels. Topological information around atoms can also be stored.
 *
 * @author John May
 */
public final class ChemicalGraph {

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
            for (int i = order; i < edges.length; i++)
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
    public Atom atom(int i) {
        return atoms[checkRange(i)];
    }

    /**
     * Add an labelled edge to the graph.
     *
     * @param e new edge
     * @throws IllegalArgumentException attempt to create an edge between atoms
     *                                  which do not exist
     */
    void addEdge(Edge e) {
        int u = checkRange(e.either()), v = checkRange(e.other(u));
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
        return edges[checkRange(u)].size();
    }

    /**
     * Access the edges of which vertex 'u' is an endpoint.
     *
     * @param u a vertex
     * @return edges incident to 'u'
     * @throws IllegalArgumentException attempting to access the edges of an
     *                                  atom which does not exist
     */
    public List<Edge> edges(int u) {
        return Collections.unmodifiableList(edges[checkRange(u)]);
    }

    /**
     * Determine if the vertices 'u' and 'v' are adjacent and there is an edge
     * which connects them.
     *
     * @param u a vertex
     * @param v another vertex
     * @return whether they are adjacent
     */
    public boolean adjacent(int u, int v) {
        checkRange(v);
        for (Edge e : edges(u))
            if (e.other(u) == v)
                return true;
        return false;
    }

    /**
     * Access the edge connecting two adjacent vertices.
     *
     * @param u a vertex
     * @param v another vertex (adjacent to u)
     * @return the edge connected u and v
     * @throws IllegalArgumentException u and v are not adjacent
     */
    public Edge edge(int u, int v) {
        for (Edge e : edges(u))
            if (e.other(u) == v)
                return e;
        throw new IllegalArgumentException(u + ", " + v + " are not adjacent");
    }

    /**
     * Add a topology description to the graph. The topology describes the
     * configuration around a given atom.
     *
     * @param t topology
     * @return whether the topology replaced an existing configuration
     */
    boolean addTopology(Topology t) {
        if (t == Topology.unknown())
            return false;
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
    public int order() {
        return order;
    }

    /**
     * The size is the number edges in the graph, |E|.
     *
     * @return number of edges
     */
    public int size() {
        return size;
    }

    /**
     * Permute the vertices of a graph using a given permutation.
     *
     * <blockquote><pre>
     * g = CNCO
     * h = g.permuate(new int[]{1, 0, 3, 2});
     * h = NCOC
     * </pre></blockquote>
     *
     * @param p a permutation mapping indicate the new index of each atom
     * @return a new chemical graph with the vertices permuted by the given
     *         ordering
     */
    ChemicalGraph permute(int[] p) {

        if (p.length != order)
            throw new IllegalArgumentException("permuation size should equal |V| (order)");

        ChemicalGraph g = new ChemicalGraph(order);
        g.order = order;
        g.size = size;

        for (int u = 0; u < order; u++) {
            g.atoms[p[u]] = atoms[u];
            g.addTopology(topologyOf(u).transform(p));
        }

        for (int u = 0; u < order; u++) {
            for (Edge e : edges[u]) {
                if (u < e.other(u)) {
                    int v = p[u], w = p[e.other(u)];
                    Edge f = new Edge(v, w, e.bond(u));
                    g.edges[v].add(f);
                    g.edges[w].add(f);
                }
            }
        }

        // ensure edges are in sorted order
        return g.sort();
    }

    /**
     * Access the atoms of the chemical graph.
     *
     * <blockquote><pre>
     * for (Atom a : g.atoms()) {
     *
     * }
     * </pre></blockquote>
     *
     * @return iterable of atoms
     */
    public Iterable<Atom> atoms() {
        return Arrays.asList(atoms).subList(0, order);
    }

    /**
     * Access the edges of the chemical graph.
     *
     * @return iterable of edges
     */
    public Iterable<Edge> edges() {
        List<Edge> es = new ArrayList<Edge>(size);
        for (int u = 0; u < order; u++) {
            for (Edge e : this.edges[u]) {
                if (e.other(u) > u)
                    es.add(e);
            }
        }
        return Collections.unmodifiableCollection(es);
    }

    /**
     * Apply a function to the chemical graph.
     *
     * @param f   a function which transforms a graph into something.
     * @param <T> output type of the function
     * @return the output of the function
     */
    public <T> T apply(Function<ChemicalGraph, T> f) {
        return f.apply(this);
    }

    /**
     * Sort the edges of each vertex in the chemical graph. Ensures that when
     * invoking {@link #edges(int)} the connected vertices will be in natural
     * order. The actual order of atoms does not change. The atom order can be
     * rearranged using {@link #permute(int[])}.
     *
     * <blockquote><pre>
     * g.edges(5) = {5-1, 5-0, 6-5, 5-2}  // unsorted
     * g.edges(5) = {5-0, 5-1, 5-2, 6-5}  // sorted
     * </pre></blockquote>
     *
     * @return self-reference for fluent invocation
     * @see #permute(int[])
     */
    ChemicalGraph sort() {
        for (int u = 0; u < order; u++) {
            Collections.sort(edges[u], EdgeComparator.forVertex(u));
        }
        return this;
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

    private int checkRange(int u) {
        if (u < 0 || u >= order)
            throw new IllegalArgumentException("invalid atom index: " + u);
        return u;
    }

    private static final class EdgeComparator implements Comparator<Edge> {
        private int u;

        private EdgeComparator(int u) {
            this.u = u;
        }

        static EdgeComparator forVertex(int u) {
            return new EdgeComparator(u);
        }

        @Override public int compare(Edge e, Edge f) {
            int v = e.other(u), w = f.other(u);
            if (v > w)
                return +1;
            else if (v < w)
                return -1;
            return 0;
        }
    }
}
