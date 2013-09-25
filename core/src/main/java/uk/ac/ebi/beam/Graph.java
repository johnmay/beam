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

package uk.ac.ebi.beam;

import java.io.IOException;
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
public final class Graph {

    /** The vertex labels, atoms. */
    private Atom[] atoms;

    /** Incidence list storage of edges with attached bond labels. * */
    private List<Edge>[] edges;

    /** Topologies indexed by the atom which they describe. */
    private Map<Integer, Topology> topologies;

    /** Vertex and edge counts. */
    private int order, size;
    
    /** Indicates at least part of the molecule is delocalised. */
    private boolean delocalised;

    /**
     * Create a new chemical graph with expected size.
     *
     * @param expSize expected size
     */
    @SuppressWarnings("unchecked") Graph(int expSize) {
        this.order = 0;
        this.size = 0;
        this.edges = new List[expSize];
        for (int i = 0; i < expSize; i++)
            edges[i] = new ArrayList<Edge>(4);
        this.atoms = new Atom[expSize];
        this.topologies = new HashMap<Integer, Topology>(10);
    }

    /**
     * (internal) - set the atom label at position 'i'.
     *
     * @param i index
     * @param a atom
     */
    void setAtom(int i, Atom a) {
        atoms[i] = a;
    }

    /**
     * Resize the graph if we are at maximum capacity.
     */
    private void ensureCapacity() {
        if (order >= atoms.length) {
            atoms = Arrays.copyOf(atoms, order * 2);
            edges = Arrays.copyOf(edges, order * 2);
            for (int i = order; i < edges.length; i++)
                edges[i] = new ArrayList<Edge>(4);
        }    
    }    
    
    /**
     * Add an atom to the graph and return the index to which the atom was
     * added.
     *
     * @param a add an atom
     * @return index of the atom in the graph (vertex)
     */
    int addAtom(Atom a) {
        ensureCapacity();
        atoms[order++] = a;
        delocalised = delocalised | a.aromatic();
        return order - 1;
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
     * Access the vertices adjacent to 'u' in <b>sorted</b> order. This
     * convenience method is provided to assist in configuring atom-based stereo
     * using the {@link #configurationOf(int)} method. For general purpose
     * access to the neighbors of a vertex the {@link #edges(int)} is
     * preferred.
     *
     * @param u a vertex
     * @return fixed-size array of vertices
     * @see #configurationOf(int)
     */
    public int[] neighbors(int u) {
        List<Edge> es = edges[checkRange(u)];
        int[] vs = new int[es.size()];
        int deg = es.size();
        for (int i = 0; i < deg; i++)
            vs[i] = es.get(i).other(u);
        Arrays.sort(vs);
        return vs;
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
     * The number of implied (or labelled) hydrogens for the vertex 'u'. Note
     * the count does not include any bonded vertices which may also be
     * hydrogen.
     *
     * @param u the vertex to access the implicit h count for.
     * @return the number of implicit hydrogens
     */
    public int implHCount(int u) {
        return atom(checkRange(u)).hydrogens(this, u);
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
     * Replace an edge in the graph.
     *
     * @param org the original edge
     * @param rep the replacement
     */
    void replace(Edge org, Edge rep) {

        int u = org.either();
        int v = org.other(u);

        for (int i = 0; i < edges[u].size(); i++) {
            if (edges[u].get(i) == org) {
                edges[u].set(i, rep);
            }
        }

        for (int i = 0; i < edges[v].size(); i++) {
            if (edges[v].get(i) == org) {
                edges[v].set(i, rep);
            }
        }

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
     * Provides the stereo-configuration of the atom label at vertex 'u'. The
     * configuration describes the relative-stereo as though the atoms were
     * arranged by atom number. <br/><br/>
     *
     * <b>Further Explanation for Tetrahedral Centres</b> As an example the
     * molecule {@code O[C@]12CCCC[C@@]1(O)CCCC2} has two tetrahedral centres.
     * <br/> 1. The first one is on vertex '1' and looking from vertex '0' the
     * other neighbors [6, 11, 2] proceed anti-clockwise ('@') - note ring
     * bonds. It is easy to see that if we use the natural order of the molecule
     * and order the neighbor [2, 6, 11] the winding is still anti-clockwise and
     * '@TH1' is returned. 2. The second centre is on vertex '6' and looking
     * from vertex '5' the ordering proceeds as [1, 7, 8] with clockwise
     * winding. When we arrange the atoms by their natural order we will now be
     * looking from vertex '1' as it is the lowest. The other neighbors then
     * proceed in the order [5, 7, 8]. Drawing out the configuration it's clear
     * that we look from vertex '1' instead of '5' the winding is now
     * anti-clockwise and the configuration is also '@TH1'.
     *
     * @param u a vertex in the graph
     * @return The configuration around
     */
    public Configuration configurationOf(int u) {

        Topology t = topologyOf(u);

        if (t == Topology.unknown())
            return t.configuration();

        // identity permutation
        int[] p = new int[order];
        for (int i = 0; i < order; i++)
            p[i] = i;

        return t.orderBy(p).configuration();
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
     * Convenience method to create a graph from a provided SMILES string.
     *
     * @param smi string containing SMILES line notation.
     * @return graph instance from the SMILES
     * @throws InvalidSmilesException thrown if there was a syntax error while
     *                                parsing the SMILES.
     */
    public static Graph fromSmiles(String smi) throws
                                               IOException {
        if (smi == null)
            throw new NullPointerException("no SMILES provided");
        return Parser.parse(smi);
    }

    /**
     * Convenience method to write a SMILES string for the current configuration
     * of the molecule.
     *
     * @return the SMILES string for the molecule.
     */
    public String toSmiles() {
        return Generator.generate(this);
    }

    /**
     * Delocalise a kekulé graph representation to one with <i>aromatic</i>
     * bonds. The original graph remains unchanged.
     * 
     * TODO: more explanation
     *
     * @return aromatic representation
     */
    public Graph aromatic() {
        // note Daylight use SSSR - should update and use that by default but
        // provide the AllCycles method
        return AllCycles.daylightModel(this).aromaticForm();
    }
  
    /**
     * Localise delocalized (aromatic) bonds in this molecule producing the
     * Kekulé form. The original graph is not modified.
     *
     * <blockquote><pre>
     * Graph furan        = Graph.fromSmiles("o1cccc1");
     * Graph furan_kekule = furan.kekule();
     * </pre></blockquote>
     *
     * If the graph could not be converted to a kekulé representation then a
     * checked exception is thrown. Graphs cannot be converted if their
     * structures are erroneous and there is no valid way to assign the
     * delocalised electrons. <p/>
     *
     * Some reasons are shown below.
     *
     * <blockquote><pre>
     * n1cncc1             pyrole (incorrect) could be either C1C=NC=N1 or
     * N1C=CN=C1
     * n1c[nH]cc1          pyrole (correct)
     *
     * [Hg+2][c-]1ccccc1   mercury(2+) ion benzenide (incorrect)
     * [Hg+2].[c-]1ccccc1  mercury(2+) ion benzenide (correct)
     * </pre></blockquote>
     *
     * @return kekulé representation
     * @throws InvalidSmilesException molecule exploded on contact with reality
     */
    public Graph kekule() throws InvalidSmilesException {
        return Localise.localise(this);
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
    Graph permute(int[] p) {

        if (p.length != order)
            throw new IllegalArgumentException("permuation size should equal |V| (order)");

        Graph g = new Graph(order);
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
                if (e.other(u) < u)
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
    <T> T apply(Function<Graph, T> f) throws Exception {
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
    Graph sort() {
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
    
    void markDelocalised() {
        delocalised = true;
    }
    
    boolean isDelocalised() {
        return delocalised;
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
