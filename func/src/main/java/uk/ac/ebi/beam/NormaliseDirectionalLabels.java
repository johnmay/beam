package uk.ac.ebi.beam;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Normalise directional labels such that the first label is always a '/'. Given
 * a molecule with directional bonds {@code F\C=C\F} the labels are normalised
 * to be {@code F/C=C/F}.
 *
 * @author John May
 */
final class NormaliseDirectionalLabels
        extends AbstractFunction<Graph, Graph> {

    @Override public Graph apply(Graph g) {
        Traversal traversal = new Traversal(g);
        Graph h = new Graph(g.order());

        // copy atom/topology information this is unchanged
        for (int u = 0; u < g.order(); u++) {
            h.addAtom(g.atom(u));
            h.addTopology(g.topologyOf(u));
        }


        // change edges (only changed added to replacement)
        for (int u = 0; u < g.order(); u++) {
            for (final Edge e : g.edges(u)) {
                if (e.other(u) > u) {
                    if (traversal.acc.containsKey(e)) {
                        h.addEdge(traversal.acc.get(e));
                    } else {
                        h.addEdge(e);
                    }
                }
            }
        }

        return h.sort(new Graph.CanOrderFirst());
    }

    private static final class Traversal {

        private final Graph     g;
        private final boolean[] visited;
        private final int[]     ordering;
        private       int       i;
        private Map<Edge, Edge> acc = new HashMap<Edge, Edge>();

        private List<Edge>   doubleBonds = new ArrayList<Edge>();
        private Set<Integer> adj         = new HashSet<Integer>();

        private Traversal(Graph g) {
            this.g = g;
            this.visited = new boolean[g.order()];
            this.ordering = new int[g.order()];

            BitSet dbAtoms = new BitSet();
            for (int u = 0; u < g.order(); u++) {
                if (!visited[u])
                    dbAtoms.or(visit(u, u));
            }

            for (Edge e : doubleBonds) {
                flip(g, e, dbAtoms);
            }
        }

        private BitSet visit(int p, int u) {
            visited[u] = true;
            ordering[u] = i++;
            BitSet dbAtoms = new BitSet();
            for (Edge e : g.edges(u)) {
                int v = e.other(u);                             
              
                if (!visited[v]) {
                    if (e.bond() == Bond.DOUBLE && hasAdjDirectionalLabels(g, e)) {

                        dbAtoms.set(u);
                        dbAtoms.set(v);
                        
                        // only the first bond we encounter in an isolated system
                        // is marked - if we need to flip the other we propagate
                        // this down the chain
                        boolean newSystem = !adj.contains(u) && !adj.contains(v);                               

                        // to stop adding other we mark all vertices adjacent to the
                        // double bond
                        for (Edge f : g.edges(u))
                            adj.add(f.other(u));
                        for (Edge f : g.edges(v))
                            adj.add(f.other(v));

                        if (newSystem)
                            doubleBonds.add(e);
                    }
                    dbAtoms.or(visit(u, v));
                }
            }
            
            return dbAtoms;
        }

        private boolean hasAdjDirectionalLabels(Graph g, Edge e) {
            int u = e.either();
            int v = e.other(u);
            return hasAdjDirectionalLabels(g, u) && hasAdjDirectionalLabels(g, v);
        }
        
        private boolean hasAdjDirectionalLabels(Graph g, int u) {           
            for (Edge f : g.edges(u))
                if (f.bond().directional())
                    return true;               
            return false;
        }

        private void flip(Graph g, Edge e, BitSet dbAtoms) {

            int u = e.either();
            int v = e.other(u);

            if (ordering[u] < ordering[v]) {
                Edge first = firstDirectionalLabel(g, u);
                if (first != null) {
                    flip(first, u, dbAtoms);
                } else {
                    first = firstDirectionalLabel(g, v);
                    flip(first, v, dbAtoms);
                }
            } else {
                Edge first = firstDirectionalLabel(g, v);
                if (first != null) {
                    flip(first, v, dbAtoms);
                } else {
                    first = firstDirectionalLabel(g, u);
                    flip(first, u, dbAtoms);
                }
            }
        }

        private void flip(Edge first, int u, BitSet dbAtoms) {
            if (ordering[first.other(u)] < ordering[u]) {
                if (first.bond(u) == Bond.UP)
                    invertExistingDirectionalLabels(g,
                                                    new BitSet(),
                                                    acc,
                                                    dbAtoms,
                                                    u);
            } else {
                if (first.bond(u) == Bond.DOWN)
                    invertExistingDirectionalLabels(g,
                                                    new BitSet(),
                                                    acc,
                                                    dbAtoms,
                                                    u);
            }
        }

        Edge firstDirectionalLabel(Graph g, int u) {
            Edge first = null;
            for (Edge f : g.edges(u)) {
                if (f.bond() == Bond.UP || f.bond() == Bond.DOWN) {
                    if (first == null || ordering[f.other(u)] < ordering[first
                            .other(u)])
                        first = f;
                }
            }
            return first;
        }

        private void invertExistingDirectionalLabels(Graph g,
                                                     BitSet visited,
                                                     Map<Edge, Edge> replacement,
                                                     BitSet dbAtoms,
                                                     int u) {
            visited.set(u);
            for (Edge e : g.edges(u)) {
                int v = e.other(u);
                if (!visited.get(v)) {
                    Edge f = replacement.get(e);
                    if (f != null) {
                        replacement.put(e,
                                        f.inverse());
                    } else {
                        replacement.put(e,
                                        e.inverse());
                    }
                    if (dbAtoms.get(v))
                        invertExistingDirectionalLabels(g, visited, replacement, dbAtoms, v);
                }
            }
        }
    }
}
