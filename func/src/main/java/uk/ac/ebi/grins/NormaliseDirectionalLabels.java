package uk.ac.ebi.grins;

import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Normalise directional labels such that the first label is always a '/'. Given
 * a molecule with directional bonds {@code F\C=C\F} the labels are normalised
 * to be {@code F/C=C/F}.
 *
 * @author John May
 */
public class NormaliseDirectionalLabels
        extends AbstractFunction<ChemicalGraph, ChemicalGraph> {

    @Override public ChemicalGraph apply(ChemicalGraph g) {
        Traversal traversal = new Traversal(g);
        ChemicalGraph h = new ChemicalGraph(g.order());

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

        return h;
    }

    private static final class Traversal {

        private final ChemicalGraph g;
        private final boolean[]     visited;
        private final int[]         ordering;
        private       int           i;
        private Map<Edge, Edge> acc = new HashMap<Edge, Edge>();

        private Set<Edge>    doubleBonds = new HashSet<Edge>();
        private Set<Integer> adj         = new HashSet<Integer>();

        private Traversal(ChemicalGraph g) {
            this.g = g;
            this.visited = new boolean[g.order()];
            this.ordering = new int[g.order()];

            for (int u = 0; u < g.order(); u++) {
                if (!visited[u])
                    visit(u, u);
            }
            for (Edge e : doubleBonds) {
                flip(g, e, acc);
            }
        }

        private void visit(int p, int u) {
            visited[u] = true;
            ordering[u] = i++;
            for (Edge e : g.edges(u)) {
                int v = e.other(u);

                if (!visited[v]) {
                    if (e.bond() == Bond.DOUBLE && hasAdjDirectionalLabels(g, e)) {

                        // only the first bond we encounter in an isolated system
                        // is marked - if we need to flip the other we propagate
                        // this down the chain
                        if (!adj.contains(u) && !adj.contains(v))
                            doubleBonds.add(e);

                        // to stop adding other we mark all vertices adjacent to the
                        // double bond
                        for (Edge f : g.edges(u))
                            adj.add(f.other(u));
                        for (Edge f : g.edges(v))
                            adj.add(f.other(v));
                    }
                    visit(u, v);
                }
            }
        }

        private boolean hasAdjDirectionalLabels(ChemicalGraph g, Edge e) {
            int u = e.either();
            int v = e.other(u);

            for (Edge f : g.edges(u))
                if (f.bond() == Bond.UP || f.bond() == Bond.DOWN)
                    return true;
            for (Edge f : g.edges(v))
                if (f.bond() == Bond.UP || f.bond() == Bond.DOWN)
                    return true;

            return false;
        }

        private void flip(ChemicalGraph g, Edge e, Map<Edge, Edge> acc) {

            int u = e.either();
            int v = e.other(u);

            if (ordering[u] < ordering[v]) {
                Edge first = firstDirectionalLabel(g, u);
                if (first != null) {
                    flip(first, u);
                } else {
                    first = firstDirectionalLabel(g, v);
                    flip(first, v);
                }
            } else {
                Edge first = firstDirectionalLabel(g, v);
                if (first != null) {
                    flip(first, v);
                } else {
                    first = firstDirectionalLabel(g, u);
                    flip(first, u);
                }
            }

        }

        private void flip(Edge first, int u) {
            if (ordering[first.other(u)] < ordering[u]) {
                if (first.bond(u) == Bond.UP)
                    invertExistingDirectionalLabels(g,
                                                    new BitSet(),
                                                    acc,
                                                    u);
            } else {
                if (first.bond(u) == Bond.DOWN)
                    invertExistingDirectionalLabels(g,
                                                    new BitSet(),
                                                    acc,
                                                    u);
            }
        }

        Edge firstDirectionalLabel(ChemicalGraph g, int u) {
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

        private void invertExistingDirectionalLabels(ChemicalGraph g,
                                                     BitSet visited,
                                                     Map<Edge, Edge> replacement,
                                                     int u) {
            visited.set(u);
            if (g.topologyOf(u) == null)
                return;
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
                    invertExistingDirectionalLabels(g, visited, replacement, v);
                }
            }
        }
    }
}
