package uk.ac.ebi.grins;

import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static uk.ac.ebi.grins.Configuration.Type.DoubleBond;

/**
 * Given a chemical graph with atom-centric double bond stereo configurations
 * (trigonal topology) - remove the topology but add in direction up/down edge
 * labels.
 *
 * @author John May
 */
final class FromTrigonalTopology {

    ChemicalGraph transform(ChemicalGraph g) {
        ChemicalGraph h = new ChemicalGraph(g.order());

        // copy atom/topology information this is unchanged
        for (int u = 0; u < g.order(); u++) {
            if (g.topologyOf(u).configuration().type() == DoubleBond) {
                h.addAtom(reducedAtom(g, u));
            } else {
                h.addAtom(g.atom(u));
                h.addTopology(g.topologyOf(u));
            }
        }

        Map<Edge, Edge> replacements = new Traversal(g).replacement;


        // append the edges, replacing any which need to be changed
        for (int u = 0; u < g.order(); u++) {
            for (Edge e : g.edges(u)) {
                if (e.other(u) > u) {
                    Edge replacement = replacements.get(e);
                    if (replacement != null)
                        e = replacement;
                    h.addEdge(e);
                }
            }
        }

        return h;
    }

    private Atom reducedAtom(ChemicalGraph g, int u) {
        Atom a = g.atom(u);
        if (a.isotope() < 0 && a.atomClass() == 0 && a.charge() == 0) {
            int nElectrons = 0;
            for (Edge e : g.edges(u)) {
                nElectrons += e.bond().electrons();
            }
            if (a.hydrogens() == a.element()
                                  .implicitHydrogens(nElectrons / 2)) {
                if (a.aromatic()) {
                    // preliminary
                    switch (a.element()) {
                        case Carbon:
                            return Atom.AromaticSubset.Carbon;
                    }
                } else {
                    switch (a.element()) {
                        case Carbon:
                            return Atom.AliphaticSubset.Carbon;
                    }
                }
            }
        }
        return a;
    }

    private static final class Traversal {

        private final ChemicalGraph g;
        private final boolean[]     visited;
        private final int[]         ordering;
        private       int           i;
        private Map<Edge, Edge> replacement = new HashMap<Edge, Edge>();

        private static final Bond[] labels = new Bond[]{Bond.DOWN, Bond.UP};

        private Traversal(ChemicalGraph g) {
            this.g = g;
            this.visited = new boolean[g.order()];
            this.ordering = new int[g.order()];


            for (int u = 0; u < g.order(); u++) {
                if (!visited[u])
                    visit(u, u);
            }
        }

        private void visit(int p, int u) {

            visited[u] = true;

            // offset - the index of the edge with a double bond label
            int offset = -1;

            List<Edge> es = g.edges(u);
            for (int i = 0; i < es.size(); i++) {
                Edge e = es.get(i);
                int v = e.other(u);
                if (!visited[v])
                    visit(u, v);
                ordering[v] = 2 + i;
                if (e.bond() == Bond.DOUBLE)
                    offset = i;
            }

            ordering[p] = 0;
            ordering[u] = 1;


            Topology t = g.topologyOf(u);

            if (t.configuration().type() == DoubleBond) {

                if (offset < 0)
                    throw new IllegalArgumentException("found atom-centric double bond" +
                                                               "specifiation but no double bond label.");

                // order the topology to ensure it matches the traversal order
                Topology topology = t.orderBy(ordering);

                // labelling start depends on configuration ...
                int j = topology.configuration()
                                .shorthand() == Configuration.ANTI_CLOCKWISE ? 0
                                                                             : 1;

                // ... and which end of the double bond we're looking from
                if (ordering[es.get(offset).other(u)] < ordering[u]) {

                } else if (es.size() == 2 &&
                        ordering[u] < ordering[es.get((offset + 1) % es.size())
                                                 .other(u)]) {
                    j++;
                }

                // now create the new labels for the non-double bond atoms
                for (int i = 1; i < es.size(); i++) {
                    Edge e = es.get((offset + i) % es.size());
                    Bond label = labels[j++ % 2];

                    Edge f = new Edge(u,
                                      e.other(u),
                                      label);
                    Edge existing = replacement.get(e);

                    // check for conflict - need to rewrite existing labels
                    if (existing != null && existing.bond(u) != label) {
                        BitSet visited = new BitSet();
                        visited.set(u);
                        invertExistingDirectionalLabels(visited, e.other(u));
                    }
                    replacement.put(e, f);
                }
            }
        }

        private void invertExistingDirectionalLabels(BitSet visited, int u) {
            visited.set(u);
            if (g.topologyOf(u) == null)
                return;
            for (Edge e : g.edges(u)) {
                int v = e.other(u);
                if (!visited.get(v)) {
                    Edge f = replacement.get(e);
                    if (f != null) {
                        replacement.put(e,
                                        new Edge(u, v, f.bond(u).inverse()));
                    }
                    invertExistingDirectionalLabels(visited, v);
                }
            }
        }

    }
}
