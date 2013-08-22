package uk.ac.ebi.grins;

import java.util.List;

import static uk.ac.ebi.grins.Configuration.Type.None;

/**
 * Given a chemical graph with 0 or more atoms. Convert that graph to one where
 * fully specified bracket atoms which can be specified as organic subsets.
 *
 * @author John May
 */
final class ToSubsetAtoms extends AbstractFunction<Graph,Graph> {

    public Graph apply(Graph g) {

        Graph h = new Graph(g.order());

        for (int u = 0; u < g.order(); u++) {

            // only attempt subset conversion if no known topology
            Topology t = g.topologyOf(u);

            if (t.type() == None) {
                h.addAtom(toSubset(g.atom(u),
                                   electronSum(g.edges(u), g)));
            } else {
                h.addAtom(g.atom(u));
                h.addTopology(t);
            }
        }

        // edges are unchanged
        for (Edge e : g.edges())
            h.addEdge(e);

        return h;
    }

    static Atom toSubset(Atom a, int nElectrons) {

        // atom is already a subset atom
        if (a.subset())
            return a;

        // element is not organic and thus cannot be part of the subset
        if (!a.element().organic())
            return a;


        // if any of these values are set the atom cannot be a subset atom
        if (a.charge() != 0 || a.atomClass() != 0 || a.isotope() >= 0)
           return a;

        // does the implied availableElectrons from the bond order sum match that
        // which was stored - if aromatic we only check the lowest valence state
        int impliedHCount = a.aromatic() ? a.element().availableDelocalisedElectrons(nElectrons) / 2
                                         : a.element().availableElectrons(nElectrons) / 2;

        // mismatch in number of hydrogens we must write this as a bracket atom
        if (impliedHCount != a.hydrogens())
            return a;

        if (a.aromatic()) {
            return AtomImpl.AromaticSubset.ofElement(a.element());
        } else {
            return AtomImpl.AliphaticSubset.ofElement(a.element());
        }
    }

    private int electronSum(final List<Edge> es,
                            final Graph g) {
        int nElectrons = 0;
        for (Edge e : es) {
            int u = e.either();
            int v = e.other(u);
            nElectrons += e.bond().electrons(g.atom(u), g.atom(v));
        }
        return nElectrons;
    }
}
