package uk.ac.ebi.grins;

import java.util.Arrays;
import java.util.List;

import static uk.ac.ebi.grins.Configuration.UNKNOWN;

/**
 * Given a chemical graph with 0 or more atoms. Convert that graph to one where
 * fully specified bracket atoms which can be specified as organic subsets.
 *
 * @author John May
 */
final class ToSubsetAtoms {

    ChemicalGraph transform(ChemicalGraph g) {

        ChemicalGraph h = new ChemicalGraph(g.order());

        for (int u = 0; u < g.order(); u++) {

            // only attempt subset conversion if no known topology
            Topology t = g.topologyOf(u);

            if (t.configuration() == UNKNOWN) {
                h.addAtom(toSubset(g.atom(u),
                                   bondOrderSum(g.edges(u))));
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

    static Atom toSubset(Atom a, int bondOrderSum) {

        // atom is already a subset atom
        if (a.subset())
            return a;

        // element is not organic and thus cannot be part of the subset
        if (!a.element().organic())
            return a;


        // if any of these values are set the atom cannot be a subset atom
        if (a.charge() != 0 || a.atomClass() != 0 || a.isotope() >= 0)
           return a;

        // does the implied hydrogens from the bond order sum match that
        // which was stored
        int impliedHCount = a.element().implicitHydrogens(bondOrderSum);
        if (impliedHCount != a.hydrogens())
            return a;

        if (a.aromatic()) {
            return AtomImpl.AromaticSubset.ofElement(a.element());
        } else {
            return AtomImpl.AliphaticSubset.ofElement(a.element());
        }
    }

    private int bondOrderSum(final List<Edge> es) {
        int nElectrons = 0;
        for (Edge e : es)
            nElectrons += e.bond().electrons();
        return nElectrons / 2;
    }
}
