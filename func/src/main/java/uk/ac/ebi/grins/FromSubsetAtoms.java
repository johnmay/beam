package uk.ac.ebi.grins;

import java.util.List;

/**
 * Given a chemical graph with 0 or more atoms. Convert that graph to one where
 * all atoms are fully specified bracket atoms.
 *
 * @author John May
 */
final class FromSubsetAtoms extends AbstractFunction<ChemicalGraph,ChemicalGraph> {

    public ChemicalGraph apply(ChemicalGraph g) {

        ChemicalGraph h = new ChemicalGraph(g.order());

        for (int u = 0; u < g.order(); u++) {
            h.addAtom(fromSubset(g.atom(u),
                                 bondOrderSum(g.edges(u))));
            h.addTopology(g.topologyOf(u));
        }

        // edges are unchanged
        for (Edge e : g.edges())
            h.addEdge(e);

        return h;
    }

    static Atom fromSubset(Atom a, int bondOrderSum) {

        // atom is already a non-subset atom
        if (!a.subset())
            return a;

        int hCount = a.element().implicitHydrogens(bondOrderSum);

        return new AtomImpl.BracketAtom(-1,
                                        a.element(),
                                        hCount,
                                        0,
                                        0,
                                        a.aromatic());
    }

    private int bondOrderSum(final List<Edge> es) {
        int nElectrons = 0;
        for (Edge e : es)
            nElectrons += e.bond().electrons();
        return nElectrons / 2;
    }

}
