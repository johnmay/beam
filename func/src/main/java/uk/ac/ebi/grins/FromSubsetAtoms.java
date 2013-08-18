package uk.ac.ebi.grins;

/**
 * Given a chemical graph with 0 or more atoms. Convert that graph to one where
 * all atoms are fully specified bracket atoms.
 *
 * @author John May
 */
final class FromSubsetAtoms
        extends AbstractFunction<ChemicalGraph, ChemicalGraph> {

    public ChemicalGraph apply(ChemicalGraph g) {

        ChemicalGraph h = new ChemicalGraph(g.order());

        for (int u = 0; u < g.order(); u++) {

            int nElectrons = 0;
            for (Edge e : g.edges(u))
                nElectrons += e.bond().electrons();

            h.addAtom(fromSubset(g.atom(u),
                                 nElectrons));
            h.addTopology(g.topologyOf(u));
        }

        // edges are unchanged
        for (Edge e : g.edges())
            h.addEdge(e);

        return h;
    }

    static Atom fromSubset(Atom a, int bondElectronSum) {

        // atom is already a non-subset atom
        if (!a.subset())
            return a;

        Element  e = a.element();
        int hCount = a.aromatic() ? e.delocalisedElectrons(bondElectronSum) / 2
                                  : e.electrons(bondElectronSum)            / 2;

        return new AtomImpl.BracketAtom(-1,
                                        a.element(),
                                        hCount,
                                        0,
                                        0,
                                        a.aromatic());
    }
}
