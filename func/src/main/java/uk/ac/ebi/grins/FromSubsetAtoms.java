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
                                 nElectrons,
                                 g.degree(u)));
            h.addTopology(g.topologyOf(u));
        }

        // edges are unchanged
        for (Edge e : g.edges())
            h.addEdge(e);

        return h;
    }

    static Atom fromSubset(Atom a, int bondElectronSum, int deg) {

        // atom is already a non-subset atom
        if (!a.subset())
            return a;

        Element e = a.element();
        int electrons = a.aromatic() ? e.delocalisedElectrons(bondElectronSum)
                                     : e.electrons(bondElectronSum);

        // XXX: if there was an odd number of electrons there was an odd number
        // or aromatic bonds (usually 1 or 3) - if there was one it was
        // only a single bond it's likely a spouting from a ring - otherwise
        // someones making our life difficult (e.g. c1=cc=cc=c1) in which we
        // 'give' back 2 free electrons for use indeterminacy the hCount
        int hCount = (electrons & 0x1) == 1 ? deg > 1 ? (electrons + 2) / 2
                                                      : electrons / 2
                                            : electrons / 2;


        return new AtomImpl.BracketAtom(-1,
                                        a.element(),
                                        hCount,
                                        0,
                                        0,
                                        a.aromatic());
    }
}
