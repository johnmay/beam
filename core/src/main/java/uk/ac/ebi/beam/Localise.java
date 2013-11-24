package uk.ac.ebi.beam;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

/**
 * Utility to localise aromatic bonds.
 *
 * @author John May
 */
final class Localise {

    private final Graph delocalised, localised;
    private final BitSet subset;
    private final Map<Edge, Edge> edgeAssignments = new HashMap<Edge, Edge>();

    private Localise(Graph delocalised, BitSet subset, BitSet aromatic) throws InvalidSmilesException {

        this.delocalised = delocalised;
        this.localised = new Graph(delocalised.order());
        this.subset = subset;

        // make initial matching - then improve it
        Matching m = MaximumMatching.maximise(delocalised,
                                              ArbitraryMatching.of(delocalised, subset),
                                              IntSet.fromBitSet(subset));


        // the edge assignments have all aromatic bonds going to single bonds
        // we now use the matching to update these 
        for (int u = subset.nextSetBit(0); u >= 0; u = subset.nextSetBit(u + 1)) {
            if (m.unmatched(u)) {
                throw new InvalidSmilesException("A valid kekule structure could not be assigned to: " +
                                                         delocalised.toSmiles());
            }
            int v = m.other(u);
            subset.clear(v);
            edgeAssignments.put(delocalised.edge(u, v),
                                Bond.DOUBLE_AROMATIC.edge(u, v));
        }

        // create the new graph
        for (int v = 0; v < delocalised.order(); v++) {
            localised.addAtom(delocalised.atom(v).toAliphatic());
            localised.addTopology(delocalised.topologyOf(v));
        }

        for (Edge orgEdge : delocalised.edges()) {
            Edge newEdge = edgeAssignments.get(orgEdge);
            if (newEdge != null) {
                localised.addEdge(newEdge);
            }
            else {
                switch (orgEdge.bond()) {
                    case SINGLE:
                        localised.addEdge(Bond.IMPLICIT.edge(orgEdge.either(),
                                                             orgEdge.other(orgEdge.either())));
                        break;
                    case AROMATIC:
                        localised.addEdge(Bond.IMPLICIT_AROMATIC.edge(orgEdge.either(),
                                                                      orgEdge.other(orgEdge.either())));
                        break;
                    case IMPLICIT:
                        int u = orgEdge.either();
                        int v = orgEdge.other(u);
                        if (aromatic.get(u) && aromatic.get(v))
                            localised.addEdge(Bond.IMPLICIT_AROMATIC.edge(orgEdge.either(),
                                                                          orgEdge.other(orgEdge.either())));
                        else
                            localised.addEdge(orgEdge);
                        break;
                    default:
                        localised.addEdge(orgEdge);
                        break;
                }
            }
        }
    }


    static BitSet buildSet(Graph g, BitSet aromatic) {

        BitSet undecided = new BitSet(g.order());

        for (int v = 0; v < g.order(); v++) {
            if (g.atom(v).aromatic()) {
                aromatic.set(v);
                if (!predetermined(g, v))
                    undecided.set(v);
            }
        }

        return undecided;
    }

    static boolean predetermined(Graph g, int v) {

        Atom a = g.atom(v);

        int q = a.charge();
        int deg = g.degree(v) + g.implHCount(v);

        for (Edge e : g.edges(v)) {
            if (e.bond() == Bond.DOUBLE) {
                if (q == 0 && (a.element() == Element.Nitrogen || (a.element() == Element.Sulfur && deg > 3))
                        && g.atom(e.other(v)).element() == Element.Oxygen)
                    return false;
                return true;
            }
            // triple or quadruple bond - we don't need to assign anymore p electrons
            else if (e.bond().order() > 2) {
                return true;
            }
        }

        // no pi bonds does the degree and charge indicate that
        // there can be no other pi bonds
        switch (a.element()) {
            case Carbon:
                return (q == 1 || q == -1) && deg == 3;
            case Silicon:
            case Germanium:
                return q < 0 ? true : false;
            case Nitrogen:
            case Phosphorus:
            case Arsenic:
            case Antimony:
                if (q == 0)
                    return deg == 3 ? true : false;
                else if (q == 1)
                    return deg == 3 ? false : true;
                else
                    return true;
            case Oxygen:
            case Sulfur:
            case Selenium:
            case Tellurium:
                if (q == 0)
                    return deg == 2 ? true : false;
                else if (q < 0)
                    return true;
                else
                    return false;
        }

        return false;
    }

    static Graph localise(Graph delocalised) throws InvalidSmilesException {

        // nothing to do
        if (!delocalised.isDelocalised())
            return delocalised;

        BitSet aromatic = new BitSet();
        BitSet subset   = buildSet(delocalised, aromatic);
        if (hasOddCardinality(subset))
            throw new InvalidSmilesException("A valid kekule structure could not be assigned to: " +
                                                     delocalised.toSmiles());
        return new Localise(delocalised, subset, aromatic).localised;
    }

    private static boolean hasOddCardinality(BitSet s) {
        return (s.cardinality() & 0x1) == 1;
    }
}
