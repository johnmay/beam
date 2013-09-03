package uk.ac.ebi.beam;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

/**
 * Utility to localise aromatic bonds.
 *
 * @author John May
 */
public final class Localise {

    private final Graph delocalised, localised;
    private final IntSet subset;
    private final Map<Edge, Edge> edgeAssignments = new HashMap<Edge, Edge>();

    private Localise(Graph delocalised) throws InvalidSmilesException {

        this.delocalised = delocalised;
        this.localised = new Graph(delocalised.order());
        this.subset = buildSet(delocalised);

        // make initial matching - then improve it
        Matching m = MaximumMatching.maximise(delocalised,
                                              ArbitraryMatching.of(delocalised, subset),
                                              subset);


        // the edge assignments have all aromatic bonds going to single bonds
        // we now use the matching to update these 
        for (Tuple tuple : m.matches()) {
            int u = tuple.first();
            int v = tuple.second();
            edgeAssignments.put(delocalised.edge(u, v),
                                Bond.DOUBLE.edge(u, v));
        }

        // create the new graph
        for (int v = 0; v < delocalised.order(); v++) {
            localised.addAtom(delocalised.atom(v).toAliphatic());
            localised.addTopology(delocalised.topologyOf(v));
        }
        for (Edge orgEdge : delocalised.edges()) {
            Edge newEdge = edgeAssignments.get(orgEdge);
            if (newEdge != null)
                localised.addEdge(newEdge);
            else if (orgEdge.bond() == Bond.AROMATIC)
                localised.addEdge(Bond.SINGLE.edge(orgEdge.either(),
                                                   orgEdge.other(orgEdge.either())));
            else
                localised.addEdge(orgEdge);
        }

        // verify hydrogen count
        for (int v = 0; v < delocalised.order(); v++) {
            if (localised.implHCount(v) != delocalised.implHCount(v)) {
                throw new InvalidSmilesException("cannot assign localised bonds to structure - valence error");
            }
        }
    }


    IntSet buildSet(Graph g) {

        BitSet undecided = new BitSet(g.order());
        for (int v = 0; v < g.order(); v++) {

            if (!g.atom(v).aromatic())
                continue;

            // if the all aromatic bonds are single
            if (!predetermined(g, v)) {
                undecided.set(v);
            }
        }

        return IntSet.fromBitSet(undecided);
    }

    boolean predetermined(Graph g, int v) {

        Atom a = delocalised.atom(v);

        int q   = a.charge();
        int deg = delocalised.degree(v) + delocalised.implHCount(v);
        
        for (Edge e : delocalised.edges(v)) {
            if (e.bond() == Bond.DOUBLE) {
                if (q == 0 
                        && (a.element() == Element.Nitrogen || (a.element() == Element.Sulfur && deg > 3 )) 
                        && g.atom(e.other(v)).element() == Element.Oxygen)
                    return false;
                return true;
            }
        }

        

        switch (a.element()) {
            case Carbon:
                if (q == 1 && deg == 3) return true;
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
                else
                    return false;
        }

        return false;
    }


    static Graph localise(Graph delocalised) throws InvalidSmilesException {
        return new Localise(delocalised).localised;
    }
}
