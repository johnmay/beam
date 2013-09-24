package uk.ac.ebi.beam;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

/**
 * Utility to localise aromatic bonds.
 *
 * @author John May
 */
final class Localise {

    private final Graph delocalised, localised;
    private final BitSet subset;
    private final Map<Edge, Edge> edgeAssignments = new HashMap<Edge, Edge>();

    private Localise(Graph delocalised, BitSet subset) throws InvalidSmilesException {

        this.delocalised = delocalised;
        this.localised   = new Graph(delocalised.order());
        this.subset      = subset;

        // make initial matching - then improve it
        Matching m = MaximumMatching.maximise(delocalised,
                                              ArbitraryMatching.of(delocalised, subset),
                                              IntSet.fromBitSet(subset));

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
            else if (orgEdge.bond() == Bond.AROMATIC || orgEdge.bond() == Bond.SINGLE)
                localised.addEdge(Bond.IMPLICIT.edge(orgEdge.either(),
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


    static BitSet buildSet(Graph g) {

        BitSet undecided = new BitSet(g.order());
        
        for (int v = 0; v < g.order(); v++) {
            if (g.atom(v).aromatic()) {
                undecided.set(v, !predetermined(g, v));
            }
        }

        return undecided;
    }

    static boolean predetermined(Graph g, int v) {

        Atom a = g.atom(v);

        int q   = a.charge();
        int deg = g.degree(v) + g.implHCount(v);

        for (Edge e : g.edges(v)) {
            if (e.bond() == Bond.DOUBLE) {
                if (q == 0
                        && (a.element() == Element.Nitrogen || (a.element() == Element.Sulfur && deg > 3))
                        && g.atom(e.other(v)).element() == Element.Oxygen)
                    return false;
                return true;
            }
            // triple or quadruple bond - we don't need to assign anymore p electrons
            else if (e.bond().order() > 2) {
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
        BitSet subset = buildSet(delocalised);         
        if (hasOddCardinality(subset))
            throw new InvalidSmilesException("No localised structure can be assigned.");
        return subset.isEmpty() ? delocalised 
                                : new Localise(delocalised, subset).localised;
    }

    private static boolean hasOddCardinality(BitSet s) {
        return (s.cardinality() & 0x1) == 1;
    }
}
