package uk.ac.ebi.beam;

import java.util.BitSet;

/**
 * Simple matching greedily chooses edges and matches. The produced matching is
 * not guaranteed to be maximum but provides a starting point for improvement
 * through augmentation.
 *
 * @author John May
 */
final class ArbitraryMatching {

    /**
     * Create an arbitrary matching on the subset of vertices ('s') of provided
     * graph. The provided matching should be empty.  
     *
     * @param g graph to match
     * @param m empty matching (presumed)
     * @param s subset of vertices
     * @return number of vertices matched
     */
    static int initial(final Graph g, Matching m, final BitSet s) {

        int nMatched = 0;
        
        for (int v = s.nextSetBit(0); v >= 0; v = s.nextSetBit(v + 1)) {

            // skip if already matched
            if (m.matched(v))
                continue;

            // find a single edge which is not matched and match it
            for (final Edge e : g.edges(v)) {
                int w = e.other(v);
                if (m.unmatched(w) && s.get(w)) {
                    m.match(v, w);
                    nMatched += 2;
                    break;
                }
            }
        }

        return nMatched;
    }
}
