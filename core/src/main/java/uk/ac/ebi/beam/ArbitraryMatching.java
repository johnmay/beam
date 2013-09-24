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
     * graph.
     *
     * @param g graph to match
     * @param s subset of vertices
     * @return non-maximal matching
     */
    static Matching of(final Graph g, final BitSet s) {

        final Matching m = Matching.empty(g);

        for (int v = s.nextSetBit(0); v >= 0; v = s.nextSetBit(v + 1)) {

            // skip if already matched
            if (!m.unmatched(v))
                continue;

            // find a single edge which is not matched and match it
            for (Edge e : g.edges(v)) {
                int w = e.other(v);

                if (!s.get(w) || !m.unmatched(w))
                    continue;

                m.match(v, w);
                break;
            }
        }

        return m;
    }

}
