package uk.ac.ebi.beam;

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
    static Matching of(final Graph g, final boolean[] s) {

        final Matching m = Matching.empty(g);

        for (int v = 0; v < g.order(); v++) {

            // skip if not in the subset of vertices or already matched
            if (!s[v] || !m.unmatched(v))
                continue;

            // find a single edge which is not matched and match it
            for (Edge e : g.edges(v)) {
                int w = e.other(v);

                if (!s[w] || !m.unmatched(w))
                    continue;

                m.match(v, w);
                break;
            }
        }

        return m;
    }

}
