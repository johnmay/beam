package uk.ac.ebi.beam;

import java.util.Arrays;
import java.util.BitSet;

/**
 * Simple matching greedily chooses edges and matches. The produced matching is not guaranteed to be
 * maximum but provides a starting point for improvement through augmentation.
 *
 * @author John May
 */
final class ArbitraryMatching {

    /**
     * Create an arbitrary matching on the subset of vertices ('s') of provided graph. The provided
     * matching should be empty.
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
                if ((e.bond() != Bond.SINGLE) && m.unmatched(w) && s.get(w)) {
                    m.match(v, w);
                    nMatched += 2;
                    break;
                }
            }
        }

        return nMatched;
    }

    /**
     * When precisely two vertices are unmatched we only need to find a single augmenting path. 
     * Rather than run through edmonds with blossoms etc we simple do a targest DFS for the path.
     * 
     * @param g graph
     * @param m matching
     * @param nMatched current matching cardinality must be |s|-nMathced == 2
     * @param s subset size
     * @return new match cardinality 
     */
    static int augmentOnce(final Graph g, final Matching m, int nMatched, final BitSet s) {

        int vStart = s.nextSetBit(0);
        while (vStart >= 0) {
            if (!m.matched(vStart)) break;
            vStart = s.nextSetBit(vStart + 1);
        }
        int vEnd = s.nextSetBit(vStart + 1);
        while (vEnd >= 0) {
            if (!m.matched(vEnd)) break;
            vEnd = s.nextSetBit(vEnd + 1);
        }

        // find an augmenting path between vStart and vEnd
        int[] path = new int[g.order()];
        int len = findPath(g, vStart, vEnd, s, path, 0, m, false);
        if (len > 0) {
            // augment
            for (int i = 0; i < len; i += 2) {
                m.match(path[i], path[i + 1]);
            }
            nMatched += 2;
        }

        return nMatched;
    }

    static int findPath(Graph g, int v, int end, BitSet unvisited, int[] path, int len, Matching m, boolean matchNeeded) {
        unvisited.clear(v);
        path[len++] = v;
        int l;
        for (final Edge e : g.edges(v)) {
            int w = e.other(v);
            if (unvisited.get(w)) {
                if (w == end) {
                    path[len] = w;
                    len++;
                    unvisited.set(v);
                    // odd length path no good
                    return ((len & 0x1) == 1) ? 0 : len;
                }
                else if ((m.other(w) == v) == matchNeeded) {
                    if ((l = findPath(g, w, end, unvisited, path, len, m, !matchNeeded)) > 0) {
                        unvisited.set(v);
                        return l;
                    }
                }
            }
        }
        unvisited.set(v);
        return 0;
    }
}
