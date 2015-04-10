package uk.ac.ebi.beam;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;

/**
 * see. http://en.wikipedia.org/wiki/Biconnected_component
 *
 * @author John May
 */
final class BiconnectedComponents {

    private int[] depth, low;

    private final Graph  g;
    private final Edge[] stack;
    private int nstack = 0;

    private final List<List<Edge>> components = new ArrayList<List<Edge>>(2);

    private final BitSet cyclic = new BitSet();

    int count = 0;

    BiconnectedComponents(Graph g) {
        this.depth = new int[g.order()];
        this.low   = new int[g.order()];
        this.g     = g;
        this.stack = new Edge[g.size()];
        for (int u = 0; u < g.order(); u++) {
            if (depth[u] == 0)
                visit(u, null);
        }
        low = null;
        depth = null;
    }

    private void visit(final int u, final Edge from) {
        low[u] = depth[u] = ++count;
        int j = g.degree(u);
        while (--j>=0) {
            
            final Edge e = g.edgeAt(u, j);
            if (e==from) continue;
            
            final int v = e.other(u);
            if (depth[v] == 0) {
                stack[nstack] = e;
                ++nstack;
                visit(v, e);
                if (low[v] == depth[u])
                    store(e);
                else if (low[v] > depth[u])
                    --nstack;
                if (low[v] < low[u])
                    low[u] = low[v];
            }
            else if (depth[v] < depth[u]) {
                // back edge
                stack[nstack] = e;
                ++nstack;
                if (depth[v] < low[u])
                    low[u] = depth[v];
            }
        }
    }

    private void store(Edge e) {
        List<Edge> component = new ArrayList<Edge>(6);
        Edge f;
        do {
            f = stack[--nstack];
            markCyclic(f);
            component.add(f);
        } while (f != e);
        components.add(Collections.unmodifiableList(component));
    }

    private void markCyclic(Edge f) {
        cyclic.set(f.either());
        cyclic.set(f.other(f.either()));
    }

    public List<List<Edge>> components() {
        return Collections.unmodifiableList(components);
    }

    BitSet cyclic() {
        return cyclic;
    }
}
