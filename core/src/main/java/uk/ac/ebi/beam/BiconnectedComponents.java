package uk.ac.ebi.beam;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Deque;
import java.util.List;

/**
 * see. http://en.wikipedia.org/wiki/Biconnected_component
 *
 * @author John May
 */
final class BiconnectedComponents {

    private int[] d, low;

    private final Graph  g;
    private final Edge[] stack;
    private int nstack = 0;

    private final List<List<Edge>> components = new ArrayList<List<Edge>>(2);

    private final BitSet cyclic = new BitSet();

    int count = 0;

    BiconnectedComponents(Graph g) {
        this.d     = new int[g.order()];
        this.low   = new int[g.order()];
        this.g     = g;
        this.stack = new Edge[g.size()];
        for (int u = 0; u < g.order(); u++) {
            if (d[u] == 0)
                visit(u, u);
        }
        low = null;
        d = null;
    }

    private void visit(int u, int p) {
        low[u] = d[u] = ++count;
        for (final Edge e : g.edges(u)) {
            final int v = e.other(u);
            if (d[v] == 0) {
                stack[nstack++] = e;
                visit(v, u);
                if (low[v] == d[u]) {
                    store(e);
                }
                else if (low[v] > d[u]) {
                    --nstack;
                }
                if (low[v] < low[u])
                    low[u] = low[v];
            }
            else if (v != p && d[v] < d[u]) {
                // back edge
                stack[nstack++] = e;
                if (d[v] < low[u])
                    low[u] = d[v];
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
