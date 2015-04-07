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

    private boolean[] visited;
    private int[]     d, low;

    private final Graph       g;
    private final Deque<Edge> stack;

    private final List<List<Edge>> components = new ArrayList<List<Edge>>(2);

    private final BitSet cyclic = new BitSet();

    int count = 0;

    BiconnectedComponents(Graph g) {
        this.visited = new boolean[g.order()];
        this.d = new int[g.order()];
        this.low = new int[g.order()];
        this.g = g;
        this.stack = new ArrayDeque<Edge>();
        for (int u = 0; u < g.order(); u++) {
            if (!visited[u])
                visit(u, u);
        }
        low = null;
        d = null;
    }

    private void visit(int u, int p) {
        visited[u] = true;
        low[u] = d[u] = ++count;
        for (Edge e : g.edges(u)) {
            int v = e.other(u);
            if (!visited[v]) {
                stack.push(e);
                visit(v, u);
                if (low[v] >= d[u])
                    store(e);
                low[u] = Math.min(low[u], low[v]);
            }
            else if (v != p && d[v] < d[u]) {
                // back edge
                stack.push(e);
                low[u] = Math.min(low[u], d[v]);
            }
        }
    }

    private void store(Edge e) {
        if (stack.peek() != e) {
            List<Edge> component = new ArrayList<Edge>(6);
            Edge f = null;
            do {
                f = stack.pop();
                markCyclic(f);
                component.add(f);
            } while (f!=e);
            components.add(Collections.unmodifiableList(component));
        } else {
            stack.pop();
        }
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
