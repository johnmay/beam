package uk.ac.ebi.beam;

import java.util.ArrayList;
import java.util.Arrays;
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
    private final BitSet simple = new BitSet();

    int count = 0;

    private int mark = 1;

    BiconnectedComponents(Graph g) {
        this(g, true);
    }

    BiconnectedComponents(Graph g, boolean storeComponents) {
        this.depth = new int[g.order()];
        this.low = new int[g.order()];
        this.g = g;
        this.stack  = new Edge[g.size()];
        if (storeComponents) {
            for (int u = 0; u < g.order(); u++) {
                if (depth[u] == 0)
                    visitWithComp(u, null);
            }
        } else {
            for (int u = 0; u < g.order(); u++) {
                if (depth[u] == 0)
                    visit(u, null);
            }
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

    private void visitWithComp(final int u, final Edge from) {
        low[u] = depth[u] = ++count;
        int j = g.degree(u);
        while (--j>=0) {

            final Edge e = g.edgeAt(u, j);
            if (e==from) continue;

            final int v = e.other(u);
            if (depth[v] == 0) {
                stack[nstack] = e;
                ++nstack;
                visitWithComp(v, e);
                if (low[v] == depth[u])
                    storeWithComp(e);
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
    
    boolean inCycle(int v) {
        return cyclic.get(v);
    }

    boolean inSimpleCycle(int v) {
        return simple.get(v);
    }

    private void store(Edge e) {
        Edge f;

        final BitSet tmp = new BitSet();
        
        // count the number of unique vertices and edges
        int numEdges = 0;
        
        do {
            f = stack[--nstack];
            int v = f.either();
            int w = f.other(v);

            tmp.set(v);
            tmp.set(w);
            
            numEdges++;
        } while (f != e);
        
        cyclic.or(tmp);
        
        if (tmp.cardinality() == numEdges)
            simple.or(tmp);
    }

    private void storeWithComp(Edge e) {
        List<Edge> component = new ArrayList<Edge>(6);
        Edge f;

        final BitSet tmp = new BitSet();

        // count the number of unique vertices and edges
        int numEdges = 0;

        do {
            f = stack[--nstack];
            int v = f.either();
            int w = f.other(v);

            tmp.set(v);
            tmp.set(w);

            component.add(f);
            numEdges++;
        } while (f != e);

        cyclic.or(tmp);

        if (tmp.cardinality() == numEdges)
            simple.or(tmp);

        components.add(Collections.unmodifiableList(component));
    }

    public List<List<Edge>> components() {
        return Collections.unmodifiableList(components);
    }

    BitSet cyclic() {
        return cyclic;
    }
}
