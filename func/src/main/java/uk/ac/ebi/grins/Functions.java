package uk.ac.ebi.grins;

import java.util.Random;

/**
 * Collection of utilities for transforming chemical graphs.
 *
 * @author John May
 */
public final class Functions {

    /// non-instantiable
    private Functions() {
    }

    /**
     * Randomise the atom order of the provided chemical graph.
     *
     * @param g chemical graph
     * @return a copy of the original graph with the order of the atoms
     *         randomised
     */
    public static ChemicalGraph randomise(ChemicalGraph g) {
        return g.permute(random(g.order()));
    }

    /**
     * Invert the atom order of the provided chemical graph.
     *
     * @param g chemical graph
     * @return a copy of the original graph with the order of the atoms
     *         inverted
     */
    public static ChemicalGraph inverse(ChemicalGraph g) {
        return g.permute(inv(ident(g.order())));
    }

    private static int[] ident(int n) {
        int[] p = new int[n];
        for (int i = 0; i < n; i++)
            p[i] = i;
        return p;
    }

    private static int[] random(int n) {
        int[] p = ident(n);
        Random rnd = new Random();
        for (int i = n; i > 1; i--)
            swap(p, i - 1, rnd.nextInt(i));
        return p;
    }

    private static int[] inv(int[] p) {
        int[] q = p.clone();
        for (int i = 0; i < p.length; i++)
            q[p[i]] = i;
        return q;
    }

    private static void swap(int[] p, int i, int j) {
        int tmp = p[i];
        p[i] = p[j];
        p[j] = tmp;
    }
}
