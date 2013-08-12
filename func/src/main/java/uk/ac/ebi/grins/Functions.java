package uk.ac.ebi.grins;

import java.util.Random;

/**
 * Collection of utilities for transforming chemical graphs.
 *
 * @author John May
 */
public final class Functions {

    // convert to atom-based double-bond configurations
    private static final ToTrigonalTopology ttt = new ToTrigonalTopology();

    // convert to bond-based double-bond configuration
    private static final FromTrigonalTopology ftt = new FromTrigonalTopology();

    // bond label conversion -> to implicit
    private static final ExplicitToImplicit eti = new ExplicitToImplicit();

    // bond label conversion -> to explicit
    private static final ImplicitToExplicit ite = new ImplicitToExplicit();

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
     * Reverse the atom order of the provided chemical graph.
     *
     * @param g chemical graph
     * @return a copy of the original graph with the order of the atoms
     *         reversed
     */
    public static ChemicalGraph reverse(ChemicalGraph g) {
        return g.permute(reverse(g.order()));
    }

    /**
     * Convert any directional bond based stereo configuration to atom-based
     * specification.
     *
     * @param g chemical graph graph
     * @return a copy of the original graph but with directional bonds removed
     *         and atom-based double-bond stereo configruation.
     */
    public static ChemicalGraph atomBasedDBStereo(ChemicalGraph g) {
        return eti.apply(ttt.apply(ite.apply(g)));
    }

    /**
     * Convert a graph with atom-based double-bond stereo configuration to
     * bond-based specification (direction UP and DOWN bonds).
     *
     * @param g chemical graph graph
     * @return a copy of the original graph but with bond-based stereo-chemistry
     */
    public static ChemicalGraph bondBasedDBStereo(ChemicalGraph g) {
        return eti.apply(ftt.apply(ite.apply(g)));
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

    private static int[] reverse(int n) {
        int[] p = new int[n];
        for (int i = 0; i < n; i++)
            p[i] = n - i - 1;
        return p;
    }

    // inverse of permutation
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
