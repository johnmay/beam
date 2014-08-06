package uk.ac.ebi.beam;

/** @author John May */
public final class PartitionOld {

    /* Nil value (-1) */
    private final int Nil = -1;

    /** The color of each element. */
    final int[] colors;

    /** Ordering of cells. */
    final int[] ordering, aux;

    /** Indices of non-discrete cells in 'order'. */
    final int[] cells;

    /** Indicates whether an index is discrete. */
    boolean[] discrete;

    /** Number of cells. */
    int nCells;

    PartitionOld(int n) {
        colors   = new int[n];
        ordering = new int[n];
        aux      = new int[n];
        cells    = new int[n];
        discrete = new boolean[n];

        // initial all entries are equivalent (labeled 1)
        nCells = 1;
        for (int i = 0; i < n; i++)
            ordering[i] = i;
        
        if (n < 2)
            return;
        
        cells[0] = 0;
        cells[1] = n - 1;
        if (n > 2)
            cells[2] = Nil;
    }

    /**
     * Obtain the color of the value.
     * 
     * @param v value
     * @return 
     */
    int colorOf(int v) {
        return colors[v];
    }

    /**
     * Refine the non-discrete cells in the partition using the provided comparator.
     *
     * @param cmp values
     * @return the number of ranks
     */
    boolean refine(Comparison cmp) {

        int n = 0;
        int nCellsOld = nCells;

        for (int i = 0; i < cells.length && cells[i] != Nil; i += 2) {
            int lo = cells[i];
            int hi = cells[i + 1];

            sort(ordering, lo, 1 + hi - lo, cmp);

            aux[n++] = lo;

            // assign new cells (separate array - ordering_aux)
            int u = ordering[lo];
            int base = colors[u];
            int col = base;
            for (int j = lo + 1; j <= hi; j++) {
                int v = ordering[j];

                // the values u and v (at index j-1 and j) are different
                if (cmp.less(u, v)) {

                    if (j - 1 - aux[n - 1] == 0) {
                        discrete[ordering[j - 1]] = true;
                        aux[n - 1] = j;
                    }
                    else {
                        aux[n++] = j - 1;
                        aux[n++] = j;
                    }

                    nCells++;
                    col = base + (j - lo);
                }

                colors[v] = col;
                u = v;
            }

            if (aux[n - 1] == hi) {
                discrete[ordering[hi]] = true;
                n--;
            }
            else {
                aux[n++] = hi;
            }

        }

        // early termination
        if (n < aux.length)
            aux[n++] = Nil;

        System.arraycopy(aux, 0, cells, 0, n);

        return nCells > nCellsOld;
    }
    
    /**
     * A partition is the unit partition if all elements are equivalent.
     * 
     * @return the partition is the unit partition (or not).
     */
    boolean unit() {
        return len(0) == colors.length;
    }

    /**
     * A partition is discrete if the each element is in a separate.
     *
     * @return the partition is discrete (or not).
     */
    boolean discrete() {
        return nCells == colors.length;
    }

    /**
     * Is the value is in a discrete.
     *
     * @param v value
     * @return the value is in a discrete cell
     */
    boolean discrete(int v) {
        return discrete[v];
    }

    /**
     * Obtain the values in the cell.
     *
     * @param cell index of non-discrete cell
     * @return the items in the cell
     */
    int[] cell(int cell) {
        int lo = cells[2 * cell];
        int hi = cells[2 * cell + 1];
        int[] xs = new int[1 + hi - lo];
        System.arraycopy(ordering, lo, xs, 0, 1 + hi - lo);
        return xs;
    }

    /**
     * Length of the non-discrete {@code cell}.
     *
     * @param cell index of cell
     * @return the number of items in the cell (they are equivalent)
     */
    int len(int cell) {
        return 1 + cells[2 * cell + 1] - cells[2 * cell];
    }

    /**
     * Split a non-discrete cell at a given position. Only non-discrete cells
     * are index and so '0' is always the lowest ranking non-discrete cell.
     *
     * @param cell cell index (start at 0)
     * @param i    the index of the element to split (0 <= i < len(cell))
     */
    void split(int cell, int i) {

        // lo and hi indices into the ordering array 
        int lo = cells[2 * cell];
        int hi = cells[2 * cell + 1];

        int col = 1 + colors[ordering[lo]];
        int split = lo + i;

        if (split > hi)
            throw new IllegalStateException("can not split at index " + i + " not in cell");

        for (int j = hi; j > split; j--) {
            colors[ordering[j]] = col;
        }

        for (int j = split; j > lo; j--) {
            exch(ordering, j - 1, j);
            colors[ordering[j]] = col;
        }

        cells[2 * cell]++;
        nCells++;
    }
    
    int[] permutation() {
        return ordering;
    }

    /**
     * Exchence the element at index i and j.
     *
     * @param xs array of value
     * @param i  first index
     * @param j  second index
     */
    private static void exch(int[] xs, int i, int j) {
        int tmp = xs[i];
        xs[i] = xs[j];
        xs[j] = tmp;
    }

    @Override public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[{").append(ordering[0]);
        for (int i = 1; i < ordering.length; i++) {
            if (colors[ordering[i]] != colors[ordering[i - 1]])
                sb.append("}, {");
            else
                sb.append(", ");
            sb.append(ordering[i]);

        }
        sb.append("}]");
        return sb.toString();
    }

    /**
     * Sort the values (using merge sort) in {@code vs} from {@code lo} (until
     * {@code len}) by the {@code prev[]} and then {@code curr[]} invariants to
     * determine rank. The values in {@code vs} are indices into the invariant
     * arrays.
     *
     * @param vs  values (indices)
     * @param lo  the first value to start sorting from
     * @param len the len of values to consider
     * @param cmp comparison
     */
    void sort(int[] vs, int lo, int len, Comparison cmp) {

        if (len < 20) {
            insertionSort(vs, lo, len, cmp);
            return;
        }

        int split = len / 2;

        sort(vs, lo, split, cmp);
        sort(vs, lo + split, len - split, cmp);

        // sub arrays already sorted, no need to merge
        if (!cmp.less(vs[lo + split], vs[lo + split - 1]))
            return;

        merge(vs, lo, split, len, cmp);
    }

    /**
     * Merge the values which are sorted between {@code lo} - {@code split} and
     * {@code split} - {@code len}.
     *
     * @param vs    vertices
     * @param lo    start index
     * @param split the middle index (partition)
     * @param len   the range to merge
     * @param cmp   comparison
     */
    private void merge(int[] vs, int lo, int split, int len, Comparison cmp) {
        System.arraycopy(vs, lo, aux, lo, len);

        int i = lo, j = lo + split;
        int iMax = lo + split, jMax = lo + len;
        for (int k = lo, end = lo + len; k < end; k++) {
            if (i == iMax)
                vs[k] = aux[j++];
            else if (j == jMax)
                vs[k] = aux[i++];
            else if (cmp.less(aux[i], aux[j]))
                vs[k] = aux[i++];
            else
                vs[k] = aux[j++];
        }
    }

    static void insertionSort(int[] vs, int lo, int len, Comparison cmp) {
        for (int j = lo + 1, hi = lo + len; j < hi; j++) {
            int v = vs[j];
            int i = j - 1;
            while ((i >= lo) && cmp.less(v, vs[i]))
                vs[i + 1] = vs[i--];
            vs[i + 1] = v;
        }
    }

    static Comparison cmpWithInv(final long[] inv) {
        return new InvComparison(inv);
    }

    private static final class InvComparison implements Comparison {
        long[] inv;

        private InvComparison(long[] inv) {
            this.inv = inv;
        }

        @Override public boolean less(int i, int j) {
            return inv[i] < inv[j];
        }
    }

    static interface Comparison {
        boolean less(int i, int j);
    }
}
