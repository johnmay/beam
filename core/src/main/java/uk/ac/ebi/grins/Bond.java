package uk.ac.ebi.grins;

/**
 * Enumeration of valid connections between atoms. The connections include all
 * the valid undirected and directed bond types and {@link #DOT}. Opposed to the
 * other types, {@link #DOT} indicates that two atoms are not connected. <p/>
 *
 * @author John May
 */
public enum Bond {

    /** Atoms are not bonded. */
    DOT(".", 0),

    /** Atoms are bonded by either a single or aromatic bond. */
    IMPLICIT("", 0) {
        @Override public int electrons() {
            throw new IllegalArgumentException("unknown number of electrons in implied bond");
        }
    },

    /** Atoms are bonded by a single pair of electrons. */
    SINGLE("-", 2),

    /** Atoms are bonded by two pairs of electrons. */
    DOUBLE("=", 4),

    /** Atoms are bonded by three pairs of electrons. */
    TRIPLE("#", 6),

    /** Atoms are bonded by four pairs of electrons. */
    QUADRUPLE("$", 8),

    /** Atoms are bonded by a delocalized bond of an aromatic system. */
    AROMATIC(":", 3),

    /**
     * Directional, single or aromatic bond (currently always single). The bond
     * is relative to each endpoint such that the second endpoint is
     * <i>above</i> the first or the first end point is <i>below</i> the
     * second.
     */
    UP("/", 2) {
        @Override public Bond inverse() {
            return DOWN;
        }
    },

    /**
     * Directional, single or aromatic bond (currently always single). The bond
     * is relative to each endpoint such that the second endpoint is
     * <i>below</i> the first or the first end point is <i>above</i> the
     * second.
     */
    DOWN("\\", 2) {
        @Override public Bond inverse() {
            return UP;
        }
    };

    /** The symbol for the bond in the SMILES grammar. */
    private final String symbol;

    /** The total number of electrons shared, i.e. not the number of pairs. */
    private final int electrons;

    private Bond(String symbol, int electrons) {
        this.symbol = symbol;
        this.electrons = electrons;
    }


    /**
     * The symbol of the bond in the SMILES grammar.
     *
     * @return bond symbol
     */
    public final String symbol() {
        return symbol;
    }

    /**
     * The total number electrons shared between atoms.
     *
     * @return number of electrons
     */
    public int electrons() {
        return electrons;
    }

    /**
     * Access the inverse of a directional bond ({@link #UP}, {@link #DOWN}). If
     * a bond is non-directional the same bond is returned.
     *
     * @return inverse of the bond
     */
    public Bond inverse() {
        return this;
    }

    /** @inheritDoc */
    public final String toString() {
        return symbol;
    }
}
