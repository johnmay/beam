package uk.ac.ebi.grins;

/**
 * Defines a function which can be applied to one type to produce another.
 *
 * @author John May
 */
public interface Function<S, T> {

    /**
     * Apply the function to an instance of 's' and producing and instance 't'.
     *
     * @param s input instance
     * @return output instance
     */
    T apply(S s);

}
