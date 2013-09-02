package uk.ac.ebi.beam;

import org.hamcrest.collection.IsIterableWithSize;
import org.junit.Test;

import static org.hamcrest.CoreMatchers.hasItems;
import static org.hamcrest.CoreMatchers.not;
import static org.junit.Assert.assertThat;

/** @author John May */
public class MatchingTest {

    @Test public void empty() throws Exception {
        Matching matching = Matching.empty(Graph.fromSmiles("CCCCC"));
        assertThat(matching.matches(),
                   IsIterableWithSize.<Tuple>iterableWithSize(0));
    }

    @Test public void basic() throws Exception {
        Matching matching = Matching.empty(Graph.fromSmiles("CCCCC"));
        matching.match(0, 1);
        matching.match(2, 3);
        assertThat(matching.matches(),
                   IsIterableWithSize.<Tuple>iterableWithSize(2));
        assertThat(matching.matches(),
                   hasItems(Tuple.of(0, 1),
                            Tuple.of(2, 3)));
    }

    @Test public void adjusted() throws Exception {
        Matching matching = Matching.empty(Graph.fromSmiles("CCCCC"));
        matching.match(0, 1);
        matching.match(2, 3);
        matching.match(1, 2); // 0-1 and 2-3 should not be 

        assertThat(matching.matches(),
                   not(hasItems(Tuple.of(0, 1),
                                Tuple.of(2, 3))));
        assertThat(matching.matches(),
                   IsIterableWithSize.<Tuple>iterableWithSize(1));
        assertThat(matching.matches(),
                   hasItems(Tuple.of(1, 2)));
    }

}
