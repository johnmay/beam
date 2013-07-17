package uk.ac.ebi.grins;

import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.sameInstance;
import static org.junit.Assert.assertThat;

/** @author John May */
public class TopologyTest {

    @Test
    public void unknown() throws Exception {
        assertThat(Topology.unknown()
                           .configuration(),
                   is(Configuration.UNKNOWN));
    }

    @Test(expected = IllegalArgumentException.class)
    public void unknownAtom() throws Exception {
        Topology.unknown()
                .atom();
    }

    @Test
    public void unknownTransform() {
        assertThat(Topology.unknown().transform(new int[0]),
                   is(sameInstance(Topology.unknown())));
    }

    @Test
    public void unknownOrderBy() {
        assertThat(Topology.unknown().orderBy(new int[0]),
                   is(sameInstance(Topology.unknown())));
    }
}
