package uk.ac.ebi.grins;

import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.not;
import static org.hamcrest.CoreMatchers.sameInstance;
import static org.junit.Assert.assertThat;
import static uk.ac.ebi.grins.Configuration.TH1;
import static uk.ac.ebi.grins.Configuration.TH2;

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

    @Test
    public void permutationParity() {
        assertThat(Topology.parity(new int[]{0, 1, 2, 3},
                                   new int[]{0, 1, 2, 3}), is(1));   // even
        assertThat(Topology.parity(new int[]{0, 1, 2, 3},
                                   new int[]{0, 1, 3, 2}), is(-1));  // swap 2,3 = odd
        assertThat(Topology.parity(new int[]{0, 1, 2, 3},
                                   new int[]{1, 0, 3, 2}), is(1));   // swap 0,1 = even
        assertThat(Topology.parity(new int[]{0, 1, 2, 3},
                                   new int[]{2, 0, 3, 1}), is(-1));  // swap 0,3 = odd
    }

    @Test
    public void sort() {
        int[] org = new int[]{1, 2, 3, 4};
        assertThat(Topology.sort(org, new int[]{0, 1, 2, 3, 4}),
                   is(not(sameInstance(org))));
        assertThat(Topology.sort(org, new int[]{0, 1, 2, 3, 4}),
                   is(new int[]{1, 2, 3, 4}));
        assertThat(Topology.sort(org, new int[]{0, 2, 1, 3, 4}),
                   is(new int[]{2, 1, 3, 4}));
        // non-sequential
        assertThat(Topology.sort(org, new int[]{0, 2, 1, 7, 4}),
                   is(new int[]{2, 1, 4, 3}));
    }

    @Test public void tetrahedralAtom() {
        Topology t1 = Topology.tetrahedral(1, new int[]{0, 2, 3, 4}, TH1);
        assertThat(t1.atom(), is(1));
    }

    @Test public void tetrahedralOrderBy() {
        // test shows the first example of tetrahedral configuration from the
        // OpenSMILES specification

        // N=1, Br=2, O=3, C=4
        Topology t1 = Topology.tetrahedral(0, new int[]{1, 2, 3, 4}, TH1);

        // N, Br, O, C
        assertThat(t1.orderBy(new int[]{0, 1, 2, 4, 3})
                     .configuration(), is(TH2));
        // O, Br, C, N
        assertThat(t1.orderBy(new int[]{0, 4, 2, 1, 3})
                     .configuration(), is(TH1));
        // C, Br, N, O
        assertThat(t1.orderBy(new int[]{0, 3, 2, 4, 1})
                     .configuration(), is(TH1));
        // C, Br, O, N
        assertThat(t1.orderBy(new int[]{0, 4, 2, 3, 1})
                     .configuration(), is(TH2));
        // Br, O, N, C
        assertThat(t1.orderBy(new int[]{0, 3, 1, 2, 4})
                     .configuration(), is(TH1));
        // Br, C, O, N
        assertThat(t1.orderBy(new int[]{0, 4, 1, 3, 2})
                     .configuration(), is(TH1));
        // Br, N, C, O
        assertThat(t1.orderBy(new int[]{0, 2, 1, 4, 3})
                     .configuration(), is(TH1));
        // Br, N, O, C
        assertThat(t1.orderBy(new int[]{0, 2, 1, 3, 4})
                     .configuration(), is(TH2));
    }

}
