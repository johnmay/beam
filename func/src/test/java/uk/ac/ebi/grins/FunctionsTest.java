package uk.ac.ebi.grins;

import org.junit.Test;

import java.io.IOException;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

/** @author John May */
public class FunctionsTest {

    @Test public void reverse_ethanol() throws Exception {
        Graph g = Graph.fromSmiles("CCO");
        assertThat(Functions.reverse(g).toSmiles(),
                   is("OCC"));
    }

    @Test public void reverse_withBranch() throws Exception {
        Graph g = Graph.fromSmiles("CC(CC(CO)C)CCO");
        assertThat(Functions.reverse(g).toSmiles(),
                   is("OCCC(CC(C)CO)C"));
    }

    @Test public void atomBasedDBStereo() throws Exception {
        Graph g = Graph.fromSmiles("F/C=C/F");
        assertThat(Functions.atomBasedDBStereo(g).toSmiles(),
                   is("F[C@H]=[C@@H]F"));
    }

    @Test public void bondBasedDBStereo() throws Exception {
        Graph g = Graph.fromSmiles("F[C@H]=[C@@H]F");
        assertThat(Functions.bondBasedDBStereo(g).toSmiles(),
                   is("F/C=C/F"));
    }

    @Test public void canoncalise() throws IOException {
        Graph g = Graph.fromSmiles("CCOCC");
        Graph h = Functions.canonicalize(g,
                                         new long[]{56, 67, 3, 67, 56});
        assertThat(h.toSmiles(),
                   is("O(CC)CC"));
    }
}
