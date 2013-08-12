package uk.ac.ebi.grins;

import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

/** @author John May */
public class FunctionsTest {

    @Test public void reverse_ethanol() throws Exception {
        ChemicalGraph g = ChemicalGraph.fromSmiles("CCO");
        assertThat(Functions.reverse(g).toSmiles(),
                   is("OCC"));
    }

    @Test public void reverse_withBranch() throws Exception {
        ChemicalGraph g = ChemicalGraph.fromSmiles("CC(CC(CO)C)CCO");
        assertThat(Functions.reverse(g).toSmiles(),
                   is("OCCC(CC(C)CO)C"));
    }
}