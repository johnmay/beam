package uk.ac.ebi.grins;

import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;
import static uk.ac.ebi.grins.Bond.AROMATIC;
import static uk.ac.ebi.grins.Bond.DOT;
import static uk.ac.ebi.grins.Bond.DOUBLE;
import static uk.ac.ebi.grins.Bond.DOWN;
import static uk.ac.ebi.grins.Bond.IMPLICIT;
import static uk.ac.ebi.grins.Bond.QUADRUPLE;
import static uk.ac.ebi.grins.Bond.SINGLE;
import static uk.ac.ebi.grins.Bond.TRIPLE;
import static uk.ac.ebi.grins.Bond.UP;

/** @author John May */
public class BondTest {

    @Test public void dotElectrons() throws Exception {
        assertThat(DOT.electrons(), is(0));
    }

    @Test public void singleElectrons() throws Exception {
        assertThat(SINGLE.electrons(), is(2));
    }

    @Test public void doubleElectrons() throws Exception {
        assertThat(DOUBLE.electrons(), is(4));
    }

    @Test public void tripleElectrons() throws Exception {
        assertThat(TRIPLE.electrons(), is(6));
    }

    @Test public void quadrupleElectrons() throws Exception {
        assertThat(QUADRUPLE.electrons(), is(8));
    }

    @Test public void aromaticElectrons() throws Exception {
        assertThat(AROMATIC.electrons(), is(3));
    }

    @Test public void upElectrons() throws Exception {
        assertThat(UP.electrons(), is(2));
    }

    @Test public void downElectrons() throws Exception {
        assertThat(DOWN.electrons(), is(2));
    }

    @Test(expected = IllegalArgumentException.class)
    public void implicitElectrons() throws Exception {
        IMPLICIT.electrons();
    }

    @Test public void dotInverse() throws Exception {
        assertThat(DOT.inverse(), is(DOT));
    }

    @Test public void singleInverse() throws Exception {
        assertThat(SINGLE.inverse(), is(SINGLE));
    }

    @Test public void doubleInverse() throws Exception {
        assertThat(DOUBLE.inverse(), is(DOUBLE));
    }

    @Test public void tripleInverse() throws Exception {
        assertThat(TRIPLE.inverse(), is(TRIPLE));
    }

    @Test public void quadrupleInverse() throws Exception {
        assertThat(QUADRUPLE.inverse(), is(QUADRUPLE));
    }

    @Test public void aromaticInverse() throws Exception {
        assertThat(AROMATIC.inverse(), is(AROMATIC));
    }

    @Test public void upInverse() throws Exception {
        assertThat(UP.inverse(), is(DOWN));
    }

    @Test public void downInverse() throws Exception {
        assertThat(DOWN.inverse(), is(UP));
    }

    @Test public void implicitInverse() throws Exception {
        assertThat(IMPLICIT.inverse(), is(IMPLICIT));
    }

    @Test public void dotSymbol() throws Exception {
        assertThat(DOT.symbol(), is("."));
    }

    @Test public void singleSymbol() throws Exception {
        assertThat(SINGLE.symbol(), is("-"));
    }

    @Test public void doubleSymbol() throws Exception {
        assertThat(DOUBLE.symbol(), is("="));
    }

    @Test public void tripleSymbol() throws Exception {
        assertThat(TRIPLE.symbol(), is("#"));
    }

    @Test public void quadrupleSymbol() throws Exception {
        assertThat(QUADRUPLE.symbol(), is("$"));
    }

    @Test public void aromaticSymbol() throws Exception {
        assertThat(AROMATIC.symbol(), is(":"));
    }

    @Test public void upSymbol() throws Exception {
        assertThat(UP.symbol(), is("/"));
    }

    @Test public void downSymbol() throws Exception {
        assertThat(DOWN.symbol(), is("\\"));
    }

    @Test public void implicitSymbol() throws Exception {
        assertThat(IMPLICIT.symbol(), is(""));
    }
}
