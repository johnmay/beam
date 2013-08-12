package uk.ac.ebi.grins;

import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

/** @author John May */
public class AtomTest {

    @Test public void aliphaticSubsetFromElement() {
        for (Atom a : Atom.AliphaticSubset.values()) {
            assertThat(Atom.AliphaticSubset.ofElement(a.element()), is(a));
        }
    }

    @Test(expected = IllegalArgumentException.class)
    public void aliphaticSubsetInvalidElement() {
        Atom.AliphaticSubset.ofElement(Element.Californium);
    }

    @Test public void aromaticSubsetFromElement() {
        for (Atom a : Atom.AromaticSubset.values()) {
            assertThat(Atom.AromaticSubset.ofElement(a.element()), is(a));
        }
    }

    @Test(expected = IllegalArgumentException.class)
    public void aromaticSubsetInvalidElement() {
        Atom.AromaticSubset.ofElement(Element.Unknown);
    }
}
