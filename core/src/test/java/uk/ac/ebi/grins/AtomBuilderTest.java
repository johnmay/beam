package uk.ac.ebi.grins;

import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

/** @author John May */
public class AtomBuilderTest {

    @Test public void create_element_c() throws Exception {
        Atom a = AtomBuilder.create(Element.Carbon)
                            .build();
        assertThat(a.element(), is(Element.Carbon));
        assertThat(a.isotope(), is(-1));
        assertThat(a.charge(), is(0));
        assertThat(a.atomClass(), is(0));
        assertThat(a.aromatic(), is(false));
    }

    @Test public void create_element_n() throws Exception {
        Atom a = AtomBuilder.create(Element.Nitrogen)
                            .build();
        assertThat(a.element(), is(Element.Nitrogen));
        assertThat(a.isotope(), is(-1));
        assertThat(a.charge(), is(0));
        assertThat(a.atomClass(), is(0));
        assertThat(a.aromatic(), is(false));
    }

    @Test(expected = NullPointerException.class)
    public void create_element_null() throws Exception {
        Atom a = AtomBuilder.create((Element) null)
                            .build();
    }

    @Test public void create_symbol_aliphatic_c() throws Exception {
        Atom a = AtomBuilder.create("C")
                            .build();
        assertThat(a.element(), is(Element.Carbon));
        assertThat(a.isotope(), is(-1));
        assertThat(a.charge(), is(0));
        assertThat(a.atomClass(), is(0));
        assertThat(a.aromatic(), is(false));
    }

    @Test public void create_symbol_aromatic_c() throws Exception {
        Atom a = AtomBuilder.create("c")
                            .build();
        assertThat(a.element(), is(Element.Carbon));
        assertThat(a.isotope(), is(-1));
        assertThat(a.charge(), is(0));
        assertThat(a.atomClass(), is(0));
        assertThat(a.aromatic(), is(true));
    }

    @Test(expected = IllegalArgumentException.class)
    public void create_symbol_non_aromatic() throws Exception {
        Atom a = AtomBuilder.create("cl")
                            .build();
    }

    @Test
    public void create_symbol_defaultToUnknown() throws Exception {
        Atom a = AtomBuilder.create("N/A")
                            .build();
        assertThat(a.element(), is(Element.Unknown));
    }

    @Test
    public void create_symbol_null() throws Exception {
        Atom a = AtomBuilder.create((String) null)
                            .build();
        assertThat(a.element(), is(Element.Unknown));
    }
}
