package uk.ac.ebi.grins;

import org.hamcrest.CoreMatchers;
import org.junit.Assert;
import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.sameInstance;
import static org.junit.Assert.assertThat;

/** @author John May */
public class FromSubsetAtomsTest {

    @Test public void unknown() throws Exception {
        transform("*", "[*]");
    }

    @Test public void inorganic() throws Exception {
        transform("[Ne]", "[Ne]");
    }

    @Test public void methane() throws Exception {
        transform("C", "[CH4]");
    }

    @Test public void ethane_withAtomClass() throws Exception {
        transform("[CH3:1]C", "[CH3:1][CH3]");
    }

    @Test public void ethanol() throws InvalidSmilesException {
        transform("CCO", "[CH3][CH2][OH]");
    }

    @Test public void stereoSpecification() throws InvalidSmilesException {
        transform("[C@H](N)(O)C", "[C@H]([NH2])([OH])[CH3]");
        transform("[C@@H](N)(O)C", "[C@@H]([NH2])([OH])[CH3]");
    }

    @Test public void noStereoSpecification() throws InvalidSmilesException {
        transform("C(N)(O)C", "[CH]([NH2])([OH])[CH3]");
    }

    @Test public void bracketAtom() {
        // should provide identity of bracket atom
        Atom input  = new AtomImpl.BracketAtom(Element.Carbon, 1, 0);
        Atom output = FromSubsetAtoms.fromSubset(input, 0);
        Assert.assertThat(input, CoreMatchers.is(CoreMatchers
                                                         .sameInstance(output)));
    }

    @Test public void aliphatic_carbon() {
        Atom actual = FromSubsetAtoms.fromSubset(AtomImpl.AliphaticSubset.Carbon, 3);
        Atom expect = new AtomImpl.BracketAtom(Element.Carbon, 1, 0);
        Assert.assertThat(expect, CoreMatchers.is(actual));
    }

    @Test public void aromatic_carbon() {
        Atom actual = FromSubsetAtoms.fromSubset(AtomImpl.AromaticSubset.Carbon, 3);
        Atom expect = new AtomImpl.BracketAtom(-1, Element.Carbon, 1, 0, 0, true);
        Assert.assertThat(expect, CoreMatchers.is(actual));
    }

    private void transform(String input, String expected) throws
                                                          InvalidSmilesException {
        ChemicalGraph g = Parser.parse(input);
        ImplicitToExplicit ite = new ImplicitToExplicit();
        FromSubsetAtoms    fsa = new FromSubsetAtoms();
        ExplicitToImplicit eti = new ExplicitToImplicit();
        String actual = Generator.generate(eti.transform(
                                           fsa.transform(
                                           ite.transform(g))));
        Assert.assertThat(actual, CoreMatchers.is(expected));
    }

}
