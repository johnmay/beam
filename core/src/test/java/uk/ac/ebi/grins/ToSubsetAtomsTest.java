package uk.ac.ebi.grins;

import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.sameInstance;
import static org.junit.Assert.assertThat;

/** @author John May */
public class ToSubsetAtomsTest {

    @Test public void inorganic() throws Exception {
        transform("[*]", "[*]");
    }

    @Test public void inorganic2() throws Exception {
        transform("[Ne]", "[Ne]");
    }

    @Test public void methane() throws Exception {
        transform("[CH4]", "C");
    }

    @Test public void monovalent_carbon() throws Exception {
        transform("[CH3]", "[CH3]");
    }

    @Test public void divalent_carbon() throws Exception {
        transform("[CH2]", "[CH2]");
    }

    @Test public void trivalent_carbon() throws Exception {
        transform("[CH]", "[CH]");
        transform("[CH1]", "[CH]");
    }

    @Test public void carbon_12() throws Exception {
        // note the isotope is specified and so must be a bracket atom
        transform("[12C]", "[12C]");
    }

    @Test public void carbon_13() throws Exception {
        transform("[13C]", "[13C]");
    }

    @Test public void carbon_14() throws Exception {
        transform("[14C]", "[14C]");
    }

    @Test public void oxidanide() throws Exception {
        transform("[OH-]", "[OH-]");
    }

    @Test public void azanium() throws Exception {
        transform("[NH4+]", "[NH4+]");
    }

    @Test public void ethane_withAtomClass() throws Exception {
        transform("[CH3:1][CH3:0]", "[CH3:1]C");
    }

    @Test public void ethanol() throws InvalidSmilesException {
        transform("[CH3][CH2][OH]", "CCO");
    }

    @Test public void stereoSpecification() throws InvalidSmilesException {
        transform("[C@H]([NH2])([OH])[CH3]", "[C@H](N)(O)C");
        transform("[C@@H]([NH2])([OH])[CH3]", "[C@@H](N)(O)C");
    }

    @Test public void noStereoSpecification() throws InvalidSmilesException {
        transform("[CH]([NH2])([OH])[CH3]", "C(N)(O)C");
    }

    @Test public void aliphaticSubset() throws Exception {
        for (Atom a : Atom.AliphaticSubset.values()){
            assertThat(toSubset(a, 0), is(sameInstance(a)));
        }
    }

    @Test public void aromaticSubset() throws Exception {
        for (Atom a : Atom.AromaticSubset.values()){
            assertThat(toSubset(a, 0), is(sameInstance(a)));
        }
    }

    @Test public void atomWithCharge() throws Exception {
        Atom a = new Atom.BracketAtom(-1, Element.Oxygen, 0, -1, 0, false);
        assertThat(toSubset(a, 0), is(sameInstance(a)));
    }

    @Test public void atomWithClassLabel() throws Exception {
        Atom a = new Atom.BracketAtom(-1, Element.Oxygen, 1, 0, 1, false);
        assertThat(toSubset(a, 0), is(sameInstance(a)));
    }

    @Test public void atomWithIsotopeLabel() throws Exception {
        Atom a = new Atom.BracketAtom(0, Element.Sulfur, 0, 0, 0, false);
        assertThat(toSubset(a, 0), is(sameInstance(a)));
    }

    @Test public void atomWithRequiredHydrogens() throws Exception {
        Atom a = new Atom.BracketAtom(0, Element.Oxygen, 2, 0, 0, false);
        assertThat(toSubset(a, 0), is(sameInstance(a)));
    }

    private Atom toSubset(Atom a, int bondOrderSum) {
        return new ToSubsetAtoms().toSubset(a, bondOrderSum);
    }

    private void transform(String input, String expected) throws
                                                          InvalidSmilesException {
        ChemicalGraph g = Parser.parse(input);
        ImplicitToExplicit ite = new ImplicitToExplicit();
        ToSubsetAtoms tsa = new ToSubsetAtoms();
        ExplicitToImplicit eti = new ExplicitToImplicit();
        String actual = Generator.generate(eti.transform(
                                           tsa.transform(
                                           ite.transform(g))));
        assertThat(actual, is(expected));
    }

}
