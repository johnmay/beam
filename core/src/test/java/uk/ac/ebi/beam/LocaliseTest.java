package uk.ac.ebi.beam;

import org.junit.Test;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

/** @author John May */
public class LocaliseTest {

    @Test public void furan() throws Exception {
        test("o1cccc1", "O1C=CC=C1");
    }

    @Test public void benzen() throws Exception {
        test("c1ccccc1", "C1=CC=CC=C1");
    }

    @Test public void quinone() throws Exception {
        test("oc1ccc(o)cc1", "O=C1C=CC(=O)C=C1");
        test("O=c1ccc(=O)cc1", "O=C1C=CC(=O)C=C1");
    }

    @Test(expected = InvalidSmilesException.class)
    public void methane() throws Exception {
        test("c", "C"); // note daylight makes it 'CH3' but we say - valence error
    }

    @Test public void ethene() throws Exception {
        test("cc", "C=C");
    }

    @Test(expected = InvalidSmilesException.class)
    public void invalid_acyclic_chain() throws Exception {
        test("ccc", "n/a");
    }

    @Test public void buta_1_3_diene() throws Exception {
        test("cccc", "C=CC=C");
    }

    // some allow lower-case to be radical, this should throw an exception
    @Test(expected = InvalidSmilesException.class)
    public void carbon_radical() throws Exception {
        test("C1CCcCC1", "n/a");
    }

    @Test public void _hexa_1_3_5_triene() throws Exception {
        test("cccccc", "C=CC=CC=C");
    }

    @Test public void _4H_pyran_4_one() throws Exception {
        test("oc1ccocc1", "O=C1C=COC=C1");
    }

    @Test public void pyrole() throws Exception {
        test("[nH]1cccc1", "[NH]1C=CC=C1");
    }

    @Test public void imidazole() throws Exception {
        test("c1c[nH]cn1", "C1=C[NH]C=N1");
    }

    @Test public void benzimidazole() throws Exception {
        test("c1nc2ccccc2[nH]1", "C1=NC2=CC=CC=C2[NH]1");
    }

    @Test public void napthalene() throws Exception {
        test("c1ccc2ccccc2c1", "C1=CC=C2C=CC=CC2=C1");
    }

    @Test public void anthracene() throws Exception {
        test("c1ccc2cc3ccccc3cc2c1", "C1=CC=C2C=C3C=CC=CC3=CC2=C1");
    }

    @Test public void thiophene() throws Exception {
        test("s1cccc1", "S1C=CC=C1");
    }

    @Test public void imidazol_3_ium() throws Exception {
        test("c1c[nH+]c[nH]1", "C1=C[NH+]=C[NH]1");
    }

    @Test public void exocyclic_NO_bond() throws Exception {
        test("Nc1c2nc[nH]c2ncn1=O", "NC1=C2N=C[NH]C2=NC=N1=O");
    }

    @Test public void biphenyl() throws Exception {
        test("c1ccccc1c1ccccc1", "C1=CC=CC=C1C2=CC=CC=C2");
        test("c1ccccc1-c1ccccc1", "C1=CC=CC=C1C2=CC=CC=C2");
    }

    @Test public void phospho_nitro_ring() throws Exception {
        test("n1pnpnp1", "N1=PN=PN=P1");
    }

    @Test public void phospho_nitro_ring_exocyclic_oxygen() throws Exception {
        test("n1p(O)(O)np(O)(O)np1(O)(O)", "N1=P(O)(O)N=P(O)(O)N=P1(O)O");
    }

    @Test public void hexamethylidenecyclohexane() throws Exception {
        test("cc1c(c)c(c)c(c)c(c)c1c", "C=C1C(=C)C(=C)C(=C)C(=C)C1=C");
        test("C=c1c(=C)c(=C)c(=C)c(=C)c1=C", "C=C1C(=C)C(=C)C(=C)C(=C)C1=C");
    }

    @Test(expected = InvalidSmilesException.class)
    public void tryptophanyl_radical() throws Exception {
        test("NC(Cc1c[n]c2ccccc12)C(O)=O",
             "n/a");
    }

    @Test public void thiophene_oxide() throws Exception {
        test("O=s1cccc1",
             "O=S1C=CC=C1");
    }

    @Test public void tellurophene() throws Exception {
        test("[Te]1cccc1", "[Te]1C=CC=C1");
        test("[te]1cccc1", "[Te]1C=CC=C1");
    }

    // Sulphur with two double bonds
    @Test public void chembl_1188068() throws Exception {
        test("COc1cc2nc(ns(=O)(C)c2cc1OC)N3CCN(CC3)C(=O)c4oc(SC)nn4",
             "COC1=CC2=NC(=NS(=O)(C)=C2C=C1OC)N3CCN(CC3)C(=O)C=4OC(SC)=NN4");
    }

    // Sulphur cation with exo cyclic double bond
    @Test public void chembl_1434989() throws Exception {
        test("[O-][s+]1(=O)[nH]c2c(cc(Cl)c3ccc(Cl)nc23)c4ccccc14",
             "[O-][S+]1(=O)[NH]C2=C(C=C(Cl)C3=CC=C(Cl)N=C23)C4=CC=CC=C14");
    }

    @Test public void chembl_423544() throws Exception {
        test("CCc1n[c]#[c]n1CC2CC(C(=O)O2)(c3ccccc3)c4ccccc4",
             "CCC1=N[C]#[C]N1CC2CC(C(=O)O2)(C3=CC=CC=C3)C4=CC=CC=C4");
    }

    @Test public void tropylium() throws Exception {
        test("[cH+]1cccccc1", "[CH+]1C=CC=CC=C1");
    }

    @Test(expected = InvalidSmilesException.class)
    public void pyrole_invalid() throws Exception {
        test("n1cncc1", "n/a");
    }

    @Test(expected = InvalidSmilesException.class)
    public void imidazole_invalid() throws Exception {
        test("c1nc2ccccc2n1", "n/a");
    }

    @Test
    public void mixing_aromatic_and_aliphatic() throws Exception {
        test("c1=cc=cc=c1", "C1=CC=CC=C1");
        test("c-1c-cc-cc1", "C1=CC=CC=C1");
        test("C:1:C:C:C:C:C1", "C1CCCCC1");
    }

    // http://sourceforge.net/mailarchive/forum.php?thread_name=60825b0f0709302037g2d68f2eamdb5ebecf3baea6d1%40mail.gmail.com&forum_name=blueobelisk-smiles
    @Test public void bezene_inconsistent() throws Exception {
        test("c1=ccccc1", "C1=CC=CC=C1");
        test("c1=cc=ccc1", "C1=CC=CC=C1");
        test("c1=cc=cc=c1", "C1=CC=CC=C1");
        test("c1=c:c:c:c:c1", "C1=CC=CC=C1");
        test("c1=c:c=c:c:c1", "C1=CC=CC=C1");
        test("c1=c-c=c:c:c1", "C1=CC=CC=C1");
    }

    @Test public void fluorene() throws Exception {
        test("C1c2ccccc2-c3ccccc13", "C1C2=CC=CC=C2C3=CC=CC=C13");
        test("C1c2ccccc2c3ccccc13", "C1C2=CC=CC=C2C3=CC=CC=C13");
    }

    @Test public void hexaoxane() throws Exception {
        test("o1ooooo1", "O1OOOOO1");
    }

    @Test public void pyrole_aliphatic_n() throws Exception {
        test("c1cNcc1", "C1=CNC=C1");
    }

    @Test public void furan_aliphatic_o() throws Exception {
        test("c1cOcc1", "C1=COC=C1");
    }

    @Test public void bo_25756() throws Exception {
        test("Nc1c2c3ccccc3c4cccc(cc1)c24",
             "NC1=C2C3=CC=CC=C3C4=CC=CC(C=C1)=C24");
    }
    
    /* Examples from http://www.daylight.com/dayhtml_tutorials/languages/smiles/smiles_examples.html */

    @Test public void viagra() throws Exception {
        test("CCc1nn(C)c2c(=O)[nH]c(nc12)c3cc(ccc3OCC)S(=O)(=O)N4CCN(C)CC4",
             "CCC1=NN(C)C=2C(=O)[NH]C(=NC12)C3=CC(=CC=C3OCC)S(=O)(=O)N4CCN(C)CC4");
    }

    @Test public void xanax() throws Exception {
        test("Cc1nnc2CN=C(c3ccccc3)c4cc(Cl)ccc4-n12",
             "CC1=NN=C2CN=C(C3=CC=CC=C3)C4=CC(Cl)=CC=C4N12");
    }

    @Test public void phentermine() throws Exception {
        test("CC(C)(N)Cc1ccccc1",
             "CC(C)(N)CC1=CC=CC=C1");
    }

    @Test public void valium() throws Exception {
        test("CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13",
             "CN1C(=O)CN=C(C2=CC=CC=C2)C3=CC(Cl)=CC=C13");
    }

    @Test public void ambien() throws Exception {
        test("CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c3ccc(C)cc3",
             "CN(C)C(=O)CC1=C(N=C2C=CC(C)=CN12)C3=CC=C(C)C=C3");
    }

    @Test public void nexium() throws Exception {
        test("COc1ccc2[nH]c(nc2c1)S(=O)Cc3ncc(C)c(OC)c3C",
             "COC1=CC=C2[NH]C(=NC2=C1)S(=O)CC3=NC=C(C)C(OC)=C3C");
    }

    @Test public void vioxx() throws Exception {
        test("CS(=O)(=O)c1ccc(cc1)C2=C(C(=O)OC2)c3ccccc3",
             "CS(=O)(=O)C1=CC=C(C=C1)C2=C(C(=O)OC2)C3=CC=CC=C3");
    }

    @Test public void paxil() throws Exception {
        test("Fc1ccc(cc1)C2CCNCC2COc3ccc4OCOc4c3",
             "FC1=CC=C(C=C1)C2CCNCC2COC3=CC=C4OCOC4=C3");
    }

    @Test public void lipitor() throws Exception {
        test("CC(C)c1c(C(=O)Nc2ccccc2)c(c(c3ccc(F)cc3)n1CC[C@@H]4C[C@@H](O)CC(=O)O4)c5ccccc5",
             "CC(C)C1=C(C(=O)NC2=CC=CC=C2)C(=C(C3=CC=C(F)C=C3)N1CC[C@@H]4C[C@@H](O)CC(=O)O4)C5=CC=CC=C5");
    }

    @Test public void cialis() throws Exception {
        test("CN1CC(=O)N2[C@@H](c3[nH]c4ccccc4c3C[C@@H]2C1=O)c5ccc6OCOc6c5",
             "CN1CC(=O)N2[C@@H](C=3[NH]C4=CC=CC=C4C3C[C@@H]2C1=O)C5=CC=C6OCOC6=C5");
    }

    @Test public void strychnine() throws Exception {
        test("O=C1C[C@H]2OCC=C3CN4CC[C@@]56[C@H]4C[C@H]3[C@H]2[C@H]6N1c7ccccc75",
             "O=C1C[C@H]2OCC=C3CN4CC[C@]56[C@H]4C[C@H]3[C@H]2[C@H]5N1C7=CC=CC=C76");
    }

    @Test public void cocaine() throws Exception {
        test("COC(=O)[C@H]1[C@@H]2CC[C@H](C[C@@H]1OC(=O)c3ccccc3)N2C",
             "COC(=O)[C@H]1[C@@H]2CC[C@H](C[C@@H]1OC(=O)C3=CC=CC=C3)N2C");
    }

    @Test public void quinine() throws Exception {
        test("COc1ccc2nccc([C@@H](O)[C@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
             "COC1=CC=C2N=CC=C([C@@H](O)[C@H]3C[C@@H]4CCN3C[C@@H]4C=C)C2=C1");
    }

    @Test public void lysergicAcid() throws Exception {
        test("CN1C[C@@H](C=C2[C@H]1Cc3c[nH]c4cccc2c34)C(=O)O",
             "CN1C[C@@H](C=C2[C@H]1CC3=C[NH]C4=CC=CC2=C34)C(=O)O");
    }

    @Test public void LSD() throws Exception {
        test("CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3c[nH]c4cccc(C2=C1)c34",
             "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2CC3=C[NH]C4=CC=CC(C2=C1)=C34");
    }

    @Test public void morphine() throws Exception {
        test("CN1CC[C@]23[C@H]4Oc5c3c(C[C@@H]1[C@@H]2C=C[C@@H]4O)ccc5O",
             "CN1CC[C@@]23[C@H]4OC5=C2C(C[C@@H]1[C@@H]3C=C[C@@H]4O)=CC=C5O");
    }

    @Test public void heroin() throws Exception {
        test("CN1CC[C@]23[C@H]4Oc5c3c(C[C@@H]1[C@@H]2C=C[C@@H]4OC(=O)C)ccc5OC(=O)C",
             "CN1CC[C@@]23[C@H]4OC5=C2C(C[C@@H]1[C@@H]3C=C[C@@H]4OC(=O)C)=CC=C5OC(=O)C");
    }

    @Test public void nicotine() throws Exception {
        test("CN1CCC[C@H]1c2cccnc2",
             "CN1CCC[C@H]1C2=CC=CN=C2");
    }

    @Test public void caffeine() throws Exception {
        test("Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
             "CN1C=NC=2N(C)C(=O)N(C)C(=O)C12");
    }

    // N,N-Diallylmelamine 
    @Test public void ncs4420() throws Exception {
        test("[nH2]c1nc(nc(n1)n(Ccc)Ccc)[nH2]",
             "[NH2]C1=NC(=NC(=N1)N(CC=C)CC=C)[NH2]");
    }

    @Test public void carbon_anion() throws Exception {
        test("O=c1cc[cH-]cc1",
             "O=C1C=C[CH-]C=C1");
        test("oc1cc[cH-]cc1",
             "O=C1C=C[CH-]C=C1");
    }

    static void test(String delocalised, String localised) throws Exception {
        Graph g = Graph.fromSmiles(delocalised);
        Graph h = Localise.localise(g);
        assertThat(h.toSmiles(), is(localised));
    }

}
