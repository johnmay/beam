package uk.ac.ebi.beam;

import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;

public class RelaxedParsingTest {

    public static void assertRelaxed(String exp, String smi) throws IOException {
        Assert.assertThrows(smi + " should throw and exception when not relaxed parsing",
                            InvalidSmilesException.class, () -> {Parser.parse(smi);});
        String act = Parser.relaxed(smi).toSmiles();
        Assert.assertEquals(smi + " was not parsed expected", exp, act);
    }

    @Test
    public void testRelaxedParsing_unclosedRings() throws IOException {
        assertRelaxed("CCCC", "C1CCC");
        assertRelaxed("CCCC", "C%CCC");
        assertRelaxed("C.CCC", "C.1CCC");
        assertRelaxed("CC1CCCC1", "CC1CCCC12");
        assertRelaxed("CC1CCCC1", "CC12CCCC1");
        assertRelaxed("CC1CCCC1", "CC12CCCC1");
        assertRelaxed("CC", "C1C1");
        assertRelaxed("C#1CCC1", "C=1CCC#1");
        assertRelaxed("C=1CCC1", "C#1CCC=1");
    }

    @Test
    public void testRelaxedParsing_unclosedBranches() throws IOException {
        assertRelaxed("CCCC", "C(CCC");
        assertRelaxed("CCCC", "C((CCC");
        assertRelaxed("CCCC", "C)CCC");
        assertRelaxed("CCCC", "C))CCC");
        assertRelaxed("CCC(CO)O", "CCC(CO)O(");
    }

    @Test
    public void testRelaxedParsing_bondTypes() throws IOException {
        assertRelaxed("CC#CC", "CC-=#CC");
        assertRelaxed("CC=CC", "CC-==CC");
        assertRelaxed("CC=CC", "CC#=CC");
        assertRelaxed("CC.CC", "CC..CC");
        assertRelaxed("CC.CC", "CC..CC");
        // assertRelaxed("CC", ".CC");
        assertRelaxed("CC", "=CC");
    }

    @Test
    public void testRelaxedParsing_atoms() throws IOException {
        assertRelaxed("CC[NH2]", "CC[NH2");
        assertRelaxed("CC[*]", "CC[NBoc");
        assertRelaxed("CC*", "CC[");
        assertRelaxed("CC", "]CC");
        assertRelaxed("CC", "@]CC");
        assertRelaxed("[H]CC", "@H]CC");
        assertRelaxed("CC", "+]CC");
    }

}
