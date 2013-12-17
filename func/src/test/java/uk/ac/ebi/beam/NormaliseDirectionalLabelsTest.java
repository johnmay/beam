package uk.ac.ebi.beam;

import org.hamcrest.CoreMatchers;
import org.junit.Assert;
import org.junit.Test;

/** @author John May */
public class NormaliseDirectionalLabelsTest {

    @Test public void simple() throws InvalidSmilesException {
        transform("F\\C=C\\F",
                  "F/C=C/F");
    }

    @Test public void ordering2() throws InvalidSmilesException {
        transform("C(\\F)(/C)=C\\F",
                  "C(/F)(\\C)=C/F");
    }

    @Test public void simple2() throws InvalidSmilesException {
        transform("C(\\F)=C\\F",
                  "C(/F)=C/F");
    }

    @Test public void partial() throws InvalidSmilesException {
        transform("FC=C(\\F)/C=C/F",
                  "FC=C(\\F)/C=C/F");
    }

    @Test public void partial2() throws InvalidSmilesException {
        transform("FC=C(F)C=C(F)\\C=C\\F",
                  "FC=C(F)C=C(F)/C=C/F");
    }

    @Test public void conjugated() throws InvalidSmilesException {
        transform("F\\C=C(\\F)/C(/F)=C\\F",
                  "F/C=C(/F)\\C(\\F)=C/F");
    }
    
    @Test public void cyclic() throws InvalidSmilesException {
        transform("C/C=C\\1/C\\C(=C\\C)\\C1",
                  "C/C=C\\1/C/C(=C/C)/C1");
        transform("C/C=C\\1/C/C(=C/C)/C1",
                  "C/C=C\\1/C/C(=C/C)/C1");
    }
    
    @Test public void chebi15617() throws InvalidSmilesException {
        transform("C/C=C\\1/[C@@H](C)C(=O)N/C1=C\\C/2=N/C(=C\\C3=C(CCC(=O)O)C(=C(\\C=C/4\\C(=C(CC)C(=O)N4)C)N3)C)/C(=C2C)CCC(=O)O",
                  "C/C=C\\1/[C@@H](C)C(=O)N/C1=C\\C/2=N/C(=C\\C3=C(CCC(=O)O)C(=C(/C=C\\4/C(=C(CC)C(=O)N4)C)N3)C)/C(=C2C)CCC(=O)O");
        transform("C/C=C\\1/[C@@H](C)C(=O)N/C1=C\\C/2=N/C(=C\\C3=C(CCC(=O)O)C(=C(/C=C\\4/C(=C(CC)C(=O)N4)C)N3)C)/C(=C2C)CCC(=O)O",
                  "C/C=C\\1/[C@@H](C)C(=O)N/C1=C\\C/2=N/C(=C\\C3=C(CCC(=O)O)C(=C(/C=C\\4/C(=C(CC)C(=O)N4)C)N3)C)/C(=C2C)CCC(=O)O");
    }
    

    static void transform(String smi, String exp) throws
                                                  InvalidSmilesException {
        Assert.assertThat(Generator.generate(new NormaliseDirectionalLabels()
                                                     .apply(Parser.parse(smi))), CoreMatchers
                                  .is(exp));
    }

}
