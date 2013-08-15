package uk.ac.ebi.grins;

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
                  "FC=C(/F)\\C=C\\F");
    }

    @Test public void partial2() throws InvalidSmilesException {
        transform("FC=C(F)C=C(F)\\C=C\\F",
                  "FC=C(F)C=C(F)/C=C/F");
    }

    @Test public void conjugated() throws InvalidSmilesException {
        transform("F\\C=C(\\F)/C(/F)=C\\F",
                  "F/C=C(/F)\\C(\\F)=C/F");
    }

    static void transform(String smi, String exp) throws
                                                  InvalidSmilesException {
        Assert.assertThat(Generator.generate(new NormaliseDirectionalLabels()
                                                     .apply(Parser.parse(smi))), CoreMatchers
                                  .is(exp));
    }

}
