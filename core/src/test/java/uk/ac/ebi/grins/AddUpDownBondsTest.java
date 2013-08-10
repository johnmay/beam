/*
 * Copyright (c) 2013, European Bioinformatics Institute (EMBL-EBI)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the FreeBSD Project.
 */

package uk.ac.ebi.grins;

import org.junit.Assert;
import org.junit.Test;

import static org.hamcrest.CoreMatchers.is;

/** @author John May */
public class AddUpDownBondsTest {

    @Test public void e_butene_implH() throws Exception {
        transform("C\\C=C\\C",
                  "C\\C=C\\C");
        transform("C/C=C/C",
                  "C/C=C/C");
        transform("C\\C([H])=C\\C",
                  "C\\C(\\[H])=C\\C");
        transform("C/C([H])=C/C",
                  "C/C(/[H])=C/C");
        transform("C\\C([H])=C([H])\\C",
                  "C\\C(\\[H])=C(/[H])\\C");
        transform("C/C([H])=C([H])/C",
                  "C/C(/[H])=C(\\[H])/C");
    }

    @Test public void z_butene_implH() throws Exception {
        transform("C\\C=C/C",
                  "C\\C=C/C");
        transform("C/C=C\\C",
                  "C/C=C\\C");
        transform("C\\C([H])=C/C",
                  "C\\C(\\[H])=C/C");
        transform("C/C([H])=C\\C",
                  "C/C(/[H])=C\\C");
        transform("C\\C([H])=C([H])/C",
                  "C\\C(\\[H])=C(\\[H])/C");
        transform("C/C([H])=C(/[H])\\C",
                  "C/C(/[H])=C(/[H])\\C");
    }

    @Test public void e_butene_expH() throws Exception {
        transform("C\\C([H])=C\\C",
                  "C\\C(\\[H])=C\\C");
        transform("C/C([H])=C/C",
                  "C/C(/[H])=C/C");
        transform("C\\C([H])=C([H])\\C",
                  "C\\C(\\[H])=C(/[H])\\C");
        transform("C/C([H])=C([H])/C",
                  "C/C(/[H])=C(\\[H])/C");
    }

    @Test public void z_butene_expH() throws Exception {
        transform("C\\C([H])=C/C",
                  "C\\C(\\[H])=C/C");
        transform("C/C([H])=C\\C",
                  "C/C(/[H])=C\\C");
        transform("C\\C([H])=C([H])/C",
                  "C\\C(\\[H])=C(\\[H])/C");
        transform("C/C([H])=C(/[H])\\C",
                  "C/C(/[H])=C(/[H])\\C");
    }

    @Test(expected = InvalidSmilesException.class)
    public void invalid() throws Exception {
        transform("C/C=C(/C)/([H])",
                  "n/a");
    }

    @Test
    public void conjugated2() throws Exception {
        transform("C/C=C(/C)C(\\[H])=C/C",
                  "C/C=C(/C)\\C(\\[H])=C/C");
    }

    @Test public void e_e_hexadiene_implH() throws InvalidSmilesException {
        transform("C\\C=C\\C=C\\C",
                  "C\\C=C\\C=C\\C");
        transform("C/C=C/C=C/C",
                  "C/C=C/C=C/C");
    }

    @Test public void z_z_hexadiene_implH() throws InvalidSmilesException {
        transform("C\\C=C/C=C\\C",
                  "C\\C=C/C=C\\C");
        transform("C/C=C\\C=C/C",
                  "C/C=C\\C=C/C");
    }

    @Test public void e_z_hexadiene_implH() throws InvalidSmilesException {
        transform("C\\C=C\\C=C/C",
                  "C\\C=C\\C=C/C");
        transform("C/C=C/C=C\\C",
                  "C/C=C/C=C\\C");
    }

    @Test public void z_e_hexadiene_implH() throws InvalidSmilesException {
        transform("C/C=C\\C=C\\C",
                  "C/C=C\\C=C\\C");
        transform("C\\C=C/C=C/C",
                  "C\\C=C/C=C/C");
    }

    @Test public void e_e_hexadiene_expH() throws InvalidSmilesException {
        transform("C\\C([H])=C([H])\\C([H])=C([H])\\C",
                  "C\\C(\\[H])=C(/[H])\\C(\\[H])=C(/[H])\\C");
        transform("C/C([H])=C([H])/C([H])=C([H])/C",
                  "C/C(/[H])=C(\\[H])/C(/[H])=C(\\[H])/C");
    }

    @Test public void z_z_hexadiene_expH() throws InvalidSmilesException {
        transform("C\\C([H])=C([H])/C([H])=C([H])\\C",
                  "C\\C(\\[H])=C(\\[H])/C(/[H])=C(/[H])\\C");
        transform("C/C([H])=C([H])\\C([H])=C([H])/C",
                  "C/C(/[H])=C(/[H])\\C(\\[H])=C(\\[H])/C");
    }

    @Test public void e_z_hexadiene_expH() throws InvalidSmilesException {
        transform("C\\C([H])=C([H])\\C([H])=C([H])/C",
                  "C\\C(\\[H])=C(/[H])\\C(\\[H])=C(\\[H])/C");
        transform("C/C([H])=C([H])/C([H])=C([H])\\C",
                  "C/C(/[H])=C(\\[H])/C(/[H])=C(/[H])\\C");
    }

    @Test public void z_e_hexadiene_expH() throws InvalidSmilesException {
        transform("C/C(/[H])=C(/[H])\\C(\\[H])=C(/[H])\\C",
                  "C/C(/[H])=C(/[H])\\C(\\[H])=C(/[H])\\C");
        transform("C\\C(\\[H])=C(\\[H])/C(/[H])=C(\\[H])/C",
                  "C\\C(\\[H])=C(\\[H])/C(/[H])=C(\\[H])/C");
    }

    /** Ensure cumulated double bonds don't break the transformation */
    @Test public void allene() throws InvalidSmilesException {
        transform("CC=C=CC",
                  "CC=C=CC");
    }

    /** Ensure cumulated double bonds don't break the transformation */
    @Test public void cumullene() throws InvalidSmilesException {
        transform("CC=C=C=CC",
                  "CC=C=C=CC");
    }

    @Test(expected = InvalidSmilesException.class)
    public void invalidConjugated() throws Exception {
        transform("C/C=C(/C)C(/[H])=C/C",
                  "n/a");
    }

    /**
     * One double bond has an end point with no directional bonds. Fully
     * specifying the other double bond adds one directional bond giving it a
     * readable configuration.
     */
    @Test public void partial1() throws InvalidSmilesException {
        transform("[H]\\C(C([H])=C(\\C)[H])=C(/C)[H]",
                  "[H]\\C(\\C([H])=C(\\C)/[H])=C(/C)\\[H]");
    }


    /** The directional label between two double bonds is missing. */
    @Test public void middle() throws InvalidSmilesException {
        transform("C(\\[H])(=C(\\C)[H])C(=C(/C)[H])/[H]",
                  "C(\\[H])(=C(\\C)/[H])/C(=C(/C)\\[H])/[H]");
    }

    static void transform(String smi, String exp) throws
                                                  InvalidSmilesException {
        Assert.assertThat(Generator.generate(new AddUpDownBonds()
                                                     .transform(Parser.parse(smi))), is(exp));
    }
}
