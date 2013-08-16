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
public class GraphBuilderTest {

    @Test
    public void clockwise_parity() {

        GraphBuilder gb = GraphBuilder.create(5);
        ChemicalGraph g = gb.add(AtomBuilder.aliphatic("C").build())
                            .add(AtomImpl.AliphaticSubset.Nitrogen)
                            .add(AtomImpl.AliphaticSubset.Oxygen)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomBuilder.explicitHydrogen())
                            .add(0, 1)
                            .add(0, 2)
                            .add(0, 3)
                            .add(0, 4)
                            .tetrahedral(0).lookingFrom(1)
                            .neighbors(2, 3, 4)
                            .parity(1)
                            .build()
                            .build();

        Assert.assertThat(g.toSmiles(), is("[C@@](N)(O)(C)[H]"));
    }

    @Test
    public void anticlockwise_parity() {

        GraphBuilder gb = GraphBuilder.create(5);
        ChemicalGraph g = gb.add(AtomBuilder.aliphatic("C").build())
                            .add(AtomImpl.AliphaticSubset.Nitrogen)
                            .add(AtomImpl.AliphaticSubset.Oxygen)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomBuilder.explicitHydrogen())
                            .add(0, 1)
                            .add(0, 2)
                            .add(0, 3)
                            .add(0, 4)
                            .tetrahedral(0).lookingFrom(1)
                            .neighbors(2, 3, 4)
                            .parity(-1)
                            .build()
                            .build();

        Assert.assertThat(g.toSmiles(), is("[C@](N)(O)(C)[H]"));
    }

    @Test
    public void e_1_2_difluroethene() {
        GraphBuilder gb = GraphBuilder.create(5);
        ChemicalGraph g = gb.add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(0, 1)
                            .doubleBond(1, 2)
                            .add(2, 3)
                            .geometric(1, 2).opposite(0, 3)
                            .build();
        Assert.assertThat(g.toSmiles(), is("F/C=C/F"));
    }

    @Test
    public void z_1_2_difluroethene() {
        GraphBuilder gb = GraphBuilder.create(5);
        ChemicalGraph g = gb.add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(0, 1)
                            .doubleBond(1, 2)
                            .add(2, 3)
                            .geometric(1, 2).together(0, 3)
                            .build();
        Assert.assertThat(g.toSmiles(), is("F/C=C\\F"));
    }


    @Test
    public void conjugated_consider_existing() {
        // the second configuration considers the existing configuration
        GraphBuilder gb = GraphBuilder.create(5);
        ChemicalGraph g = gb.add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(0, 1)
                            .doubleBond(1, 2)
                            .add(2, 3)
                            .doubleBond(3, 4)
                            .add(4, 5)
                            .geometric(1, 2).together(0, 3)
                            .geometric(3, 4).together(2, 5)
                            .build();
        Assert.assertThat(g.toSmiles(), is("F/C=C\\C=C/F"));
    }

    @Test
    public void conjugated_resolve_conflict() {
        // assigning the second one first means we have to consider this
        // on the first one
        GraphBuilder gb = GraphBuilder.create(5);
        ChemicalGraph g = gb.add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(0, 1)
                            .doubleBond(1, 2)
                            .add(2, 3)
                            .doubleBond(3, 4)
                            .add(4, 5)
                            .geometric(3, 4).together(2, 5)
                            .geometric(1, 2).together(0, 3)
                            .build();
        Assert.assertThat(g.toSmiles(), is("F\\C=C/C=C\\F"));
    }

    @Test
    public void conjugated_resolve_conflict2() {
        // we assign the first, third then second - the second one cause
        // a conflict and we must flip one of the others
        GraphBuilder gb = GraphBuilder.create(5);
        ChemicalGraph g = gb.add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Fluorine)
                            .add(0, 1)
                            .doubleBond(1, 2)
                            .add(2, 3)
                            .doubleBond(3, 4)
                            .add(4, 5)
                            .doubleBond(5, 6)
                            .add(6, 7)
                            .geometric(1, 2).opposite(0, 3)
                            .geometric(5, 6).together(4, 7)
                            .geometric(3, 4).together(2, 5)
                            .build();
        Assert.assertThat(g.toSmiles(), is("F/C=C/C=C\\C=C/F"));
    }

    @Test
    public void all_trans_octatetraene() {
        GraphBuilder gb = GraphBuilder.create(5);
        ChemicalGraph g = gb.add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(0, 1)
                            .doubleBond(1, 2)
                            .add(2, 3)
                            .doubleBond(3, 4)
                            .add(4, 5)
                            .doubleBond(5, 6)
                            .add(6, 7)
                            .doubleBond(7, 0)
                            .geometric(1, 2).together(0, 3)
                            .geometric(3, 4).together(2, 5)
                            .geometric(5, 6).together(4, 7)
                            .geometric(7, 0).together(6, 1)
                            .build();
        Assert.assertThat(g.toSmiles(), is("C=1/C=C\\C=C/C=C\\C1"));
    }

    @Test(expected = IllegalArgumentException.class)
    public void impossible_octatetraene() {
        GraphBuilder gb = GraphBuilder.create(5);
        ChemicalGraph g = gb.add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(AtomImpl.AliphaticSubset.Carbon)
                            .add(0, 1)
                            .doubleBond(1, 2)
                            .add(2, 3)
                            .doubleBond(3, 4)
                            .add(4, 5)
                            .doubleBond(5, 6)
                            .add(6, 7)
                            .doubleBond(7, 0)
                            .geometric(1, 2).together(0, 3)
                            .geometric(3, 4).opposite(2, 5)
                            .geometric(5, 6).together(4, 7)
                            .geometric(7, 0).together(6, 1)
                            .build();
        Assert.assertThat(g.toSmiles(), is("C=1/C=C\\C=C/C=C\\C1"));
    }

}
