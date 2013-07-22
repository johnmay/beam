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

import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertNull;
import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;
import static uk.ac.ebi.grins.Element.Arsenic;
import static uk.ac.ebi.grins.Element.Boron;
import static uk.ac.ebi.grins.Element.Bromine;
import static uk.ac.ebi.grins.Element.Carbon;
import static uk.ac.ebi.grins.Element.Chlorine;
import static uk.ac.ebi.grins.Element.Fluorine;
import static uk.ac.ebi.grins.Element.Iodine;
import static uk.ac.ebi.grins.Element.Nitrogen;
import static uk.ac.ebi.grins.Element.Oxygen;
import static uk.ac.ebi.grins.Element.Phosphorus;
import static uk.ac.ebi.grins.Element.Selenium;
import static uk.ac.ebi.grins.Element.Sulfur;
import static uk.ac.ebi.grins.Element.Unknown;

/** @author John May */
public class ElementTest {

    @Test public void organicSymbols() throws Exception {
        assertThat(Element.ofSymbol("B"), is(Boron));
        assertThat(Element.ofSymbol("C"), is(Carbon));
        assertThat(Element.ofSymbol("N"), is(Nitrogen));
        assertThat(Element.ofSymbol("O"), is(Oxygen));
        assertThat(Element.ofSymbol("P"), is(Phosphorus));
        assertThat(Element.ofSymbol("S"), is(Sulfur));
        assertThat(Element.ofSymbol("F"), is(Fluorine));
        assertThat(Element.ofSymbol("Br"), is(Bromine));
        assertThat(Element.ofSymbol("Cl"), is(Chlorine));
        assertThat(Element.ofSymbol("I"), is(Iodine));
    }

    @Test public void aromaticSymbols() throws Exception {
        assertThat(Element.ofSymbol("b"), is(Boron));
        assertThat(Element.ofSymbol("c"), is(Carbon));
        assertThat(Element.ofSymbol("n"), is(Nitrogen));
        assertThat(Element.ofSymbol("o"), is(Oxygen));
        assertThat(Element.ofSymbol("p"), is(Phosphorus));
        assertThat(Element.ofSymbol("s"), is(Sulfur));
        assertThat(Element.ofSymbol("se"), is(Selenium));
        assertThat(Element.ofSymbol("as"), is(Arsenic));
    }

    @Test public void symbols() {
        for (Element e : Element.values()) {
            assertThat(Element.ofSymbol(e.symbol()), is(e));
        }
    }

    @Test public void invalidSymbol() {
        assertNull(Element.ofSymbol("J"));
    }

    @Test public void organic() {
        for (Element e : Arrays.asList(Boron,
                                       Carbon,
                                       Nitrogen,
                                       Oxygen,
                                       Phosphorus,
                                       Sulfur,
                                       Fluorine,
                                       Chlorine,
                                       Bromine,
                                       Iodine)) {
            assertTrue(e.organic());
        }
    }

    @Test public void aromatic() {
        for (Element e : Arrays.asList(Boron,
                                       Carbon,
                                       Nitrogen,
                                       Oxygen,
                                       Phosphorus,
                                       Sulfur,
                                       Selenium,
                                       Arsenic)) {
            assertTrue(e.aromatic());
        }
    }

    @Test(expected = IllegalArgumentException.class)
    public void inorganicHydrogens() {
        Unknown.implicitHydrogens(0);
    }

    @Test
    public void boronHydrogens() {
        assertThat(Boron.implicitHydrogens(0), is(3));
        assertThat(Boron.implicitHydrogens(1), is(2));
        assertThat(Boron.implicitHydrogens(2), is(1));
        assertThat(Boron.implicitHydrogens(3), is(0));
        assertThat(Boron.implicitHydrogens(4), is(0));
    }

    @Test
    public void carbonHydrogens() {
        assertThat(Carbon.implicitHydrogens(0), is(4));
        assertThat(Carbon.implicitHydrogens(1), is(3));
        assertThat(Carbon.implicitHydrogens(2), is(2));
        assertThat(Carbon.implicitHydrogens(3), is(1));
        assertThat(Carbon.implicitHydrogens(4), is(0));
        assertThat(Carbon.implicitHydrogens(5), is(0));
        assertThat(Carbon.implicitHydrogens(6), is(0));
    }

    @Test
    public void nitrogenHydrogens() {
        assertThat(Nitrogen.implicitHydrogens(0), is(3));
        assertThat(Nitrogen.implicitHydrogens(1), is(2));
        assertThat(Nitrogen.implicitHydrogens(2), is(1));
        assertThat(Nitrogen.implicitHydrogens(3), is(0));
        assertThat(Nitrogen.implicitHydrogens(4), is(1));
        assertThat(Nitrogen.implicitHydrogens(5), is(0));
        assertThat(Nitrogen.implicitHydrogens(6), is(0));
    }

    @Test
    public void oxygenHydrogens() {
        assertThat(Oxygen.implicitHydrogens(0), is(2));
        assertThat(Oxygen.implicitHydrogens(1), is(1));
        assertThat(Oxygen.implicitHydrogens(2), is(0));
        assertThat(Oxygen.implicitHydrogens(3), is(0));
    }

    @Test
    public void phosphorusHydrogens() {
        assertThat(Phosphorus.implicitHydrogens(0), is(3));
        assertThat(Phosphorus.implicitHydrogens(1), is(2));
        assertThat(Phosphorus.implicitHydrogens(2), is(1));
        assertThat(Phosphorus.implicitHydrogens(3), is(0));
        assertThat(Phosphorus.implicitHydrogens(4), is(1));
        assertThat(Phosphorus.implicitHydrogens(5), is(0));
        assertThat(Phosphorus.implicitHydrogens(6), is(0));
    }

    @Test
    public void sulfurHydrogens() {
        assertThat(Sulfur.implicitHydrogens(0), is(2));
        assertThat(Sulfur.implicitHydrogens(1), is(1));
        assertThat(Sulfur.implicitHydrogens(2), is(0));
        assertThat(Sulfur.implicitHydrogens(3), is(1));
        assertThat(Sulfur.implicitHydrogens(4), is(0));
        assertThat(Sulfur.implicitHydrogens(5), is(1));
        assertThat(Sulfur.implicitHydrogens(6), is(0));
        assertThat(Sulfur.implicitHydrogens(7), is(0));
    }

    @Test
    public void halogenHydrogens() {
        for (Element e : Arrays.asList(Fluorine, Chlorine, Bromine, Iodine)) {
            assertThat(e.implicitHydrogens(0), is(1));
            assertThat(e.implicitHydrogens(1), is(0));
            assertThat(e.implicitHydrogens(2), is(0));
        }
    }
}
