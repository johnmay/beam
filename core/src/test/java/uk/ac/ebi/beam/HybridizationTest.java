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

package uk.ac.ebi.beam;

import org.junit.Test;

import java.io.IOException;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;
import static uk.ac.ebi.beam.Hybridization.Sp;
import static uk.ac.ebi.beam.Hybridization.Sp2;
import static uk.ac.ebi.beam.Hybridization.Sp3;
import static uk.ac.ebi.beam.Hybridization.hybridization;

/** @author John May */
public class HybridizationTest {

    @Test public void methane() throws IOException {
        Graph g = Graph.fromSmiles("[CH4]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp3));
    }

    @Test public void ethene() throws IOException {
        Graph g = Graph.fromSmiles("[CH2]=[CH2]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp2));
        assertThat(Hybridization.hybridization(g, 1), is(Sp2));
    }

    @Test public void ethane() throws IOException {
        Graph g = Graph.fromSmiles("[CH3][CH3]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp3));
        assertThat(Hybridization.hybridization(g, 1), is(Sp3));
    }

    @Test public void carbonDioxide() throws IOException {
        Graph g = Graph.fromSmiles("[O]=[C]=[O]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp2));
        assertThat(Hybridization.hybridization(g, 1), is(Sp));
        assertThat(Hybridization.hybridization(g, 2), is(Sp2));
    }

    @Test public void carbonCation() throws IOException {
        Graph g = Graph.fromSmiles("[CH3+]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp2));
    }

    @Test public void carbonAnion() throws IOException {
        Graph g = Graph.fromSmiles("[CH3-]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp3));
    }

    @Test public void carbonRadical() throws IOException {
        Graph g = Graph.fromSmiles("[CH3]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp2));
    }

    @Test public void ammonia() throws IOException {
        Graph g = Graph.fromSmiles("[NH3]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp3));
    }

    @Test public void nitrogenSp3() throws IOException {
        Graph g = Graph.fromSmiles("[CH3][NH][CH3]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp3));
        assertThat(Hybridization.hybridization(g, 1), is(Sp3));
        assertThat(Hybridization.hybridization(g, 2), is(Sp3));
    }

    @Test public void nitrogenSp2() throws IOException {
        Graph g = Graph.fromSmiles("[CH3][N]=[CH2]");
        assertThat(Hybridization.hybridization(g, 0), is(Sp3));
        assertThat(Hybridization.hybridization(g, 1), is(Sp2));
        assertThat(Hybridization.hybridization(g, 2), is(Sp2));
    }

    @Test public void nitrogenSp3_exception1() throws IOException {
        Graph g = Graph.fromSmiles("[CH2]=[CH2][NH][CH2]=[CH2]");

        // hybridizations are Sp2, Sp2, Sp3 and Sp2, Sp2
        assertThat(hybridization(g, 0), is(Sp2));
        assertThat(hybridization(g, 1), is(Sp2));
        assertThat(hybridization(g, 2), is(Sp3));
        assertThat(hybridization(g, 3), is(Sp2));
        assertThat(hybridization(g, 4), is(Sp2));

        // but as the nitrogen is next to an Sp2 it is also Sp2
        Hybridization[] hs = Hybridization.hybridizations(g);
        assertThat(hs, is(new Hybridization[]{Sp2, Sp2, Sp2, Sp2, Sp2}));
    }

    @Test public void nitrogenSp3_exception2() throws IOException {
        Graph g = Graph.fromSmiles("[CH2]=[CH2][NH2]");

        assertThat(hybridization(g, 0), is(Sp2));
        assertThat(hybridization(g, 1), is(Sp2));
        assertThat(hybridization(g, 2), is(Sp3));

        Hybridization[] hs = Hybridization.hybridizations(g);
        assertThat(hs, is(new Hybridization[]{Sp2, Sp2, Sp2}));
    }

    @Test public void oxygenSp3_exception() throws IOException {
        Graph g = Graph.fromSmiles("[CH3][CH3][O][CH2]=[CH2]");

        assertThat(hybridization(g, 0), is(Sp3));
        assertThat(hybridization(g, 1), is(Sp3));
        assertThat(hybridization(g, 2), is(Sp3));
        assertThat(hybridization(g, 3), is(Sp2));
        assertThat(hybridization(g, 4), is(Sp2));

        Hybridization[] hs = Hybridization.hybridizations(g);
        assertThat(hs, is(new Hybridization[]{Sp3, Sp3, Sp2, Sp2, Sp2}));
    }

    @Test public void carbonSp3_exception() throws IOException {
        Graph g = Graph.fromSmiles("[CH2]=[CH2][CH2-]");

        assertThat(hybridization(g, 0), is(Sp2));
        assertThat(hybridization(g, 1), is(Sp2));
        assertThat(hybridization(g, 2), is(Sp3));

        Hybridization[] hs = Hybridization.hybridizations(g);
        assertThat(hs, is(new Hybridization[]{Sp2, Sp2, Sp2}));
    }

    // crazy examples - shows the Sp2 propagates around the ring
    @Test public void oxygeny() throws Exception {
        Graph g = Graph.fromSmiles("[O]=[CH]1[O][O][O][CH](=[O])[CH]1=[O]");
        Hybridization[] hs = Hybridization.hybridizations(g);
        assertThat(hs, is(new Hybridization[]{Sp2, Sp2, Sp2,
                                              Sp2, Sp2, Sp2,
                                              Sp2, Sp2, Sp2}));
    }

    @Test public void methylfuranium() throws Exception {
        Graph g = Graph.fromSmiles("[CH3][O+]1[CH]=[CH][CH]=[CH]1");
        Hybridization[] hs = Hybridization.hybridizations(g);
        assertThat(hs, is(new Hybridization[]{Sp3, Sp2, Sp2,
                                              Sp2, Sp2, Sp2}));
    }

    @Test public void oxygenLonePairs() throws Exception {
        Graph g = Graph.fromSmiles("[OH2]");
        assertThat(Hybridization.lonePairs(g, 0), is(2));
    }

    @Test public void oxygenAnionLonePairs() throws Exception {
        Graph g = Graph.fromSmiles("[OH-]");
        assertThat(Hybridization.lonePairs(g, 0), is(3));
    }

    @Test public void oxygenDiAnionLonePairs() throws Exception {
        Graph g = Graph.fromSmiles("[O-2]");
        assertThat(Hybridization.lonePairs(g, 0), is(4));
    }
}
