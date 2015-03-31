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

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

/** @author John May */
public class ParserTest {

    @Test(expected = InvalidSmilesException.class)
    public void ringBondMismatch() throws InvalidSmilesException {
        Parser.decideBond(Bond.SINGLE, Bond.DOUBLE, CharBuffer.fromString(""));
    }

    @Test
    public void ringBondDecision() throws InvalidSmilesException {
        assertThat(Parser.decideBond(Bond.DOUBLE, Bond.DOUBLE, CharBuffer.fromString("")), is(Bond.DOUBLE));
        assertThat(Parser.decideBond(Bond.DOUBLE, Bond.IMPLICIT, CharBuffer.fromString("")), is(Bond.DOUBLE));
        assertThat(Parser.decideBond(Bond.IMPLICIT, Bond.DOUBLE, CharBuffer.fromString("")), is(Bond.DOUBLE));
    }

    @Test public void invalidTetrahedral() throws InvalidSmilesException {
        Graph g = Parser.parse("[C@-](N)(O)C");
        Assert.assertThat(g.topologyOf(0), is(Topology.unknown()));
    }

    @Test public void invalidTetrahedral2() throws InvalidSmilesException {
        Graph g = Parser.parse("[C@](N)(O)C");
        Assert.assertThat(g.topologyOf(0), is(Topology.unknown()));
    }

    @Test(expected = InvalidSmilesException.class)
    public void unclosedRing1() throws Exception {
        Parser.parse("C1CCCCC");
    }

    @Test(expected = InvalidSmilesException.class)
    public void unclosedRing2() throws Exception {
        Parser.parse("C1CCCCC1CCCC1CCCC");
    }

    @Test(expected = InvalidSmilesException.class)
    public void unclosedBranch1() throws Exception {
        Parser.parse("CCCC(CCCC");
    }

    @Test(expected = InvalidSmilesException.class)
    public void unclosedBranch2() throws Exception {
        Parser.parse("CCCC(CCC(CC)");
    }

    @Test(expected = InvalidSmilesException.class)
    public void unopenedBranch1() throws Exception {
        Parser.parse("CCCCCC)CCC");
    }

    @Test(expected = InvalidSmilesException.class)
    public void unopenedBranch2() throws Exception {
        Parser.parse("CCCCCC))CCC");
    }

    @Test public void tellurophene() throws InvalidSmilesException {
        Parser.parse("c1cc[te]c1");
    }

    @Test public void mixingAromaticAndKekule() throws InvalidSmilesException {
        Graph g = Parser.parse("C:1:C:C:C:C:C1");
        for (Edge e : g.edges()) {
            assertThat(e.bond(), is(Bond.AROMATIC));
        }
    }

    @Test public void hydrogen() throws IOException {
        Graph g = Parser.losse("HH");
        assertThat(g.order(), is(2));
        assertThat(g.toSmiles(), is("[H][H]"));
    }

    @Test(expected = InvalidSmilesException.class)
    public void hydrogen_strict() throws IOException {
        Graph g = Parser.strict("HH");
    }

    @Test public void deuterium() throws IOException {
        Graph g = Parser.losse("DD");
        assertThat(g.order(), is(2));
        assertThat(g.toSmiles(), is("[2H][2H]"));
    }

    @Test(expected = InvalidSmilesException.class)
    public void deuterium_strict() throws IOException {
        Graph g = Parser.strict("DD");
    }

    @Test public void tritium() throws IOException {
        Graph g = Parser.losse("TT");
        assertThat(g.order(), is(2));
        assertThat(g.toSmiles(), is("[3H][3H]"));
    }

    @Test(expected = InvalidSmilesException.class)
    public void tritium_strict() throws IOException {
        Graph g = Parser.strict("TT");
    }

    @Test public void tellurium() throws IOException {
        Graph g = Parser.losse("[te]");
        assertTrue(g.atom(0).aromatic());
        assertThat(g.atom(0).element(), is(Element.Tellurium));
    }

    @Test(expected = InvalidSmilesException.class)
    public void tellurium_strict() throws IOException {
        Graph g = Parser.strict("[te]");
    }

    @Test public void largeRnum() throws Exception {
        Graph g = Parser.parse("C%99CCCC%99");
    }

    // not part of spec
    @Test public void r_label() throws InvalidSmilesException {
        Graph g = Parser.parse("CC(C)C[R]");
        assertThat(g.atom(4).label(), is("R"));
    }

    @Test public void random_label() throws InvalidSmilesException {
        Graph g = Parser.parse("CC(C)C[Really?]");
        assertThat(g.atom(4).label(), is("Really?"));
    }
    
    @Test(expected = InvalidSmilesException.class) public void bad_label() throws Exception {
        Parser.parse("[Nope-[not]-[ok]");
    }
    
    @Test(expected = InvalidSmilesException.class) public void bad_label2() throws Exception {
        Parser.parse("[this-[is-not-okay]");
    }

    @Test public void parseCLockwiseExtendedTetrahedral() throws Exception {
        Graph g = Graph.fromSmiles("C(C)=[C@@]=CC");
        assertThat(g.topologyOf(2).configuration(), is(Configuration.AL2));
    }

    @Test public void parseAnticlockwiseExtendedTetrahedral() throws Exception {
        Graph g = Graph.fromSmiles("C(C)=[C@]=CC");
        assertThat(g.topologyOf(2).configuration(), is(Configuration.AL1));
    }
    
    // ek! what a difficult one - this example is from MetaCyc
    @Test public void nested_label() throws InvalidSmilesException {
        Graph g = Parser.parse("CCCCCCC=CCCCCCCCC=CC(=O)[a holo-[acyl-carrier protein]]");
        assertThat(g.atom(g.order() - 1).label(), is("a holo-[acyl-carrier protein]"));
    }
    
    @Test public void seleniumTh() throws Exception {
        assertThat(Parser.parse("[Se@](=O)(C)CC").toSmiles(), is("[Se@](=O)(C)CC"));
    }
    
    @Test public void sulphurIonTh() throws Exception {
        assertThat(Parser.parse("[S@+]([O-])(C)CC").toSmiles(), is("[S@+]([O-])(C)CC"));   
    }

    // chembl has some of these odditites, not sure which tool produced them
    @Test(expected = InvalidSmilesException.class)
    public void rejectChEMBLBadBonds() throws Exception {
        Parser.parse("C\\=C");
    }

    @Test(expected = InvalidSmilesException.class)
    public void rejectMultipleUpBonds() throws Exception {
        Parser.parse("C/C=C(/C)/C");
    }

    @Test
    public void acceptMultipleBonds() throws Exception {
        Parser.parse("C/C=C/C\\C=C/C");
    }
    
    // this one has been mistreated... ignore for now
    @Ignore
    @Test(expected = InvalidSmilesException.class)
    public void chembl345045Mangleded() throws Exception {
        Parser.parse("c1c(ccc(c1)F)c2/c3n/c(c(\\c4[nH]c(/c(c/5\\nc(/c(c/6\\s\\c2\\cc6)/c7ccc(cc7)F)C=C5)/c8ccc(cc8)S(=O)(=O)[O-])cc4)/c9ccc(cc9)S(=O)(=O)[O-])/C=C3.[Na+].[Na+] CHEMBL345045");
    }

    @Test(expected = InvalidSmilesException.class) 
    public void openBracketIsInvalid() throws Exception {
        Parser.parse("[");
    }
}
