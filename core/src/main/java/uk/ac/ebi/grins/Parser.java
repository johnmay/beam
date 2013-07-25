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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import static uk.ac.ebi.grins.Atom.AromaticSubset;
import static uk.ac.ebi.grins.Atom.OrganicSubset;


/** @author John May */
final class Parser {

    private final IntStack stack = new IntStack(10);
    private final ChemicalGraph g;
    private RingBond[] rings = new RingBond[10];

    private Map<Integer, LocalArrangement> arrangement
            = new HashMap<Integer, LocalArrangement>(5);

    Parser(CharBuffer buffer) throws InvalidSmilesException {
        g = new ChemicalGraph(1 + (2 * (buffer.length() / 3)));
        readSmiles(buffer);
    }

    Parser(String str) throws InvalidSmilesException {
        this(CharBuffer.fromString(str));
    }

    static ChemicalGraph parse(String str) throws InvalidSmilesException {
        return new Parser(str).molecule();
    }

    ChemicalGraph molecule() {
        return g;
    }

    private void addAtom(Atom a) {
        int v = g.addAtom(a);
        if (!stack.empty()) {
            int u = stack.pop();
            if (bond != Bond.DOT)
                g.addEdge(new Edge(u, v, bond));
            bond = Bond.IMPLICIT;
        }
        stack.push(v);
    }

    // primary dispatch table
    private void readSmiles(CharBuffer buffer) throws InvalidSmilesException {
        while (buffer.hasRemaining()) {
            char c = buffer.get();
            switch (c) {

                // organic subset
                case 'B':
                    if (buffer.getIf('r'))
                        addAtom(OrganicSubset.Bromine);
                    else
                        addAtom(OrganicSubset.Boron);
                    break;
                case 'C':
                    if (buffer.getIf('l'))
                        addAtom(OrganicSubset.Chlorine);
                    else
                        addAtom(OrganicSubset.Carbon);
                    break;
                case 'N':
                    addAtom(OrganicSubset.Nitrogen);
                    break;
                case 'O':
                    addAtom(OrganicSubset.Oxygen);
                    break;
                case 'P':
                    addAtom(OrganicSubset.Phosphorus);
                    break;
                case 'S':
                    addAtom(OrganicSubset.Sulfur);
                    break;
                case 'F':
                    addAtom(OrganicSubset.Fluorine);
                    break;
                case 'I':
                    addAtom(OrganicSubset.Iodine);
                    break;

                // aromatic subset
                case 'b':
                    addAtom(AromaticSubset.Boron);
                    break;
                case 'c':
                    addAtom(AromaticSubset.Carbon);
                    break;
                case 'n':
                    addAtom(AromaticSubset.Nitrogen);
                    break;
                case 'o':
                    addAtom(AromaticSubset.Oxygen);
                    break;
                case 'p':
                    addAtom(AromaticSubset.Phosphorus);
                    break;
                case 's':
                    addAtom(AromaticSubset.Sulfur);
                    break;

                // bracket atom
                case '[':
                    addAtom(readBracketAtom(buffer));
                    break;

                // ring bonds
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                    ring(c - '0');
                    break;
                case '%':
                    int num = buffer.getNumber();
                    if (num < 0)
                        throw new InvalidSmilesException("number (<digit>+) must follow '%'", buffer);
                    ring(num);
                    break;

                // bond/dot
                case '-':
                    bond = Bond.SINGLE;
                    break;
                case '=':
                    bond = Bond.DOUBLE;
                    break;
                case '#':
                    bond = Bond.TRIPLE;
                    break;
                case '$':
                    bond = Bond.QUADRUPLE;
                    break;
                case ':':
                    bond = Bond.AROMATIC;
                    break;
                case '/':
                    bond = Bond.UP;
                    break;
                case '\\':
                    bond = Bond.DOWN;
                    break;
                case '.':
                    bond = Bond.DOT;
                    break;

                // branching
                case '(':
                    stack.push(stack.peek());
                    break;
                case ')':
                    stack.pop();
                    break;

                // termination
                case '\t':
                case ' ':
                case '\n':
                case '\r':
                    return;

                default:
                    throw new InvalidSmilesException("unexpected character", buffer);
            }
        }
    }

    static Atom readBracketAtom(CharBuffer buffer) throws
                                                    InvalidSmilesException {

        // try to read isotope number, -1 if not read
        int isotope = buffer.getNumber();

        // lowercase indicates aromatic
        boolean aromatic = buffer.next() >= 'a' && buffer.next() <= 'z';

        Element element = Element.read(buffer);

        Configuration configuration = Configuration.read(buffer);

        int hCount = readHydrogens(buffer);
        int charge = readCharge(buffer);
        int atomClass = readClass(buffer);

        if (!buffer.getIf(']'))
            throw InvalidSmilesException.invalidBracketAtom(buffer);

        // create the atom
        return new Atom.BracketAtom(isotope < 0 ? 0 : isotope,
                                    element,
                                    hCount,
                                    charge,
                                    atomClass,
                                    aromatic);
    }

    /**
     * Read the hydrogen count and progress the provided buffer. The hydrogen
     * count is specified by a 'H' an 0 or more digits. A 'H' without digits is
     * intercepted as 'H1'. When there is no 'H' or 'H0' is specified then the
     * the hydrogen count is 0.
     *
     * @param buffer a character buffer
     * @return the hydrogen count, 0 if none
     */
    static int readHydrogens(final CharBuffer buffer) {
        if (buffer.getIf('H')) {
            // when no number is specified 'H' then there is 1 hydrogen
            int count = buffer.getNumber();
            return count < 0 ? 1 : count;
        }
        return 0;
    }

    /**
     * Read a charge value and progress the provide buffer. The charge value is
     * present in bracket atoms either directly after the symbol, the chiral
     * specification or the hydrogen count. The specification of charge by
     * concatenated signs (e.g. ++, --) and other bad form (e.g. '++-1') is
     * intercepted.
     *
     * @param buffer a character buffer
     * @return the formal charge value, 0 if none present
     * @see <a href="http://www.opensmiles.org/opensmiles.html#charge">Charge -
     *      OpenSMILES Specification</a>
     */
    static int readCharge(final CharBuffer buffer) {
        return readCharge(0, buffer);
    }

    /**
     * Internal method for parsing charge, to allow concatenated signs (--, ++)
     * the method recursively invokes increment or decrementing an accumulator.
     *
     * @param acc    accumulator
     * @param buffer a character buffer
     * @return the charge value
     */
    private static int readCharge(int acc, final CharBuffer buffer) {
        if (buffer.getIf('+'))
            return buffer.nextIsDigit() ? acc + buffer.getNumber()
                                        : readCharge(acc + 1, buffer);
        if (buffer.getIf('-'))
            return buffer.nextIsDigit() ? acc - buffer.getNumber()
                                        : readCharge(acc - 1, buffer);
        return acc;
    }

    /**
     * Read the atom class of a bracket atom and progress the buffer (if read).
     * The atom class is the last attribute of the bracket atom and is
     * identified by a ':' followed by one or more digits. The atom class may be
     * padded such that ':005' and ':5' are equivalent.
     *
     * @param buffer a character buffer
     * @return the atom class, or 0
     * @see <a href="http://www.opensmiles.org/opensmiles.html#atomclass">Atom
     *      Class - OpenSMILES Specification</a>
     */
    static int readClass(CharBuffer buffer) throws InvalidSmilesException {
        if (buffer.getIf(':')) {
            if (buffer.nextIsDigit())
                return buffer.getNumber();
            throw new InvalidSmilesException("invalid atom class, <digit>+ must follow ':'", buffer);
        }
        return 0;
    }

    private void ring(int rnum) throws InvalidSmilesException {
        if (rings.length <= rnum || rings[rnum] == null) {
            openRing(rnum);
        } else {
            closeRing(rnum);
        }
    }

    private void openRing(int rnum) {
        if (rnum >= rings.length)
            rings = Arrays.copyOf(rings, rnum + 1);
        int u = stack.peek();

        // create a ring bond
        rings[rnum] = new RingBond(u, bond);

        // keep track of arrangement (important for stereo configurations)
        createArrangement(u).add(-rnum);

        bond = Bond.IMPLICIT;
    }

    private LocalArrangement createArrangement(int u) {
        LocalArrangement la = arrangement.get(u);
        if (la == null) {
            la = new LocalArrangement();
            for (Edge e : g.edges(stack.peek()))
                la.add(e.other(u));
            arrangement.put(u, la);
        }
        return la;
    }

    private void closeRing(int rnum) throws InvalidSmilesException {
        RingBond rbond = rings[rnum];
        rings[rnum] = null;
        g.addEdge(new Edge(rbond.u, stack.peek(),
                           decideBond(rbond.bond, bond.inverse())));
        bond = Bond.IMPLICIT;
        // adjust the arrangement replacing where this ring number was openned
        arrangement.get(rbond.u).replace(-rnum, stack.peek());
    }

    /**
     * Decide the bond to use for a ring bond. The bond symbol can be present on
     * either or both bonded atoms. This method takes those bonds, chooses the
     * correct one or reports an error if there is a conflict.
     *
     * Equivalent SMILES:
     * <blockquote><pre>
     *     C=1CCCCC=1
     *     C=1CCCCC1    (preferred)
     *     C1CCCCC=1
     * </pre></blockquote>
     *
     * @param a a bond
     * @param b other bond
     * @return the bond to use for this edge
     * @throws InvalidSmilesException ring bonds did not match
     */
    static Bond decideBond(final Bond a, final Bond b) throws
                                                       InvalidSmilesException {
        if (a == b)
            return a;
        else if (a == Bond.IMPLICIT)
            return b;
        else if (b == Bond.IMPLICIT)
            return a;
        throw new InvalidSmilesException("ring bond mismatch, " + a + " and " + b);
    }

    public static void main(String[] args) throws IOException {

        String path = "/databases/zinc/zinc.smi";
        BufferedReader r = new BufferedReader(new FileReader(path));
        String line = null;
        List<String> smis = new ArrayList<String>();
        while ((line = r.readLine()) != null)
            smis.add(line);
        r.close();

        for (int i = 0; i < 10; i++) {
            int c = 0;
            long t0 = System.nanoTime();
            for (String smi : smis) {
                new Parser(smi);
            }
            long t1 = System.nanoTime();
            System.out.println(TimeUnit.NANOSECONDS.toMillis(t1 - t0) + " ms");
        }
    }

    private static class RingBond {
        int  u;
        Bond bond;

        private RingBond(int u, Bond bond) {
            this.u = u;
            this.bond = bond;
        }
    }

    private static final class LocalArrangement {

        int[] vs;
        int   n;

        private LocalArrangement() {
            this.vs = new int[4];
        }

        void add(final int v) {
            if (n == vs.length)
                vs = Arrays.copyOf(vs, n * 2);
            vs[n++] = v;
        }

        void replace(final int u, final int v) {
            for (int i = 0; i < n; i++) {
                if (vs[i] == u) {
                    vs[i] = v;
                    return;
                }
            }
        }
    }
}
