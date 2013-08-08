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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Generate a SMILES line notation for a given chemical graph.
 *
 * @author John May
 */
final class Generator {

    private final ChemicalGraph g;
    private final StringBuilder sb;

    private final int[] visitedAt;
    private int i, rnum = 1;
    private final AtomToken[]                     tokens;
    private final Map<Integer, List<RingClosure>> rings;

    /**
     * Create a new generator the given chemical graph.
     *
     * @param g chemical graph
     */
    Generator(ChemicalGraph g) {
        this.g = g;
        this.sb = new StringBuilder(g.order() * 2);
        this.visitedAt = new int[g.order()];
        this.tokens = new AtomToken[g.order()];
        this.rings = new HashMap<Integer, List<RingClosure>>();

        // prepare ring closures and topologies
        int visited = 0;
        Arrays.fill(visitedAt, -1);
        for (int u = 0; u < g.order() && visited < g.order(); u++) {
            if (visitedAt[u] < 0)
                prepare(u, u, u == 0 ? Bond.IMPLICIT : Bond.DOT);
        }

        // write notation
        visited = 0;
        Arrays.fill(visitedAt, -1);
        for (int u = 0; u < g.order() && visited < g.order(); u++) {
            if (visitedAt[u] < 0)
                write(u, u, u == 0 ? Bond.IMPLICIT : Bond.DOT);
        }
    }

    /**
     * First traversal of the molecule assigns ring bonds (numbered later) and
     * configures topologies.
     *
     * @param u the vertex to visit
     * @param p parent vertex (where we came from)
     * @param b bond type
     */
    void prepare(int u, int p, Bond b) {
        visitedAt[u] = i++;
        tokens[u] = g.atom(u).token();

        for (Edge e : g.edges(u)) {
            int v = e.other(u);
            if (visitedAt[v] < 0) {
                prepare(v, u, e.bond(u));
            } else if (v != p && visitedAt[v] < visitedAt[u]) {
                cyclicEdge(v, u, e.bond(v));
            }
        }

        // Configure the topology using the traversal order
        // XXX: need a better approach for extended tetrahedral (@AL1/@AL2)
        if (rings.containsKey(u)) {
            tokens[u].configure(g.topologyOf(u)
                                 .orderBy(localRank(u, p))
                                 .configuration());
        } else {
            tokens[u].configure(g.topologyOf(u)
                                 .orderBy(visitedAt)
                                 .configuration());
        }
    }

    /**
     * Second traversal writes the bonds and atoms to the SMILES string.
     *
     * @param u a vertex
     * @param p previous vertex
     * @param b the bond from the previous vertex to this vertex
     */
    void write(int u, int p, Bond b) {
        visitedAt[u] = i++;

        int remaining = g.degree(u);

        if (u != p)
            remaining--;

        // assign ring numbers
        if (rings.containsKey(u)) {
            for (RingClosure rc : rings.get(u)) {
                if (rc.register(rnum)) {
                    int v = rc.other(u);
                    tokens[u] = new RingNumberToken(new RingBondToken(tokens[u],
                                                                      rc.bond(u)),
                                                    rnum);
                    tokens[v] = new RingNumberToken(tokens[v],
                                                    rnum);
                    rnum++;
                }
                remaining--;
            }
        }

        sb.append(b.symbol());
        tokens[u].append(sb);

        for (Edge e : g.edges(u)) {
            int v = e.other(u);
            if (visitedAt[v] < 0) {
                if (--remaining > 0) {
                    sb.append('(');
                    write(v, u, e.bond(u));
                    sb.append(')');
                } else {
                    write(v, u, e.bond(u));
                }
            }
        }
    }

    /**
     * Local ordering around a vertex, required for stereo-configurations which
     * have ring closures. Here we need to consider the order of the ring bonds
     * in the stereo specification.
     *
     * @param u a vertex
     * @param p the vertex we came from (parent)
     * @return the local rank for the neighbors of the vertex u
     */
    private int[] localRank(int u, int p) {
        int[] localRank = new int[g.order()];
        int rank = 2;
        localRank[u] = 1;
        for (Edge e : g.edges(u)) {
            localRank[e.other(u)] = g.degree(u) + rank++;
        }
        localRank[p] = 0;
        rank = 2;
        for (RingClosure rc : rings.get(u)) {
            localRank[rc.other(u)] = rank++;
        }
        return localRank;
    }

    /**
     * Indicate that the edge connecting the vertices u and v forms a ring.
     *
     * @param u a vertex
     * @param v a vertex connected to u
     * @param b bond type connecting u to v
     */
    private void cyclicEdge(int u, int v, Bond b) {
        RingClosure r = new RingClosure(u, v, b);
        addRing(r.u, r);
        addRing(r.v, r);
    }

    /**
     * Add a ring closure to the the vertex 'u'.
     *
     * @param u  a vertex
     * @param rc ring closure
     */
    private void addRing(int u, RingClosure rc) {
        List<RingClosure> closures = rings.get(u);
        if (closures == null) {
            closures = new ArrayList<RingClosure>(2);
            rings.put(u, closures);
        }
        closures.add(rc);
    }

    /**
     * Access the generated SMILES string.
     *
     * @return smiles string
     */
    String string() {
        return sb.toString();
    }

    /**
     * Convenience method for generating a SMILES string for the specified
     * chemical graph.
     *
     * @param g the graph to generate the SMILE for
     * @return SMILES gor the provided chemical graph
     */
    static String generate(final ChemicalGraph g) {
        return new Generator(g).string();
    }

    public static void main(String[] args) throws InvalidSmilesException {
        System.out.println(Generator
                                   .generate(Parser.parse("C1/CCCCCCCCCCC\\C=1")));
        System.out.println(Generator
                                   .generate(Parser.parse("C=1/CCCCCCCCCCC\\C=1")));
        System.out.println(Generator
                                   .generate(Parser.parse("C=1/CCCCCCCCCCC\\C1")));
        System.out.println(Generator
                                   .generate(Parser.parse("C\\1CCCCCCCCCCC\\C=C1")));
        System.out.println(Generator
                                   .generate(Parser.parse("C1CCCCCCCCCCC\\C=C/1")));
        System.out.println(Generator
                                   .generate(Parser.parse("CCCC.OOOO.C[CH]C.CNO")));
    }

    static final class RingClosure {
        final int u, v;
        final Bond b;
        int rnum = -1;

        RingClosure(int u, int v, Bond b) {
            this.u = u;
            this.v = v;
            this.b = b;
        }

        int other(int x) {
            if (x == u) return v;
            if (x == v) return u;
            throw new IllegalArgumentException("non edge endpoint");
        }

        Bond bond(int x) {
            if (x == u) return b;
            else if (x == v) return b.inverse();
            throw new IllegalArgumentException("invalid endpoint");
        }

        boolean register(int rnum) {
            if (this.rnum < 0) {
                this.rnum = rnum;
                return true;
            }
            return false;
        }
    }

    static interface AtomToken {
        void configure(Configuration c);

        void append(StringBuilder sb);
    }

    static final class SubsetToken implements AtomToken {
        private final String str;

        SubsetToken(String str) {
            this.str = str;
        }

        @Override public void configure(Configuration c) {
            // do nothing
        }

        @Override public void append(StringBuilder sb) {
            sb.append(str);
        }
    }

    static final class BracketToken implements AtomToken {

        private Atom atom;
        private Configuration c = Configuration.UNKNOWN;

        BracketToken(Atom a) {
            this.atom = a;
        }

        @Override public void configure(Configuration c) {
            this.c = c;
        }

        @Override public void append(StringBuilder sb) {
            sb.append('[');
            if (atom.isotope() >= 0)
                sb.append(atom.isotope());
            sb.append(atom.aromatic() ? atom.element()
                                            .symbol()
                                            .toLowerCase(Locale.ENGLISH)
                                      : atom.element()
                                            .symbol());
            if (c != Configuration.UNKNOWN)
                sb.append(c.shorthand().symbol());
            if (atom.hydrogens() > 0)
                sb.append(Element.Hydrogen.symbol());
            if (atom.hydrogens() > 1)
                sb.append(atom.hydrogens());
            if (atom.charge() != 0) {
                sb.append(atom.charge() > 0 ? '+' : '-');
                int absCharge = Math.abs(atom.charge());
                if (absCharge > 1)
                    sb.append(absCharge);
            }
            if (atom.atomClass() != 0)
                sb.append(':').append(atom.atomClass());
            sb.append(']');
        }
    }

    static abstract class TokenAdapter implements AtomToken {

        private AtomToken parent;

        TokenAdapter(AtomToken parent) {
            this.parent = parent;
        }

        @Override public final void configure(Configuration c) {
            this.parent.configure(c);
        }

        @Override public void append(StringBuilder sb) {
            parent.append(sb);
        }
    }

    static final class RingNumberToken extends TokenAdapter {
        int rnum;

        RingNumberToken(AtomToken p, int rnum) {
            super(p);
            this.rnum = rnum;
        }

        @Override public void append(StringBuilder sb) {
            super.append(sb);
            if (rnum > 9)
                sb.append('%');
            sb.append(rnum);
        }
    }

    static final class RingBondToken extends TokenAdapter {
        Bond bond;

        RingBondToken(AtomToken p, Bond bond) {
            super(p);
            this.bond = bond;
        }

        @Override public void append(StringBuilder sb) {
            super.append(sb);
            sb.append(bond);
        }
    }
}
