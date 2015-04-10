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

import java.io.IOException;
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

    private final Graph         g;
    private final StringBuilder sb;

    private final int[]                           visitedAt;
    private       int                             i;
    private final AtomToken[]                     tokens;
    private final Map<Integer, List<RingClosure>> rings;
    private final RingNumbering                   rnums;

    /**
     * Create a new generator the given chemical graph.
     *
     * @param g chemical graph
     */
    Generator(Graph g, RingNumbering rnums) throws InvalidSmilesException {
        this(g, new int[g.order()], rnums);    
    }
    
    /**
     * Create a new generator the given chemical graph.
     *
     * @param g chemical graph
     * @param visitedAt the index of the atom in the output         
     */
    Generator(Graph g, int[] visitedAt, RingNumbering rnums) throws InvalidSmilesException {
        this.g = g;
        this.rnums = rnums;
        this.sb = new StringBuilder(g.order() * 2);
        this.visitedAt = visitedAt;
        this.tokens = new AtomToken[g.order()];
        this.rings = new HashMap<Integer, List<RingClosure>>();

        // prepare ring closures and topologies
        int visited = 0;
        Arrays.fill(visitedAt, -1);
        boolean uncofiguredStereo = false; 
        for (int u = 0; u < g.order() && visited < g.order(); u++) {
            if (visitedAt[u] < 0)
                uncofiguredStereo = prepare(u, u, u == 0 ? Bond.IMPLICIT : Bond.DOT) || uncofiguredStereo;
        }

        if (uncofiguredStereo) {
            for (int u = 0; u < g.order(); u++) {
                if (g.topologyOf(u).configuration().type() == Configuration.Type.ExtendedTetrahedral) {
                    tokens[u].configure(g.topologyOf(u)
                                         .orderBy(visitedAt)
                                         .configuration());
                }
            }
        }

        // write notation
        visited = 0;
        i = 0;
        Arrays.fill(visitedAt, -1);
        for (int u = 0; u < g.order() && visited < g.order(); u++) {
            if (visitedAt[u] < 0) {
                if (u > 0) {
                    rnums.reset();
                    write(u, u, Bond.DOT);
                } else {
                    write(u, u, Bond.IMPLICIT);
                }
            }
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
    boolean prepare(int u, int p, Bond b) {
        visitedAt[u] = i++;
        tokens[u] = g.atom(u).token();

        boolean uncofiguredStereo = false;
        
        final int d = g.degree(u);
        for (int j=0; j<d; ++j) {
            final Edge e = g.edgeAt(u,j);
            int v = e.other(u);
            if (visitedAt[v] < 0) {
                uncofiguredStereo = prepare(v, u, e.bond(u)) || uncofiguredStereo;
            } else if (v != p && visitedAt[v] < visitedAt[u]) {
                cyclicEdge(v, u, e.bond(v));
            }
        }

        // Configure the topology using the traversal order
        // extended tetrahedral (@AL1/@AL2) doesn't have all neighbors visited 
        // until the end
        if (g.topologyOf(u).configuration().type() == Configuration.Type.ExtendedTetrahedral)
            return true;
        if (rings.containsKey(u)) {
            tokens[u].configure(g.topologyOf(u)
                                 .orderBy(localRank(u, p))
                                 .configuration());
        } else {
            tokens[u].configure(g.topologyOf(u)
                                 .orderBy(visitedAt)
                                 .configuration());
        }
        
        return uncofiguredStereo;
    }

    /**
     * Second traversal writes the bonds and atoms to the SMILES string.
     *
     * @param u a vertex
     * @param p previous vertex
     * @param b the bond from the previous vertex to this vertex
     */
    void write(int u, int p, Bond b) throws InvalidSmilesException {
        visitedAt[u] = i++;

        int remaining = g.degree(u);

        if (u != p)
            remaining--;

        // assign ring numbers
        if (rings.containsKey(u)) {
            for (RingClosure rc : rings.get(u)) {
                // as we are composing tokens, make sure apply in reverse
                int rnum = rnums.next();
                if (rc.register(rnum)) {
                    int v = rc.other(u);
                    tokens[u] = new RingNumberToken(new RingBondToken(tokens[u],
                                                                      rc.bond(u)),
                                                    rnum);
                    rnums.use(rnum);
                } else {
                    tokens[u] = new RingNumberToken(tokens[u],
                                                    rc.rnum);
                    rnums.free(rc.rnum);
                }
                remaining--;
            }
        }

        sb.append(b.token());
        tokens[u].append(sb);

        final int d = g.degree(u);
        for (int j=0; j<d; ++j) {
            final Edge e = g.edgeAt(u,j);
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
        final int d = g.degree(u);
        for (int j=0; j<d; ++j) {
            final Edge e = g.edgeAt(u,j);
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
    static String generate(final Graph g) throws InvalidSmilesException {
        return new Generator(g, new IterativeRingNumbering(1)).string();
    }

    /**
     * Convenience method for generating a SMILES string for the specified
     * chemical graph.
     *
     * @param g the graph to generate the SMILE for
     * @param visitedAt store when each atom was visited
     * @return SMILES gor the provided chemical graph
     */
    static String generate(final Graph g, int[] visitedAt) throws InvalidSmilesException {
        return new Generator(g, visitedAt, new IterativeRingNumbering(1)).string();
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

    /** Defines how ring numbering proceeds. */
    static interface RingNumbering {
        /**
         * The next ring number in the sequence.
         *
         * @return ring number
         */
        int next() throws InvalidSmilesException;

        /**
         * Mark the specified ring number as used.
         *
         * @param rnum ring number
         */
        void use(int rnum);

        /**
         * Mark the specified ring number as no longer used.
         *
         * @param rnum ring number
         */
        void free(int rnum);
        
        /** Reset ring number usage */
        void reset();
    }

    /** Labelling of ring opening/closures always using the lowest ring number. */
    static final class ReuseRingNumbering implements RingNumbering {

        private boolean[] used = new boolean[100];
        private final int offset;

        ReuseRingNumbering(int first) {
            this.offset = first;
        }

        @Override public int next() throws InvalidSmilesException {
            for (int i = offset; i < used.length; i++) {
                if (!used[i]) {
                    return i;
                }
            }
            throw new InvalidSmilesException("no available ring numbers");
        }

        @Override public void use(int rnum) {
            used[rnum] = true;
        }

        @Override public void free(int rnum) {
            used[rnum] = false;
        }

        @Override public void reset() {
            // do nothing 
        }
    }

    /**
     * Iterative labelling of ring opening/closures. Once the number 99 has been
     * used the number restarts using any free numbers.
     */
    static final class IterativeRingNumbering implements RingNumbering {

        private boolean[] used = new boolean[100];
        private final int offset;
        private       int pos;

        IterativeRingNumbering(int first) {
            this.offset = first;
            this.pos = offset;
        }

        @Override public int next() throws InvalidSmilesException {
            while (pos < 100 && used[pos])
                pos++;
            if (pos < 100)
                return pos;
            pos = offset;
            while (pos < 100 && used[pos])
                pos++;
            if (pos < 100)
                return pos;
            else
                throw new InvalidSmilesException("no more ring numbers can be assigned");
        }

        @Override public void use(int rnum) {
            used[rnum] = true;
        }

        @Override public void free(int rnum) {
            used[rnum] = false;
        }

        @Override public void reset() {
            pos = 1;
        }
    }
}
