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

import java.util.HashMap;
import java.util.Map;

/**
 * Given a molecule with implicit double bond configurations add UP/DOWN bond
 * types to the implicit bonds. For example the explicit form of {@code
 * NC(/C)=C\C} is {@code N/C(/C)=C\C}.
 *
 * @author John May
 */
final class AddUpDownBonds {

    /**
     * Transform all implicit up/down to their explicit type. The original graph
     * is unmodified
     *
     * @param g a chemical graph
     * @return new chemical graph but with all explicit bonds
     */
    public ChemicalGraph transform(final ChemicalGraph g)
            throws InvalidSmilesException {

        ChemicalGraph h = new ChemicalGraph(g.order());

        // copy atom/topology information this is unchanged
        for (int u = 0; u < g.order(); u++) {
            h.addAtom(g.atom(u));
            h.addTopology(g.topologyOf(u));
        }

        Map<Edge, Edge> replacements = new HashMap<Edge, Edge>();

        // change edges (only changed added to replacement)
        for (int u = 0; u < g.order(); u++) {
            for (final Edge e : g.edges(u)) {
                if (e.other(u) > u && e.bond() == Bond.DOUBLE) {
                    replaceImplWithExpl(g, e, replacements);
                }
            }
        }

        // append the edges, replacing any which need to be changed
        for (int u = 0; u < g.order(); u++) {
            for (Edge e : g.edges(u)) {
                if (e.other(u) > u) {
                    Edge replacement = replacements.get(e);
                    if (replacement != null)
                        e = replacement;
                    h.addEdge(e);
                }
            }
        }

        return h;
    }

    /**
     * Given a double bond edge traverse the neighbors of both endpoints and
     * accumulate any explicit replacements in the 'acc' accumulator.
     *
     * @param g   the chemical graph
     * @param e   a edge in the graph ('double bond type')
     * @param acc accumulator for new edges
     * @throws InvalidSmilesException thrown if the edge could not be converted
     */
    private void replaceImplWithExpl(ChemicalGraph g,
                                     Edge e,
                                     Map<Edge, Edge> acc)
            throws InvalidSmilesException {

        int u = e.either(), v = e.other(u);


        replaceImplWithExpl(g, e, u, acc);
        replaceImplWithExpl(g, e, v, acc);

    }

    /**
     * Given a double bond edge traverse the neighbors of one of the endpoints
     * and accumulate any explicit replacements in the 'acc' accumulator.
     *
     * @param g   the chemical graph
     * @param e   a edge in the graph ('double bond type')
     * @param u   a endpoint of the edge 'e'
     * @param acc accumulator for new edges
     * @throws InvalidSmilesException thrown if the edge could not be converted
     */
    private void replaceImplWithExpl(ChemicalGraph g,
                                     Edge e,
                                     int u,
                                     Map<Edge, Edge> acc)
            throws InvalidSmilesException {

        Edge implicit = null;
        Edge explicit = null;

        for (Edge f : g.edges(u)) {
            switch (f.bond()) {
                case SINGLE:
                case IMPLICIT:
                    if (implicit != null)
                        return;
                    implicit = f;
                    break;
                case DOUBLE:
                    if (!f.equals(e))
                        return;
                    break;
                case UP:
                case DOWN:
                    if (explicit != null) {
                        if (explicit.bond(u).inverse() != f.bond(u))
                            throw new InvalidSmilesException("invalid double bond configuration");
                        return;
                    }
                    explicit = f;
                    break;
            }
        }

        // no implicit or explicit bond? don't do anything
        if (implicit == null || explicit == null)
            return;


        int v = implicit.other(u);

        Edge existing = acc.put(implicit, new Edge(u,
                                                   v,
                                                   explicit.bond(u)
                                                           .inverse()));

        if (existing != null && existing.bond() != explicit.bond(u).inverse())
            throw new InvalidSmilesException("unable to assign explict type for " + implicit);

    }

}
