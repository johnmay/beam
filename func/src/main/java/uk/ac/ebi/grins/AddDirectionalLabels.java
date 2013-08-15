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
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Given a molecule with bond-based double bond configurations - add directional
 * labels to edges which do not have it assigned. For example the molecule
 * {@code NC(/C)=C\C} has no directional label between the nitrogen and the
 * carbon. Applying this procedure will 'fill-in' missing directional
 * information on the edge - {@code N/C(/C)=C\C}. <br/>
 *
 * If required the directional labels in conjugated systems may be adjusted to
 * allow for full-specification. Attempting to assign a directional label to the
 * central carbon of {@code F/C=C(/F)C(/F)=C/F} creates a conflict. This conflict
 * will be resolved by flipping the labels on the second double-bond - {@code
 * F/C=C(/F)\C(\F)=C\F}.
 *
 * @author John May
 */
final class AddDirectionalLabels
        extends AbstractFunction<ChemicalGraph, ChemicalGraph> {

    /**
     * Transform all implicit up/down to their explicit type. The original graph
     * is unmodified
     *
     * @param g a chemical graph
     * @return new chemical graph but with all explicit bonds
     */
    public ChemicalGraph apply(final ChemicalGraph g)
            throws InvalidSmilesException {

        ChemicalGraph h = new ChemicalGraph(g.order());

        // copy atom/topology information this is unchanged
        for (int u = 0; u < g.order(); u++) {
            h.addAtom(g.atom(u));
            h.addTopology(g.topologyOf(u));
        }

        Map<Edge, Edge> replacements = new HashMap<Edge, Edge>();
        Set<Edge> remaining = new HashSet<Edge>();

        // change edges (only changed added to replacement)
        for (int u = 0; u < g.order(); u++) {
            for (final Edge e : g.edges(u)) {
                if (e.other(u) > u && e.bond() == Bond.DOUBLE) {
                    remaining.add(e);
                }
            }
        }

        boolean altered = false;
        do {
            for (Edge e : new ArrayList<Edge>(remaining)) {
                if (!replaceImplWithExpl(g, e, replacements)) {
                    remaining.remove(e);
                    altered = true;
                }
            }
        } while (altered && !remaining.isEmpty());

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
    private boolean replaceImplWithExpl(ChemicalGraph g,
                                        Edge e,
                                        Map<Edge, Edge> acc)
            throws InvalidSmilesException {

        int u = e.either(), v = e.other(u);

        boolean uDone = replaceImplWithExpl(g, e, u, acc);
        boolean vDone = replaceImplWithExpl(g, e, v, acc);

        return uDone || vDone;
    }

    /**
     * Given a double bond edge traverse the neighbors of one of the endpoints
     * and accumulate any explicit replacements in the 'acc' accumulator.
     *
     * @param g   the chemical graph
     * @param e   a edge in the graph ('double bond type')
     * @param u   a endpoint of the edge 'e'
     * @param acc accumulator for new edges
     * @return does the edge 'e' need to be reconsidered later
     * @throws InvalidSmilesException thrown if the edge could not be converted
     */
    private boolean replaceImplWithExpl(ChemicalGraph g,
                                        Edge e,
                                        int u,
                                        Map<Edge, Edge> acc)
            throws InvalidSmilesException {

        Edge implicit = null;
        Edge explicit = null;

        for (Edge f : g.edges(u)) {
            Edge f2 = acc.containsKey(f) ? acc.get(f) : f;
            switch (f2.bond(u)) {
                case SINGLE:
                case IMPLICIT:
                    if (implicit != null)
                        return true;
                    implicit = f;
                    break;
                case DOUBLE:
                    if (!f.equals(e))
                        return false;
                    break;
                case UP:
                case DOWN:
                    if (explicit != null) {

                        if (acc.containsKey(explicit))
                            explicit = acc.get(explicit);

                        // original bonds are invalid
                        if ((f.bond() == Bond.UP || f.bond() == Bond.DOWN) &&
                                explicit.bond(u).inverse() != f.bond(u)) {
                            throw new InvalidSmilesException("invalid double bond configuration");
                        }

                        if (explicit.bond(u).inverse() != f2.bond(u)) {
                            acc.put(f, f2.inverse());
                            BitSet visited = new BitSet();
                            visited.set(u);
                            invertExistingDirectionalLabels(g, visited, acc, f2
                                    .other(u));
                        }
                        return false;
                    }
                    explicit = f;
                    break;
            }
        }

        // no implicit or explicit bond? don't do anything
        if (implicit == null || explicit == null)
            return false;

        if (acc.containsKey(explicit))
            explicit = acc.get(explicit);

        int v = implicit.other(u);

        Edge existing = acc.put(implicit, new Edge(u,
                                                   v,
                                                   explicit.bond(u)
                                                           .inverse()));

        if (existing != null && existing.bond(u) != explicit.bond(u).inverse())
            throw new InvalidSmilesException("unable to assign explict type for " + implicit);

        return false;
    }

    private void invertExistingDirectionalLabels(ChemicalGraph g,
                                                 BitSet visited,
                                                 Map<Edge, Edge> replacement,
                                                 int u) {
        visited.set(u);
        if (g.topologyOf(u) == null)
            return;
        for (Edge e : g.edges(u)) {
            int v = e.other(u);
            if (!visited.get(v)) {
                Edge f = replacement.get(e);
                if (f != null) {
                    replacement.put(e,
                                    f.inverse());
                } else {
                    replacement.put(e,
                                    e.inverse());
                }
                invertExistingDirectionalLabels(g, visited, replacement, v);
            }
        }
    }

}
