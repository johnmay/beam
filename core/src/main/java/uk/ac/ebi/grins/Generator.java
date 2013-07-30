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

import java.util.Arrays;
import java.util.List;

/**
 * Generate a SMILES line notation for a given chemical graph.
 *
 * @author John May
 */
final class Generator {

    private final ChemicalGraph g;
    private final StringBuilder sb;

    private final int[] visited;
    private       int   i;

    /**
     * Create a new generator the given chemical graph.
     *
     * @param g chemical graph
     */
    Generator(ChemicalGraph g) {
        this.g = g;
        this.sb = new StringBuilder(g.order() * 2);
        this.visited = new int[g.order()];

        // prepare ring closures and topologies
        Arrays.fill(visited, -1);
        for (int u = 0; u < g.order(); u++) {
            if (visited[u] < 0)
                prepare(u, u == 0 ? Bond.IMPLICIT : Bond.DOT);
        }

        // write notation
        i = 0;
        Arrays.fill(visited, -1);
        for (int u = 0; u < g.order(); u++) {
            if (visited[u] < 0)
                write(u, u == 0 ? Bond.IMPLICIT : Bond.DOT);
        }
    }

    void prepare(int u, Bond b) {
        visited[u] = i++;
        for (Edge e : g.edges(u)) {
            int v = e.other(u);
            if (visited[v] < 0) {
                prepare(v, e.bond(u));
            } else if (visited[v] != visited[u] - 1) {
                // ring closure
            }
        }
    }

    void write(int u, Bond b) {

        visited[u] = i++;
        sb.append(b.symbol());
        sb.append(g.atom(u).element().symbol());

        List<Edge> es = g.edges(u);
        for (int i = 0; i < es.size(); i++) {
            Edge e = es.get(i);
            int v = e.other(u);
            if (visited[v] < 0) {
                write(v, e.bond(u));
            } else if (visited[v] != visited[u] - 1) {
                // ring closure
            }
        }
    }

    String string() {
        return sb.toString();
    }

    static String generate(ChemicalGraph g) {
        return new Generator(g).string();
    }
}
