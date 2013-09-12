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

import java.util.List;

/**
 * Utility to compute the orbital hybridization (Sp2) for delocalizing bonds.
 * Caveat - the values should not be taken as <i>golden</i> but rather serve the
 * purpose the demoralising bonds in SMILES molecules.
 *
 * @author John May
 */
enum Hybridization {
    /** Unknown hybridization - something went wrong. */
    Unknown,
    /** Sp orbital hybridization. */
    Sp,
    /** Sp2 orbital hybridization. */
    Sp2,
    /** Sp3 orbital hybridization. */
    Sp3,
    /** Sp3d orbital hybridization. */
    Sp3d,
    /** Sp3d2 orbital hybridization. */
    Sp3d2,
    /** Sp3d3 orbital hybridization. */
    Sp3d3;

    /**
     * Compute the hybridization of every atom in the graph {@code g}.
     *
     * @param g a graph
     * @return the hybridization of each atom
     */
    static Hybridization[] hybridizations(Graph g) {

        int ord = g.order();
        Hybridization[] hs = new Hybridization[ord];

        // Compute hybridization for each atom
        for (int u = 0; u < ord; u++)
            hs[u] = hybridization(g, u);

        // Now iteratively apply the rule that any atom with
        // one or more lone pairs is also Sp2 if it's next to
        // another Sp2 atom
        boolean changed;
        do {
            changed = false;
            for (int u = 0; u < ord; u++) {
                if (hs[u] != Sp2 && lonePairs(g, u) > 0 && attachedToSp2(hs, g, u)) {
                    hs[u] = Sp2;
                    changed = true;
                }
            }
        } while (changed);

        return hs;
    }

    /**
     * Compute orbital hybridization of the atom in.
     *
     * @param g graph
     * @param u vertex index
     * @return orbital hybridization
     */
    static Hybridization hybridization(Graph g, int u) {
        Atom atom = g.atom(u);
        if (!atom.element().aromatic())
            return Hybridization.Unknown;
        int v = valence(atom.element());
        int x = monovalant(g.edges(u)) + g.implHCount(u);
        int c = atom.charge() > 0 ? atom.charge() : 0;
        int a = atom.charge() < 0 ? -atom.charge() : 0;
        return hybridization(v, x, c, a);
    }

    /**
     * Compute the orbital hybridization.
     *
     * @param v valence electrons
     * @param x number of monovalent bonds
     * @param c cation charge
     * @param a anion charge
     * @return orbital hybridization
     */
    static Hybridization hybridization(int v, int x, int c, int a) {
        return map[(v + x - c + a) / 2];
    }

    /**
     * Number of monovalent bonds.
     *
     * @param es edges
     * @return number of monovalent bonds
     */
    static int monovalant(List<Edge> es) {
        int x = 0;
        for (final Edge e : es)
            if (e.bond() == Bond.SINGLE || e.bond() == Bond.IMPLICIT)
                x++;
        return x;
    }

    /**
     * Is the atom (vertex u) attached to an Sp2 atom.
     *
     * @param hs current hybridization values
     * @param g  graph
     * @param u  vertex
     * @return the atom label of vertex u as attached to an Sp2 atom
     */
    static boolean attachedToSp2(Hybridization[] hs, Graph g, int u) {
        for (Edge e : g.edges(u)) {
            if (hs[e.other(u)] == Sp2)
                return true;
        }
        return false;
    }

    /**
     * Number of valence electrons for possibly aromatic atoms. Note that the
     * method accepts strictly non-aromatic atoms.
     *
     * @param e element
     * @return the number of valent electrons
     */
    static int valence(Element e) {
        switch (e) {
            case Boron:
                return 3;
            case Carbon:
            case Silicon:    // nb. not 'strictly aromatic'
            case Germanium:  // nb. not 'strictly aromatic'  
            case Tin:        // nb. not 'strictly aromatic'  
                return 4;
            case Nitrogen:
            case Phosphorus:
            case Arsenic:
            case Antimony:   // nb. not 'strictly aromatic'
                return 5;
            case Oxygen:
            case Sulfur:
            case Selenium:
            case Tellurium:  // nb. not 'strictly aromatic'
                return 6;
        }
        throw new IllegalArgumentException("Unsupported element " + e);
    }

    /**
     * Number of lone pairs for a given atom.
     *
     * @param g graph
     * @param u vertex
     * @return number of lone pairs
     */
    static int lonePairs(Graph g, int u) {

        // bonded electron sum
        Atom atom = g.atom(u);
        
        if (!atom.element().aromatic())
            return 0;
        
        int sum = g.implHCount(u);
        List<Edge> es = g.edges(u);
        for (Edge e : es) {
            sum += e.bond().order();
        }
        int v = valence(g.atom(u).element()) + -atom.charge();
        return (v - sum) / 2;
    }

    /** Lookup from the computed hybridization value to enumeration type. */
    private static Hybridization[] map = new Hybridization[]{
            Unknown,  // 0
            Unknown,  // 1
            Sp,       // 2
            Sp2,      // 3
            Sp3,      // 4
            Sp3d,     // 5
            Sp3d2,    // 6
            Sp3d3,    // 7
    };
}
