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

/**
 * Defines properties of a SMILES atom. The organic and aromatic subsets are
 * provided as an enumeration for efficient reuse.
 *
 * @author John May
 */
interface Atom {

    /**
     * The element of the atom.
     *
     * @return element
     */
    Element element();

    /**
     * Whether this atom is aromatic.
     *
     * @return atom is aromatic
     */
    boolean aromatic();

    /**
     * Formal charge of the atom.
     *
     * @return formal charge
     */
    int charge();

    /**
     * Number of hydrogens this atom has. This value defines atoms with an
     * explicit hydrogen count of bracket atoms (e.g. [CH4]).
     *
     * @return hydrogen count
     * @throws IllegalArgumentException thrown if element is part of the organic
     *                                  subset and the number of hydrogens is
     *                                  implied by the bond order sum.
     */
    int hydrogens();

    /**
     * The class of the atom is defined as an integer value. The atom class is
     * specified for bracketed atoms and is prefixed by a colon.
     *
     * <blockquote><pre>
     *     [CH:1](C)([C])[H:2]
     * </pre></blockquote>
     *
     * @return class
     */
    int atomClass();

    static enum OrganicSubset implements Atom {
        Boron(Element.Boron),
        Carbon(Element.Carbon),
        Nitrogen(Element.Nitrogen),
        Oxygen(Element.Oxygen),
        Sulfur(Element.Sulfur),
        Phosphorus(Element.Phosphorus),
        Fluorine(Element.Fluorine),
        Chlorine(Element.Chlorine),
        Bromine(Element.Bromine),
        Iodine(Element.Iodine);

        private Element element;

        private OrganicSubset(Element element) {
            this.element = element;
        }

        @Override public Element element() {
            return element;
        }

        @Override public boolean aromatic() {
            return false;
        }

        @Override public int charge() {
            return 0;
        }

        @Override public int hydrogens() {
            throw new IllegalArgumentException("use bond order sum to determine implicit hydrogen count");
        }

        @Override public int atomClass() {
            return 0;
        }
    }

    static enum AromaticSubset implements Atom {
        Boron(Element.Boron),
        Carbon(Element.Carbon),
        Nitrogen(Element.Nitrogen),
        Oxygen(Element.Oxygen),
        Sulfur(Element.Sulfur),
        Phosphorus(Element.Phosphorus);

        private Element element;

        private AromaticSubset(Element element) {
            this.element = element;
        }

        @Override public Element element() {
            return element;
        }

        @Override public boolean aromatic() {
            return true;
        }

        @Override public int charge() {
            return 0;
        }

        @Override public int hydrogens() {
            throw new IllegalArgumentException("use bond order sum to determine implicit hydrogen count");
        }

        @Override public int atomClass() {
            return 0;
        }
    }


}
