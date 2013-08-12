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

/** @author John May */
public final class AtomBuilder {

    private final Element e;
    private int isotope = -1,
            hCount      = 0,
            charge      = 0,
            atomClass   = 0;
    private boolean aromatic;

    private AtomBuilder(Element e) {
        this.e = e;
    }

    public static AtomBuilder create(Element e) {
        if (e == null)
            throw new NullPointerException("no element provided");
        return new AtomBuilder(e);
    }

    public static AtomBuilder create(String symbol) {
        Element e = Element.ofSymbol(symbol);
        if (e == null)
            e = Element.Unknown;
        AtomBuilder ab = new AtomBuilder(e);
        if (symbol != null
             && !symbol.isEmpty()
             && Character.isLowerCase(symbol.charAt(0))) {
            if (!e.aromatic())
                throw new IllegalArgumentException("Attempting to create an aromatic atom for an element which cannot be aromatic");
            ab.aromatic = true;
        }
        return ab;
    }

    public Atom build() {
        return new AtomImpl.BracketAtom(isotope,
                                        e,
                                        hCount,
                                        charge,
                                        atomClass,
                                        aromatic);
    }

}
