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

import java.util.HashMap;
import java.util.Map;

/**
 * Enumeration of valid OpenSMILES elements.
 *
 * <h4>Organic subsets</h4> Several of the elements belong to the organic
 * subset. Atoms of an organic element type can be written just as their symbol
 * (see. <a href="http://www.opensmiles.org/opensmiles.html#orgsbst">Organic
 * Subset, OpenSMILES Specification</a>).
 *
 * <ul> <li>{@link #Unknown} (<code>*</code>)</li> <li>{@link #Boron}</li>
 * <li>{@link #Carbon}</li> <li>{@link #Nitrogen}</li> <li>{@link #Oxygen}</li>
 * <li>{@link #Fluorine}</li> <li>{@link #Phosphorus}</li> <li>{@link
 * #Sulfur}</li> <li>{@link #Chlorine}</li> <li>{@link #Bromine}</li> <li>{@link
 * #Iodine}</li> </ul>
 *
 * <h4>Usage</h4>
 *
 * Elements can be created by either using the value directly or by looking up
 * it's symbol. If the element may be aromatic the lower-case symbol can also be
 * used. For example the variable 'e' in the three statements below all have the
 * same value, {@link Element#Carbon}.
 *
 * <blockquote><pre>
 * Element e = Element.Carbon;
 * Element e = Element.ofSymbol("C");
 * Element e = Element.ofSymbol("c");
 * </pre></blockquote>
 *
 * When the symbol is invalid the result wil be null.
 * <blockquote><pre>
 * Element e = Element.ofSymbol("R1"); // e = null
 * </blockquote></pre>
 *
 * The {@link Element#Unknown} element can be used to represent generic/alias
 * atoms.
 * <blockquote><pre>
 * Element e = Element.Unknown;
 * Element e = Element.ofSymbol("*");
 * </blockquote></pre>
 *
 * To access the symbol of an already created element. Use {@link
 * Element#symbol()}.
 *
 * <blockquote><pre>
 * Atom    a = ...;
 * Element e = a.element();
 *
 * String  symbol = e.symbol();
 * </blockquote></pre>
 *
 * @author John May
 * @see <a href="http://www.opensmiles.org/opensmiles.html#inatoms">Atoms,
 *      OpenSMILES Specification</a>
 */
public enum Element {

    /** Unspecified/Unknown element (*) */
    Unknown("*", false, 0),

    Hydrogen("H"),
    Helium("He"),

    Lithium("Li"),
    Beryllium("Be"),
    Boron("B", true, 3),
    Carbon("C", true, 4),
    Nitrogen("N", true, 3, 5),
    Oxygen("O", true, 2),
    Fluorine("F", false, 1),
    Neon("Ne"),

    Sodium("Na"),
    Magnesium("Mg"),
    Aluminum("Al"),
    Silicon("Si", true, false),
    Phosphorus("P", true, 3, 5),
    Sulfur("S", true, 2, 4, 6),
    Chlorine("Cl", false, 1),
    Argon("Ar"),

    Potassium("K"),
    Calcium("Ca"),
    Scandium("Sc"),
    Titanium("Ti"),
    Vanadium("V"),
    Chromium("Cr"),
    Manganese("Mn"),
    Iron("Fe"),
    Cobalt("Co"),
    Nickel("Ni"),
    Copper("Cu"),
    Zinc("Zn"),
    Gallium("Ga"),
    Germanium("Ge", true, false),
    Arsenic("As", true),
    Selenium("Se", true),
    Bromine("Br", false, 1),
    Krypton("Kr"),

    Rubidium("Rb"),
    Strontium("Sr"),
    Yttrium("Y"),
    Zirconium("Zr"),
    Niobium("Nb"),
    Molybdenum("Mo"),
    Technetium("Tc"),
    Ruthenium("Ru"),
    Rhodium("Rh"),
    Palladium("Pd"),
    Silver("Ag"),
    Cadmium("Cd"),
    Indium("In"),
    Tin("Sn", true, false),
    Antimony("Sb", true, false),
    Tellurium("Te", true, false),
    Iodine("I", false, 1),
    Xenon("Xe"),

    Cesium("Cs"),
    Barium("Ba"),
    Lutetium("Lu"),
    Hafnium("Hf"),
    Tantalum("Ta"),
    Tungsten("W"),
    Rhenium("Re"),
    Osmium("Os"),
    Iridium("Ir"),
    Platinum("Pt"),
    Gold("Au"),
    Mercury("Hg"),
    Thallium("Tl"),
    Lead("Pb"),
    Bismuth("Bi", true, false),
    Polonium("Po"),
    Astatine("At"),
    Radon("Rn"),

    Francium("Fr"),
    Radium("Ra"),
    Lawrencium("Lr"),
    Rutherfordium("Rf"),
    Dubnium("Db"),
    Seaborgium("Sg"),
    Bohrium("Bh"),
    Hassium("Hs"),
    Meitnerium("Mt"),
    Darmstadtium("Ds"),
    Roentgenium("Rg"),
    Copernicium("Cn"),
    Flerovium("Fl"),
    Livermorium("Lv"),

    Lanthanum("La"),
    Cerium("Ce"),
    Praseodymium("Pr"),
    Neodymium("Nd"),
    Promethium("Pm"),
    Samarium("Sm"),
    Europium("Eu"),
    Gadolinium("Gd"),
    Terbium("Tb"),
    Dysprosium("Dy"),
    Holmium("Ho"),
    Erbium("Er"),
    Thulium("Tm"),
    Ytterbium("Yb"),

    Actinium("Ac"),
    Thorium("Th"),
    Protactinium("Pa"),
    Uranium("U"),
    Neptunium("Np"),
    Plutonium("Pu"),
    Americium("Am"),
    Curium("Cm"),
    Berkelium("Bk"),
    Californium("Cf"),
    Einsteinium("Es"),
    Fermium("Fm"),
    Mendelevium("Md"),
    Nobelium("No");

    /** The symbol of the element. */
    private final String symbol;

    /**
     * Default valence information - only present if the atom is part of the
     * organic subset.
     */
    private final int[] valence;

    private final int[] electrons;

    /**
     * Whether the element can be aromatic.  and whether the element is aromatic
     * by the OpenSMILES specification.
     */
    private final boolean aromatic, strictlyAromatic;

    /** Look up of elements by symbol */
    private static final Map<String, Element> elementMap
            = new HashMap<String, Element>();

    static {
        for (Element element : values()) {
            if (element.aromatic())
                elementMap.put(element.symbol().toLowerCase(), element);
            elementMap.put(element.symbol(), element);
        }
    }

    private Element(String symbol) {
        this(symbol, false, false, null);
    }

    private Element(String symbol,
                    boolean aromatic) {
        this(symbol, aromatic, true, null);
    }

    private Element(String symbol,
                    boolean aromatic,
                    boolean strictlyAromatic) {
        this(symbol, aromatic, strictlyAromatic, null);
    }

    private Element(String symbol,
                    boolean aromatic,
                    int... valence) {
        this(symbol, aromatic, true, valence);
    }

    private Element(String symbol,
                    boolean aromatic,
                    boolean strictlyAromatic,
                    int... valence) {
        this.symbol = symbol;
        this.valence = valence;
        if (valence != null) {
            this.electrons = new int[valence.length];
            for (int i = 0; i < valence.length; i++) {
                electrons[i] = valence[i] * 2;
            }
        } else {
            this.electrons = null;
        }
        this.aromatic = aromatic;
        this.strictlyAromatic = strictlyAromatic;
    }

    /**
     * Access the symbol of the element.
     *
     * @return element symbol
     */
    public String symbol() {
        return symbol;
    }

    /**
     * Can the element be aromatic (as defined by the OpenSMILES
     * specification).
     *
     * @return whether the element may be aromatic
     */
    boolean aromatic() {
        return aromatic;
    }

    /**
     * Is the element strictly aromatic as per OpenSMILES specification? The
     * OpenSMILES specification defines a strict set of symbols as aromatic. In
     * practise others are encountered (e.g. '[te]1ccccc1') - these can be
     * parsed but are not supported by the OpenSMILES specification. This method
     * can allows one to check if an aromatic atom is aromatic in the
     * specification.
     *
     * @return the element is aromatic
     */
    boolean strictlyAromatic() {
        return strictlyAromatic;
    }

    /**
     * Is the element a member of the organic subset and can be written without
     * brackets. If the element is both organic and aromatic is a member of the
     * aromatic subset and can still be written without brackets.
     *
     * @return the element can be written without brackets
     */
    boolean organic() {
        return valence != null;
    }

    /**
     * Determine the number of implied hydrogens an organic (or aromatic) subset
     * atom has based on it's bond order sum. The valances for the organic
     * elements (B, C, N, O, P, S, F, Cl, Br and I) are defined in the
     * OpenSMILES specification.
     *
     * @param sum bond order sum
     * @return the number of implied hydrogens
     * @throws IllegalArgumentException the element was not a member of the
     *                                  organic subset and did not have default
     *                                  valence information
     * @deprecated use {@link #availableElectrons(int)} or {@link
     *             #availableDelocalisedElectrons(int)}
     */
    @Deprecated int implicitHydrogens(int sum) {
        if (!organic())
            throw new IllegalArgumentException("inorganic atom, no preset valence: " + this);

        for (final int v : valence)
            if (sum <= v) return v - sum;

        // bond order sum exceeds or equals maximum valance
        return 0;
    }

    /**
     * Determine the number of available electrons which could be bonding to
     * implicit hydrogens. This include electrons donated from the hydrogen.
     * <br/>
     *
     * The central carbon of {@code C-C=C} 6 bonded electrons - using SMILES
     * default valence there must be 2 electrons involved in bonding an implicit
     * hydrogen (i.e. there is a single bond to a hydrogen).
     *
     * @param bondElectronSum the sum of the bonded electrons
     * @return number of electrons which could be involved with bonds to
     *         hydrogen
     */
    int availableElectrons(int bondElectronSum) {
        for (final int e : electrons)
            if (bondElectronSum <= e)
                return e - bondElectronSum;
        return 0;
    }

    /**
     * Determine the number of available electrons which could be bonding to
     * implicit hydrogens for an aromatic atom with delocalized bonds. This
     * include electrons donated from the hydrogen. <br/>
     *
     * Instead of checking higher valence states only the lowest is checked. For
     * example nitrogen has valence 3 and 5 but in a delocalized system only the
     * lowest (3) is used. The electrons which would allow bonding of implicit
     * hydrogens in the higher valence states are donated to the aromatic system
     * and thus cannot be <i>reached</i>. Using a generalisation that an 
     * aromatic bond as 3 electrons we reached the correct value for multi 
     * valence aromatic elements. <br/>
     * 
     * <blockquote><pre>
     *     c1c[nH]cn1    the aromatic subset nitrogen is bonded to two aromatic
     *                   nitrogen bond order sum of 3 (6 electrons) there are
     *                   no implicit hydrogens
     *                   
     *     c1cc2ccccn2c1 the nitrogen has three aromatic bond 4.5 bond order
     *                   (9 electrons) - as we only check the lowest valence
     *                   (3 - 4.5) < 0 so there are 0 implicit hydrogens
     *                   
     *     c1ccpcc1      the phosphorus has 2 aromatic bond (bond order sum 3)
     *                   and the lowest valence is '3' - there are no implicit
     *                   hydrogens
     *                   
     *     oc1ccscc1     the sulphur has two aromatic bonds (bond order sum 3)
     *                   the lowest valence is '2' - 3 > 2 so there are no 
     *                   implicit hydrogens
     *                   
     *     oc1ccscc1     the oxygen has a single aromatic bond, the default 
     *                   valence of oxygen in the specification is '2' there
     *                   are no hydrogens (2 - 1.5 = 0.5).
     * </pre></blockquote>
     *
     * @param bondElectronSum the sum of the bonded electrons
     * @return number of electrons which could be involved with bonds to
     *         hydrogen
     */
    int availableDelocalisedElectrons(int bondElectronSum) {
        if (bondElectronSum <= electrons[0])
            return electrons[0] - bondElectronSum;
        return 0;
    }

    /**
     * Given an element symbol, provide the element for that symbol. If no
     * symbol was found then null is returned.
     *
     * @param symbol the element symbol
     * @return element for the symbol, or null if none found
     */
    public static Element ofSymbol(final String symbol) {
        return elementMap.get(symbol);
    }

    /**
     * Read an element and progress the character buffer. If the element was not
     * read then a 'null' element is returned.
     *
     * @param buffer a character buffer
     * @return the element, or null
     */
    static Element read(final CharBuffer buffer) {
        if (!buffer.hasRemaining())
            return null;
        char c = buffer.get();
        if (buffer.hasRemaining() && buffer.next() >= 'a' && buffer
                .next() <= 'z') {
            return elementMap.get(new String(new char[]{c, buffer.get()}));
        }
        return elementMap.get(Character.toString(c));
    }
}
