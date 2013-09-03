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

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.AbstractMap;
import java.util.Arrays;
import java.util.EnumSet;
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
    Unknown("*", 0),

    Hydrogen("H"),
    Helium("He"),

    Lithium("Li"),
    Beryllium("Be"),
    Boron("B", 3),
    Carbon("C", 4),
    Nitrogen("N", 3, 5),
    Oxygen("O", 2),
    Fluorine("F", 1),
    Neon("Ne"),

    Sodium("Na"),
    Magnesium("Mg"),
    Aluminum("Al"),
    Silicon("Si"),
    Phosphorus("P", 3, 5),
    Sulfur("S", 2, 4, 6),
    Chlorine("Cl", 1),
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
    Germanium("Ge"),
    Arsenic("As"),
    Selenium("Se"),
    Bromine("Br", 1),
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
    Tin("Sn"),
    Antimony("Sb"),
    Tellurium("Te"),
    Iodine("I", 1),
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
    Bismuth("Bi"),
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

    /** Look up of elements by symbol */
    private static final Map<String, Element> elementMap
            = new HashMap<String, Element>();

    /** Provide verification of valence/charge values. */
    private ElementCheck defaults = ElementCheck.NO_CHECK;

    static {
        for (Element element : values()) {
            if (element.aromatic())
                elementMap.put(element.symbol().toLowerCase(), element);
            elementMap.put(element.symbol(), element);
        }

        // load normal ranges from 'element-defaults.txt' and set for the
        // elements
        for (Map.Entry<String, ElementCheck> e : loadDefaults().entrySet()) {
            elementMap.get(e.getKey()).defaults = e.getValue();
        }
    }

    private Element(String symbol) {
        this(symbol, null);
    }

    private Element(String symbol,
                    int... valence) {
        this.symbol = symbol;
        this.valence = valence;
        if (valence != null) {
            this.electrons = new int[valence.length];
            for (int i = 0; i < valence.length; i++) {
                electrons[i] = valence[i] * 2;
            }
        }
        else {
            this.electrons = null;
        }
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
     * Can the element be aromatic. This definition is very loose and includes
     * elements which are not part of the Daylight, OpenSMILES specification. To
     * test if ane element is aromatic by the specification use {@link
     * #aromatic(uk.ac.ebi.beam.Element.AromaticSpecification)}.
     *
     * @return whether the element may be aromatic
     */
    boolean aromatic() {
        return aromatic(AromaticSpecification.General);
    }

    /**
     * Can the element be aromatic in accordance with a given specification.
     *
     * @param spec such {@link uk.ac.ebi.beam.Element.AromaticSpecification#Daylight},
     *             {@link uk.ac.ebi.beam.Element.AromaticSpecification#OpenSmiles}
     * @return the element is accepted as being aromatic by that scheme
     */
    boolean aromatic(AromaticSpecification spec) {
        return spec.contains(this);
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
     * Verify whether the given valence and charge are 'normal' for the
     * element.
     *
     * @param v valence (bond order order sum)
     * @param q charge
     * @return whether the valence and charge are valid
     */
    boolean verify(int v, int q) {
        // table driven verification (see. element-defaults.txt)
        return defaults.verify(v, q);
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

    static Map<String, ElementCheck> loadDefaults() {
        Map<String, ElementCheck> checks = new HashMap<String, ElementCheck>(200);
        try {
            InputStream in = Element.class.getResourceAsStream("element-defaults.txt");
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String line = null;
            while ((line = br.readLine()) != null) {
                if (line.length() == 0 || line.charAt(0) == '-') // empty line or comment
                    continue;
                Map.Entry<String, ElementCheck> entry = load(line);
                checks.put(entry.getKey(), entry.getValue());
            }
            br.close();
        } catch (Exception e) {
            System.err.println("error whilst loading element-defaults.txt: " + e);
        }
        return checks;
    }

    static Map.Entry<String, ElementCheck> load(String line) {
        String[] data = line.split("\\s+");
        String symbol = data[0];
        int electrons = Integer.parseInt(data[3]);
        ValenceCheck valenceCheck = ValenceCheck.parse(data[1], electrons);
        ChargeCheck chargeCheck = ChargeCheck.parse(data[2]);
        return new AbstractMap.SimpleEntry<String, ElementCheck>(symbol,
                                                                 new ElementCheck(valenceCheck, chargeCheck));
    }

    private static final class ElementCheck {
        private final ValenceCheck valenceCheck;
        private final ChargeCheck  chargeCheck;

        private ElementCheck(ValenceCheck valenceCheck, ChargeCheck chargeCheck) {
            this.valenceCheck = valenceCheck;
            this.chargeCheck = chargeCheck;
        }

        boolean verify(int v, int q) {
            return chargeCheck.verify(q) && valenceCheck.verify(v, q);
        }

        private static final ElementCheck NO_CHECK = new ElementCheck(NoValenceCheck.INSTANCE,
                                                                      ChargeCheck.NONE);

        @Override public String toString() {
            return chargeCheck + ", " + valenceCheck;
        }
    }

    private static abstract class ValenceCheck {

        abstract boolean verify(final int v, final int q);

        static ValenceCheck parse(String line, int nElectrons) {
            String[] vs = line.split(",");
            if (vs.length == 1) {
                if (vs[0].equals("n/a")) {
                    return NoValenceCheck.INSTANCE;
                }
                else if (vs[0].charAt(0) == '(') {
                    return new FixedValence(Integer.parseInt(vs[0].substring(1, vs[0].length() - 1)));
                }
                else if (vs[0].charAt(0) == '[') {
                    return new NeutralValence(Integer.parseInt(vs[0].substring(1, vs[0].length() - 1)));
                }
                else {
                    return new ChargeAdjustedValence(Integer.parseInt(vs[0]), nElectrons);
                }
            }
            ValenceCheck[] valences = new ValenceCheck[vs.length];
            for (int i = 0; i < vs.length; i++) {
                valences[i] = parse(vs[i], nElectrons);
            }

            return new MultiValenceCheck(valences);
        }
    }

    private static final class ChargeAdjustedValence extends ValenceCheck {
        private final int valence, nElectrons;

        private ChargeAdjustedValence(int valence, int nElectrons) {
            this.valence = valence;
            this.nElectrons = nElectrons;
        }

        @Override public boolean verify(int v, int q) {
            if (nElectrons == 2 && valence + q > nElectrons - q)  // Group 2 exception
                return v == nElectrons - q;
            return valence + q == v;
        }

        @Override public String toString() {
            return "Charge(" + valence + ")";
        }
    }

    /** A valence check which is only valid at netural charge */
    private static final class NeutralValence extends ValenceCheck {
        private final int valence;

        private NeutralValence(int valence) {
            this.valence = valence;
        }

        @Override public boolean verify(int v, int q) {
            return q == 0 && v == valence;
        }

        @Override public String toString() {
            return "Neutral(" + valence + ")";
        }
    }

    private static final class FixedValence extends ValenceCheck {
        private final int valence;

        private FixedValence(int valence) {
            this.valence = valence;
        }

        @Override public boolean verify(int v, int q) {
            return valence == v;
        }

        @Override public String toString() {
            return "Fixed(" + valence + ")";
        }
    }

    private static final class MultiValenceCheck extends ValenceCheck {

        private final ValenceCheck[] valences;

        private MultiValenceCheck(ValenceCheck[] valences) {
            this.valences = valences;
        }

        @Override public boolean verify(int v, int q) {
            for (ValenceCheck vc : valences) {
                if (vc.verify(v, q)) {
                    return true;
                }
            }
            return false;
        }

        @Override public String toString() {
            return Arrays.toString(valences);
        }
    }

    private static final class NoValenceCheck extends ValenceCheck {
        @Override boolean verify(int v, int q) {
            return true;
        }

        private static final ValenceCheck INSTANCE = new NoValenceCheck();
    }

    private static final class ChargeCheck {

        private final int lo, hi;

        private ChargeCheck(int lo, int hi) {
            this.lo = lo;
            this.hi = hi;
        }

        boolean verify(final int q) {
            return lo <= q && q <= hi;
        }

        static ChargeCheck parse(String range) {
            if (range.equals("n/a"))
                return NONE;
            String[] data = range.split(",");
            int lo = Integer.parseInt(data[0]);
            int hi = Integer.parseInt(data[1]);
            return new ChargeCheck(lo, hi);
        }

        private static final ChargeCheck NONE = new ChargeCheck(Integer.MIN_VALUE, Integer.MAX_VALUE);

        @Override public String toString() {
            return lo + " < q < " + hi;
        }
    }

    /**
     * Stores which elements the Daylight and OpenSMILES specification consider
     * to be aromatic. The General scheme is what might be encountered 'in the
     * wild'.
     */
    enum AromaticSpecification {

        Daylight(Carbon,
                 Nitrogen,
                 Oxygen,
                 Sulfur,
                 Phosphorus,
                 Arsenic,
                 Selenium),

        OpenSmiles(Boron,
                   Carbon,
                   Nitrogen,
                   Oxygen,
                   Sulfur,
                   Phosphorus,
                   Arsenic,
                   Selenium),

        General(Boron,
                Carbon,
                Nitrogen,
                Oxygen,
                Sulfur,
                Phosphorus,
                Arsenic,
                Selenium,
                Silicon,
                Germanium,
                Tin,
                Antimony,
                Tellurium,
                Bismuth);

        private EnumSet<Element> elements;

        AromaticSpecification(Element... es) {
            this.elements = EnumSet.noneOf(Element.class);
            for (Element e : es)
                elements.add(e);
        }

        boolean contains(Element e) {
            return elements.contains(e);
        }
    }
}
