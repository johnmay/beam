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

import static uk.ac.ebi.grins.Configuration.Type.ExtendedTetrahedral;
import static uk.ac.ebi.grins.Configuration.Type.Implicit;
import static uk.ac.ebi.grins.Configuration.Type.None;
import static uk.ac.ebi.grins.Configuration.Type.Octahedral;
import static uk.ac.ebi.grins.Configuration.Type.SquarePlanar;
import static uk.ac.ebi.grins.Configuration.Type.Tetrahedral;
import static uk.ac.ebi.grins.Configuration.Type.TrigonalBipyramidal;

/**
 * Enumeration of valid relative configurations for atoms.
 *
 * @author John May
 */
enum Configuration {

    // atoms have unknown configuration
    UNKNOWN(None, ""),

    // shorthand for TH, AL or atom-stereo double bond
    ANTI_CLOCKWISE(Implicit, "@"),
    CLOCKWISE(Implicit, "@@"),

    // tetrahedral
    TH1(Tetrahedral, "@TH1", ANTI_CLOCKWISE),
    TH2(Tetrahedral, "@TH2", CLOCKWISE),

    // extended tetrahedral, allene-like
    AL1(ExtendedTetrahedral, "@AL1", ANTI_CLOCKWISE),
    AL2(ExtendedTetrahedral, "@AL2", CLOCKWISE),

    // square planar
    SP1(SquarePlanar, "@SP1", UNKNOWN),
    SP2(SquarePlanar, "@SP2", UNKNOWN),
    SP3(SquarePlanar, "@SP3", UNKNOWN),

    // trigonal bipyramidal
    TB1(TrigonalBipyramidal, "@TB1", UNKNOWN),
    TB2(TrigonalBipyramidal, "@TB2", UNKNOWN),
    TB3(TrigonalBipyramidal, "@TB3", UNKNOWN),
    TB4(TrigonalBipyramidal, "@TB4", UNKNOWN),
    TB5(TrigonalBipyramidal, "@TB5", UNKNOWN),
    TB6(TrigonalBipyramidal, "@TB6", UNKNOWN),
    TB7(TrigonalBipyramidal, "@TB7", UNKNOWN),
    TB8(TrigonalBipyramidal, "@TB8", UNKNOWN),
    TB9(TrigonalBipyramidal, "@TB9", UNKNOWN),
    TB10(TrigonalBipyramidal, "@TB10", UNKNOWN),
    TB11(TrigonalBipyramidal, "@TB11", UNKNOWN),
    TB12(TrigonalBipyramidal, "@TB12", UNKNOWN),
    TB13(TrigonalBipyramidal, "@TB13", UNKNOWN),
    TB14(TrigonalBipyramidal, "@TB14", UNKNOWN),
    TB15(TrigonalBipyramidal, "@TB15", UNKNOWN),
    TB16(TrigonalBipyramidal, "@TB16", UNKNOWN),
    TB17(TrigonalBipyramidal, "@TB17", UNKNOWN),
    TB18(TrigonalBipyramidal, "@TB18", UNKNOWN),
    TB19(TrigonalBipyramidal, "@TB19", UNKNOWN),
    TB20(TrigonalBipyramidal, "@TB20", UNKNOWN),

    // octahedral
    OH1(Octahedral, "@OH1", UNKNOWN),
    OH2(Octahedral, "@OH2", UNKNOWN),
    OH3(Octahedral, "@OH3", UNKNOWN),
    OH4(Octahedral, "@OH4", UNKNOWN),
    OH5(Octahedral, "@OH5", UNKNOWN),
    OH6(Octahedral, "@OH6", UNKNOWN),
    OH7(Octahedral, "@OH7", UNKNOWN),
    OH8(Octahedral, "@OH8", UNKNOWN),
    OH9(Octahedral, "@OH9", UNKNOWN),
    OH10(Octahedral, "@OH10", UNKNOWN),
    OH11(Octahedral, "@OH11", UNKNOWN),
    OH12(Octahedral, "@OH12", UNKNOWN),
    OH13(Octahedral, "@OH13", UNKNOWN),
    OH14(Octahedral, "@OH14", UNKNOWN),
    OH15(Octahedral, "@OH15", UNKNOWN),
    OH16(Octahedral, "@OH16", UNKNOWN),
    OH17(Octahedral, "@OH17", UNKNOWN),
    OH18(Octahedral, "@OH18", UNKNOWN),
    OH19(Octahedral, "@OH19", UNKNOWN),
    OH20(Octahedral, "@OH20", UNKNOWN),
    OH21(Octahedral, "@OH21", UNKNOWN),
    OH22(Octahedral, "@OH22", UNKNOWN),
    OH23(Octahedral, "@OH23", UNKNOWN),
    OH24(Octahedral, "@OH24", UNKNOWN),
    OH25(Octahedral, "@OH25", UNKNOWN),
    OH26(Octahedral, "@OH26", UNKNOWN),
    OH27(Octahedral, "@OH27", UNKNOWN),
    OH28(Octahedral, "@OH28", UNKNOWN),
    OH29(Octahedral, "@OH29", UNKNOWN),
    OH30(Octahedral, "@OH30", UNKNOWN);

    /** Type of configuration. */
    private final Type type;

    /** Symbol used to represent configuration */
    private final String symbol;

    /** Shorthand - often converted to this in output */
    private final Configuration shorthand;

    /** Lookup tables for trigonal bipyramidal and octahedral */
    private static final Configuration[] tbs = new Configuration[21];
    private static final Configuration[] ohs = new Configuration[31];

    // initialise trigonal lookup
    static {
        int i = 1;
        for (Configuration config : values()) {
            if (config.type().equals(TrigonalBipyramidal))
                tbs[i++] = config;
        }
    }

    // initialise octahedral lookup
    static {
        int i = 1;
        for (Configuration config : values()) {
            if (config.type().equals(Octahedral))
                ohs[i++] = config;
        }
    }

    private Configuration(Type type, String symbol, Configuration shorthand) {
        this.type = type;
        this.symbol = symbol;
        this.shorthand = shorthand;
    }

    private Configuration(Type type, String symbol) {
        this.type = type;
        this.symbol = symbol;
        this.shorthand = this;
    }

    /**
     * Access the shorthand for the configuration, if no shorthand is defined
     * {@link #UNKNOWN} is returned.
     *
     * @return the shorthand '@' or '@@'
     */
    public Configuration shorthand() {
        return shorthand;
    }

    /**
     * Symbol of the chiral configuration.
     *
     * @return the symbol
     */
    public String symbol() {
        return symbol;
    }

    /**
     * The general type of relative configuration this represents.
     *
     * @return type of the configuration
     * @see Type
     */
    public Type type() {
        return type;
    }

    /**
     * Read a chiral configuration from a character buffer and progress the
     * buffer. If there is no configuration then {@link Configuration#UNKNOWN}
     * is returned. Encountering an invalid permutation designator (e.g.
     * &#64;TB21) or incomplete class (e.g. &#64;T) will throw an invalid smiles
     * exception.
     *
     * @param buffer a character buffer
     * @return the configuration
     * @throws InvalidSmilesException
     */
    static Configuration read(final CharBuffer buffer) throws
                                                       InvalidSmilesException {
        if (buffer.getIf('@')) {
            if (buffer.getIf('@')) {
                return Configuration.CLOCKWISE;
            } else if (buffer.getIf('T')) {
                // TH (tetrahedral) or TB (trigonal bipyramidal)
                if (buffer.getIf('H')) {
                    if (buffer.getIf('1'))
                        return Configuration.TH1;
                    else if (buffer.getIf('2'))
                        return Configuration.TH2;
                    else
                        throw new InvalidSmilesException("invalid permutation designator for @TH, valid values are @TH1 or @TH2:",
                                                         buffer);
                } else if (buffer.getIf('B')) {
                    int num = buffer.getNumber();
                    if (num < 1 || num > 20)
                        throw new InvalidSmilesException("invalid permutation designator for @TB, valid values are '@TB1, @TB2, ... @TB20:'",
                                                         buffer);
                    return tbs[num];
                }
                throw new InvalidSmilesException("'@T' is not a valid chiral specification:", buffer);
            } else if (buffer.getIf('A')) {
                // allene (extended tetrahedral)
                if (buffer.getIf('L')) {
                    if (buffer.getIf('1'))
                        return Configuration.AL1;
                    else if (buffer.getIf('2'))
                        return Configuration.AL2;
                    else
                        throw new InvalidSmilesException("invalid permutation designator for @AL, valid values are '@AL1 or @AL2':", buffer);
                } else {
                    throw new InvalidSmilesException("'@A' is not a valid chiral specification:", buffer);
                }
            } else if (buffer.getIf('S')) {
                // square planar
                if (buffer.getIf('P')) {
                    if (buffer.getIf('1'))
                        return Configuration.SP1;
                    else if (buffer.getIf('2'))
                        return Configuration.SP2;
                    else if (buffer.getIf('3'))
                        return Configuration.SP3;
                    else
                        throw new InvalidSmilesException("invalid permutation designator for @SP, valid values are '@SP1, @SP2 or @SP3':",
                                                         buffer);
                } else {
                    throw new InvalidSmilesException("'@S' is not a valid chiral specification:", buffer);
                }
            } else if (buffer.getIf('O')) {
                if (buffer.getIf('H')) {
                    // octahedral
                    int num = buffer.getNumber();
                    if (num < 1 || num > 30)
                        throw new InvalidSmilesException("invalid permutation designator for @OH, valud values are '@OH1, @OH2, ... @OH30':", buffer);
                    return ohs[num];
                } else {
                    throw new InvalidSmilesException("'@O' is not a valid chiral specification:", buffer);
                }
            } else {
                return Configuration.ANTI_CLOCKWISE;
            }
        }
        return Configuration.UNKNOWN;
    }

    /** Types of configuration. */
    public enum Type {
        None,
        Implicit,
        Tetrahedral,
        ExtendedTetrahedral,
        SquarePlanar,
        TrigonalBipyramidal,
        Octahedral
    }
}

