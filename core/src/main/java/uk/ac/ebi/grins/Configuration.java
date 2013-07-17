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
     * The general type of relative configuration this represents.
     *
     * @return type of the configuration
     * @see Type
     */
    public Type type() {
        return type;
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

