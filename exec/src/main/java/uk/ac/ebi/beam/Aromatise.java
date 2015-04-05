/*
 * Copyright (c) 2015. John May
 */

package uk.ac.ebi.beam;

import joptsimple.OptionSet;

/**
 * Simple module simply Kekulises then emits and normalised
 * (by beam's model) aromatic form of the SMILES.
 */
public final class Aromatise extends FunctorCmdLnModule {

    public Aromatise() {
        super("arom");
    }

    @Override
    Functor createFunctor(OptionSet optionSet) {
        return new Functor() {
            @Override
            String map(String str) throws InvalidSmilesException {
                return Graph.fromSmiles(str).kekule().aromatic().toSmiles() + suffixedId(str);
            }
        };
    }
}
