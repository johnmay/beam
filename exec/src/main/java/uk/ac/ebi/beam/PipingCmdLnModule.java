/*
 * Copyright (c) 2015. John May
 */

package uk.ac.ebi.beam;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

import java.io.*;
import java.nio.charset.Charset;
import java.util.List;

/**
 * An abstract module providing much of the boiler plate to implement a
 * simple command line module that consumes from one file (or stream)
 * out produces another.
 * <p/>
 * To use the module simply extend it and implement the
 * {@link #process(BufferedReader, BufferedWriter, OptionSet)} method.
 */
public abstract class PipingCmdLnModule implements CmdLnModule {

    private final Charset UTF_8 = Charset.forName("UTF-8");

    private final String name;
    protected final OptionParser optparser = new OptionParser();

    public PipingCmdLnModule(String name) {
        this.name = name;
        this.optparser.accepts("prog-off", "no progress indicator");
    }

    /**
     * @inheritDoc
     */
    @Override
    public final String name() {
        return name;
    }

    /**
     * @inheritDoc
     */
    @Override
    public final String getHelpInfo() {
        StringWriter sw = new StringWriter();
        sw.append("\n").append(name).append(":\n");
        sw.append("\tusage: beam ").append(name).append("{OPTS}").append(" [{in.smi}] [{out.smi}]\n");
        try {
            optparser.printHelpOn(sw);
        } catch (IOException e) {
            throw new InternalError(e);
        }
        return sw.toString();
    }

    /**
     * @inheritDoc
     */
    @Override
    public final void exec(String[] args) {
        try {
            process(args);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // simply to make the try/catch in exec simpler
    private void process(String[] args) throws IOException {

        OptionSet optset = optparser.parse(args);
        List<?> nonopt = optset.nonOptionArguments();

        final InputStream in = nonopt.size() < 1 ? System.in
                                                 : new FileInputStream(nonopt.get(0).toString());
        final OutputStream out = nonopt.size() < 2 ? System.out
                                                   : new FileOutputStream(nonopt.get(1).toString());

        try (BufferedWriter bwtr = new BufferedWriter(new OutputStreamWriter(out, UTF_8));
             BufferedReader brdr = new BufferedReader(new InputStreamReader(in, UTF_8))) {
            process(brdr, bwtr, optset);
        }
    }

    /**
     * Consume input from a buffered reader and produce something in an output stream (usually line-by-line).
     *
     * @param brdr   input reader (UTF-8)
     * @param bwtr   output reader (UTF-8)
     * @param optset options for the module
     * @throws IOException low-level error, SMILES syntax errors etc
     *                     should normally be skipped and usually reported
     */
    abstract void process(final BufferedReader brdr, final BufferedWriter bwtr, OptionSet optset) throws IOException;

    /**
     * Reports a message to standard error prefixing the module name. The
     * syntax is essentially printf, note no new line.
     *
     * @param str format string
     * @param args the arguments
     */
    protected void report(String str, Object... args) {
        System.err.printf("\r[beam:" + name + "] " + str, args);
    }

    /**
     * Access the ID part of the a SMILES line. The delimiter is included allowing
     * direct concatenation to an output SMILES.
     *
     * @param smi input SMILES
     * @return the id, or empty string
     */
    protected static String suffixedId(String smi) {
        for (int i = 0; i < smi.length(); i++) {
            if (smi.charAt(i) == ' ' || smi.charAt(i) == '\t')
                return smi.substring(i);
        }
        return "";
    }
}
