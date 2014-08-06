package uk.ac.ebi.beam;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

/** @author John May */
public class Benchmark {

    public static void main(String[] args) throws IOException {

        String fname = args[0];

        InputStream in = new FileInputStream(fname);
        if (fname.endsWith(".gz")) // could check GZIP magic
            in = new GZIPInputStream(in, 2048);

        int i = 0;

        BufferedReader br = new BufferedReader(new InputStreamReader(in), 2048);
        String line = null;
        long t0 = System.nanoTime();
        while ((line = br.readLine()) != null) {
            try {
                Graph g = Graph.fromSmiles(line);
                PartitionOld partition = CanonOld.label(g);
                g.permute(partition.ordering)
                 .resonate()
                 .toSmiles();
            } catch (Exception e) {

            }
        }
        long t1 = System.nanoTime();

        System.out.println(TimeUnit.NANOSECONDS.toMillis(t1 - t0) + " ms");

        br.close();

    }

}
