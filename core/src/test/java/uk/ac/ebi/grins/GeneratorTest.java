package uk.ac.ebi.grins;

import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.util.Random;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

/** @author John May */
public class GeneratorTest {

    @Test public void permuteTH_3_nonRing() throws Exception {
        String input = "C[C@H](N)O";
        ChemicalGraph g = Parser.parse(input);
        assertThat(Generator.generate(g), is(input));
    }

    @Test public void permuteTH_4_nonRing() throws Exception {
        String input = "C[C@]([H])(N)O";
        ChemicalGraph g = Parser.parse(input);
        assertThat(Generator.generate(g), is(input));
    }

    @Test public void permuteTH_4_ring() throws Exception {
        String input = "C[C@]12CCCC[C@@]1(C)OCCC2";
        ChemicalGraph g = Parser.parse(input);
        assertThat(Generator.generate(g), is(input));
    }

    @Ignore public void test() throws InvalidSmilesException {
        System.out.println(randomPermutations("[C@]([H])(N)(C)C(=O)O", 50));
    }

    @Ignore public void test2() throws InvalidSmilesException {
        System.out.println(randomPermutations("[C@H](N)(C)C(=O)O", 50));
    }

    @Ignore public void test3() throws InvalidSmilesException {
        System.out.println(randomPermutations("[C@H]12CCCC[C@@]1(C)OCCC2", 50));
    }

    @Test public void implicitHCentre() throws InvalidSmilesException {

        roundTrip("[C@@H](N)(O)C");

        // permutations
        roundTrip("[C@@H](N)(O)C", new int[]{0, 1, 2, 3}, "[C@@H](N)(O)C");
        roundTrip("[C@@H](N)(O)C", new int[]{0, 1, 3, 2}, "[C@H](N)(C)O");
        roundTrip("[C@@H](N)(O)C", new int[]{0, 2, 1, 3}, "[C@H](O)(N)C");
        roundTrip("[C@@H](N)(O)C", new int[]{0, 2, 3, 1}, "[C@@H](C)(N)O");
        roundTrip("[C@@H](N)(O)C", new int[]{0, 3, 1, 2}, "[C@@H](O)(C)N");
        roundTrip("[C@@H](N)(O)C", new int[]{0, 3, 2, 1}, "[C@H](C)(O)N");

        roundTrip("[C@@H](N)(O)C", new int[]{1, 0, 2, 3}, "N[C@H](O)C");
        roundTrip("[C@@H](N)(O)C", new int[]{1, 0, 3, 2}, "N[C@@H](C)O");

        roundTrip("[C@@H](N)(O)C", new int[]{1, 2, 0, 3}, "O[C@@H](N)C");
        roundTrip("[C@@H](N)(O)C", new int[]{1, 3, 0, 2}, "O[C@H](C)N");

        roundTrip("[C@@H](N)(O)C", new int[]{1, 2, 3, 0}, "C[C@H](N)O");
        roundTrip("[C@@H](N)(O)C", new int[]{1, 3, 2, 0}, "C[C@@H](O)N");

        roundTrip("[C@H](N)(C)O");

        roundTrip("[C@H](N)(C)O", new int[]{0, 1, 2, 3}, "[C@H](N)(C)O");
        roundTrip("[C@H](N)(C)O", new int[]{0, 1, 3, 2}, "[C@@H](N)(O)C");
        roundTrip("[C@H](N)(C)O", new int[]{0, 2, 1, 3}, "[C@@H](C)(N)O");
        roundTrip("[C@H](N)(C)O", new int[]{0, 2, 3, 1}, "[C@H](O)(N)C");
        roundTrip("[C@H](N)(C)O", new int[]{0, 3, 1, 2}, "[C@H](C)(O)N");
        roundTrip("[C@H](N)(C)O", new int[]{0, 3, 2, 1}, "[C@@H](O)(C)N");

        roundTrip("[C@H](N)(C)O", new int[]{1, 0, 2, 3}, "N[C@@H](C)O");
        roundTrip("[C@H](N)(C)O", new int[]{1, 0, 3, 2}, "N[C@H](O)C");

        roundTrip("[C@H](N)(C)O", new int[]{1, 2, 0, 3}, "C[C@H](N)O");
        roundTrip("[C@H](N)(C)O", new int[]{1, 3, 0, 2}, "C[C@@H](O)N");

        roundTrip("[C@H](N)(C)O", new int[]{1, 2, 3, 0}, "O[C@@H](N)C");
        roundTrip("[C@H](N)(C)O", new int[]{1, 3, 2, 0}, "O[C@H](C)N");

        roundTrip("N[C@@H](C)O");
        roundTrip("N[C@@H](C)O");
        roundTrip("N[C@H](O)C");
        roundTrip("O[C@@H](N)C");
        roundTrip("O[C@H](C)N");
        roundTrip("C[C@@H](O)N");
        roundTrip("C[C@H](N)O");
    }

    @Test public void ring_closures1() throws Exception {
        roundTrip("C1=CN=CC2=NC=N[C@@H]21");
    }

    @Test public void ring_closures2() throws Exception {
        roundTrip("C1=CN=CC2=NC=N[C@H]21");
    }

    @Test public void ring_closures3() throws Exception {
        roundTrip("C1=CC(=CC2=NC(=N[C@@H]21)C(F)(F)F)N");
    }

    @Test public void ring_closures4() throws Exception {
        roundTrip("C1=CC(=CC2=NC(=N[C@H]21)C(F)(F)F)N");
    }


    @Test public void lowRingNumberOrder() throws InvalidSmilesException {
        roundTrip("C1=CC2=CC=CC=C2C=C1");
    }

    @Test public void multipleRingNumberOrder() throws InvalidSmilesException {
        roundTrip("C1=CC2=C3C4=C5C(C=CC6=C5C7=C(C=C6)C=CC(C=C2)=C37)=CC=C14");
    }

    @Test public void highRingNumberOrder() throws InvalidSmilesException {
        roundTrip("C1CC2CCC3=C4C2=C5C1CCC6=C5C7=C8C(C=C9CCC%10CCC%11CCC%12=CC(=C3)C(C%13=C8C9=C%10C%11=C%12%13)=C47)=C6");
    }

    @Test public void bondTypeOnFirstAtom1() throws InvalidSmilesException {
        String smi = "C1C=CC=CC=1";
        String exp = "C=1C=CC=CC1";
        assertThat(Generator.generate(Parser.parse(smi)), is(exp));
    }

    @Test public void bondTypeOnFirstAtom2() throws InvalidSmilesException {
        String smi = "C=1C=CC=CC1";
        String exp = "C=1C=CC=CC1";
        assertThat(Generator.generate(Parser.parse(smi)), is(exp));
    }

    @Test public void bondTypeOnFirstAtom3() throws InvalidSmilesException {
        String smi = "C=1C=CC=CC=1";
        String exp = "C=1C=CC=CC1";
        assertThat(Generator.generate(Parser.parse(smi)), is(exp));
    }

    @Test public void directionalBondTypeOnFirstAtom1() throws
                                                        InvalidSmilesException {
        String smi = "C1CCCCCCCCCCC\\C=C/1";
        String exp = "C\\1CCCCCCCCCCC\\C=C1";
        assertThat(Generator.generate(Parser.parse(smi)), is(exp));
    }

    @Test public void directionalBondTypeOnFirstAtom2() throws
                                                        InvalidSmilesException {
        String smi = "C\\1CCCCCCCCCCC\\C=C1";
        String exp = "C\\1CCCCCCCCCCC\\C=C1";
        assertThat(Generator.generate(Parser.parse(smi)), is(exp));
    }

    @Test public void directionalBondTypeOnFirstAtom3() throws
                                                        InvalidSmilesException {
        String smi = "C\\1CCCCCCCCCCC\\C=C/1";
        String exp = "C\\1CCCCCCCCCCC\\C=C1";
        assertThat(Generator.generate(Parser.parse(smi)), is(exp));
    }

    @Test public void reuseNumbering() throws IOException {
        Generator generator = new Generator(ChemicalGraph.fromSmiles("c1cc1c2ccc2"),
                                            new Generator.ReuseRingNumbering(1));
        assertThat(generator.string(), is("c1cc1c1ccc1"));
    }

    @Test public void sodiumChloride() throws InvalidSmilesException {
        roundTrip("[Na+].[Cl-]");
    }

    @Test public void disconnected() throws InvalidSmilesException {
        roundTrip("CCCC.OOOO.C[CH]C.CNO");
    }

    @Test public void reusingNumbering() {
        Generator.RingNumbering rnums = new Generator.ReuseRingNumbering(0);
        for (int i = 0; i < 50; i++) {
            int rnum = rnums.next();
            assertThat(rnum, is(i));
            rnums.use(rnum);
        }
        rnums.free(40);
        rnums.free(20);
        rnums.free(4);
        assertThat(rnums.next(), is(4));
        rnums.use(4);
        assertThat(rnums.next(), is(20));
        rnums.use(20);
        assertThat(rnums.next(), is(40));
        rnums.use(40);
        for (int i = 50; i < 100; i++) {
            int rnum = rnums.next();
            assertThat(rnum, is(i));
            rnums.use(rnum);
        }
    }

    @Test public void iterativeNumbering() {
        Generator.RingNumbering rnums = new Generator.IterativeRingNumbering(0);
        for (int i = 0; i < 50; i++) {
            int rnum = rnums.next();
            assertThat(rnum, is(i));
            rnums.use(rnum);
        }
        rnums.free(40);
        rnums.free(25);
        assertThat(rnums.next(), is(50));
        rnums.use(50);
        assertThat(rnums.next(), is(51));
        rnums.use(51);
        assertThat(rnums.next(), is(52));
        rnums.use(52);
        for (int i = 53; i < 100; i++) {
            int rnum = rnums.next();
            assertThat(rnum, is(i));
            rnums.use(rnum);
        }
        rnums.free(20);
        rnums.free(5);
        assertThat(rnums.next(), is(5));
        rnums.use(5);
        assertThat(rnums.next(), is(20));
        rnums.use(20);
        assertThat(rnums.next(), is(25));
        rnums.use(25);
        assertThat(rnums.next(), is(40));
        rnums.use(40);
    }

    @Test(expected = IllegalArgumentException.class)
    public void maxRingNumbers() {
        Generator.RingNumbering rnums = new Generator.IterativeRingNumbering(0);
        for (int i = 0; i < 101; i++) {
            int rnum = rnums.next();
            rnums.use(rnum);
        }
    }

    static void roundTrip(String smi) throws InvalidSmilesException {
        assertThat(Generator.generate(Parser.parse(smi)), is(smi));
    }

    static void roundTrip(String smi, int[] p, String res) throws
                                                           InvalidSmilesException {
        assertThat(Generator.generate(Parser.parse(smi).permute(p)), is(res));
    }

    /**
     * Generate random permutations of the molecule.
     *
     * @param input input SMILES
     * @param n     number of generations (how many molecules to produce)
     * @return a single SMILES string of disconnected molecules (input) randomly
     *         permuted
     * @throws InvalidSmilesException the input SMILES was invalid
     */
    private static String randomPermutations(String input, int n) throws
                                                                  InvalidSmilesException {
        ChemicalGraph g = Parser.parse(input);
        StringBuilder sb = new StringBuilder();
        sb.append(Generator.generate(g));
        for (int i = 0; i < n; i++) {
            sb.append('.');
            int[] p = random(g.order());
            String smi = Generator.generate(g.permute(p));
            g = Parser.parse(smi);
            sb.append(smi);
        }
        return sb.toString();
    }

    static int[] ident(int n) {
        int[] p = new int[n];
        for (int i = 0; i < n; i++)
            p[i] = i;
        return p;
    }

    static int[] random(int n) {
        int[] p = ident(n);
        Random rnd = new Random();
        for (int i = n; i > 1; i--)
            swap(p, i - 1, rnd.nextInt(i));
        return p;
    }

    static int[] inv(int[] p) {
        int[] q = p.clone();
        for (int i = 0; i < p.length; i++)
            q[p[i]] = i;
        return q;
    }

    static void swap(int[] p, int i, int j) {
        int tmp = p[i];
        p[i] = p[j];
        p[j] = tmp;
    }

}
