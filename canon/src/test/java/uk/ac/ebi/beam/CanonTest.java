package uk.ac.ebi.beam;

import org.junit.Test;
import uk.ac.ebi.beam.CanonOld;

import java.io.IOException;

/** @author John May */
public class CanonTest {

    @Test public void benzene() throws IOException {
        System.out.println(CanonOld.ofSmi("C1CCCCC1"));       
    }
                               
    @Test public void fullerene_c70() throws IOException {
        System.out.println(CanonOld.ofSmi("C1=2C3=C4C5=C1C1=C6C7=C5C5=C8C4=C4C9=C3C3=C%10C=2C2=C1C1=C%11C%12=C%13C%14=C%15C%16=C%17C%18=C%19C%20=C%16C%16=C%14C%12=C%12C%14=C%21C%22=C(C%20=C%16%14)C%14=C%19C%16=C(C4=C8C(=C%18%16)C4=C%17C%15=C(C7=C54)C%13=C61)C1=C%14C%22=C(C3=C91)C1=C%21C%12=C%11C2=C%101"));       
    }
    
}
