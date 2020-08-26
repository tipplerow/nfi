
package nfi.junit;

import jam.junit.NumericTestBase;
import jam.util.RegexUtil;

import jean.hla.Allele;
import jean.neo.NeoPeptide;
import jean.neo.PeptidePair;
import jean.neo.SelfPeptide;
import jean.peptide.Peptide;

import nfi.model.AlleleFootprintIndex;
import nfi.model.AlleleFootprintType;

import org.junit.*;
import static org.junit.Assert.*;

public class AlleleFootprintIndexTest extends NumericTestBase {
    static {
        System.setProperty(AlleleFootprintIndex.TYPE_PROPERTY, "LOG_STABILITY");
    }
    private static final Allele A0101 = Allele.instance("A0101");
    private static final Allele A0201 = Allele.instance("A0201");

    private static final NeoPeptide neo1 = NeoPeptide.instance("FQASPMHAV");
    private static final SelfPeptide self1 = SelfPeptide.instance("FLASPMHAV");
    private static final PeptidePair pair1 = PeptidePair.instance(self1, neo1);

    private static final NeoPeptide neo2 = NeoPeptide.instance("FADSPMHAL");
    private static final SelfPeptide self2 = SelfPeptide.instance("FTDSPMHAV");
    private static final PeptidePair pair2 = PeptidePair.instance(self2, neo2);

    private static final AlleleFootprintIndex logAffinity = AlleleFootprintIndex.LOG_AFFINITY;
    private static final AlleleFootprintIndex logStability = AlleleFootprintIndex.LOG_STABILITY;

    @Test public void testGlobal() {
        assertEquals(AlleleFootprintType.LOG_STABILITY, AlleleFootprintIndex.global().getFootprintType());
    }

    @Test public void testLogAffinity() {
        assertEquals(-0.622, logAffinity.compute(A0101, pair1).getIndex(), 0.001);
        assertEquals(-3.612, logAffinity.compute(A0101, pair2).getIndex(), 0.001);

        assertEquals(-1.528, logAffinity.compute(A0201, pair1).getIndex(), 0.001);
        assertEquals(-4.165, logAffinity.compute(A0201, pair2).getIndex(), 0.001);
    }

    @Test public void testLogStability() {
        assertEquals( 0.047, logStability.compute(A0101, pair1).getIndex(), 0.001);
        assertEquals(-0.396, logStability.compute(A0101, pair2).getIndex(), 0.001);

        assertEquals(-1.546, logStability.compute(A0201, pair1).getIndex(), 0.001);
        assertEquals(-1.680, logStability.compute(A0201, pair2).getIndex(), 0.001);
    }

    public static void main(String[] args) {
        org.junit.runner.JUnitCore.main("nfi.junit.AlleleFootprintIndexTest");
    }
}
