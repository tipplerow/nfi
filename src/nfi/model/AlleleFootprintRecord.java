
package nfi.model;

import jene.hla.Allele;
import jene.neo.PeptidePair;
import jene.peptide.Peptide;

/**
 * Encapsulates the results of a footprint index calculation for an
 * HLA allele and peptide pair.
 */
public final class AlleleFootprintRecord {
    private final double index;
    private final Allele allele;
    private final PeptidePair pair;
    private final AlleleFootprintType type;

    private AlleleFootprintRecord(Allele allele, PeptidePair pair, AlleleFootprintType type, double index) {
        this.allele = allele;
        this.pair   = pair;
        this.type   = type;
        this.index  = index;
    }

    /**
     * Creates a new footprint record with fixed attributes.
     *
     * @param allele the HLA allele in the calculation.
     *
     * @param pair the neo/self peptide pair in the calculation.
     *
     * @param type the enumerated footprint calculation type.
     *
     * @param index the calculated neo-antigen footprint index.
     *
     * @return a new footprint record with the specified attributes.
     */
    public static AlleleFootprintRecord create(Allele allele, PeptidePair pair, AlleleFootprintType type, double index) {
        return new AlleleFootprintRecord(allele, pair, type, index);
    }

    /**
     * Returns the HLA allele used in the calculation.
     *
     * @return the HLA allele used in the calculation.
     */
    public Allele getAllele() {
        return allele;
    }

    /**
     * Returns the neo-peptide used in the calculation.
     *
     * @return the neo-peptide used in the calculation.
     */
    public Peptide getNeoPeptide() {
        return pair.neo();
    }

    /**
     * Returns the self-peptide used in the calculation.
     *
     * @return the self-peptide used in the calculation.
     */
    public Peptide getSelfPeptide() {
        return pair.self();
    }

    /**
     * Returns the enumerated footprint calculation type.
     *
     * @return the enumerated footprint calculation type.
     */
    public AlleleFootprintType getType() {
        return type;
    }

    /**
     * Returns the computed footprint index.
     *
     * @return the computed footprint index.
     */
    public double getIndex() {
        return index;
    }
}
