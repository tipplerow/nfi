
package nfi.model;

/**
 * Enumerates the calculation types for the single-allele neo-peptide
 * footprint index.
 */
public enum AlleleFootprintType {
    /**
     * The footprint index is defined as {@code log2(a_SELF / a_NEO)},
     * where {@code a_SELF} and {@code a_NEO} are affinities for the
     * self-peptide and neo-peptide, respectively, expressed as IC50
     * concentrations.  The self-peptide affinity is in the numerator
     * because IC50 concentrations and peptide-MHC half-lifes have an
     * inverse relationship.
     */
    LOG_AFFINITY {
        @Override public AlleleFootprintIndex getAlleleFootprintIndex() {
            return AlleleFootprintIndex.LOG_AFFINITY;
        }
    },

    /**
     * The footprint index is defined as {@code log2(h_NEO / h_SELF)},
     * where {@code h_NEO} and {@code h_SELF} are the stability of the
     * neo-peptide and self-peptide, respectively, given as half-lives
     * of the peptide-MHC complex.
     */
    LOG_STABILITY {
        @Override public AlleleFootprintIndex getAlleleFootprintIndex() {
            return AlleleFootprintIndex.LOG_STABILITY;
        }
    };

    /**
     * Returns the allele footprint index calculator for this type.
     *
     * @return the allele footprint index calculator for this type.
     */
    public abstract AlleleFootprintIndex getAlleleFootprintIndex();
}

