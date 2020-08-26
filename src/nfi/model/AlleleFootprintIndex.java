
package nfi.model;

import jam.app.JamProperties;
import jam.math.DoubleUtil;

import jean.hla.Allele;
import jean.neo.PeptidePair;
import jean.peptide.Peptide;

import pepmhc.affy.AffinityCache;
import pepmhc.affy.AffinityMethod;
import pepmhc.stab.StabilityCache;
import pepmhc.stab.StabilityMethod;

/**
 * Defines an interface to calculate the neo-peptide footprint index
 * for single HLA alleles.
 */
public abstract class AlleleFootprintIndex {
    private static AlleleFootprintIndex global;

    /**
     * Name of the system property that defines the global allele
     * footprint index type.
     */
    public static final String TYPE_PROPERTY = "nfi.model.alleleFootprintType";

    /**
     * Returns the log-affinity footprint index.
     */
    public static final AlleleFootprintIndex LOG_AFFINITY = new LogAffinity();

    /**
     * Returns the log-stability footprint index.
     */
    public static final AlleleFootprintIndex LOG_STABILITY = new LogStability();

    /**
     * Returns the global allele footprint calculator with the type
     * specified by the {@code nfi.model.alleleFootprintType} system
     * property.
     *
     * @return the global allele footprint calculator with the type
     * specified by the {@code nfi.model.alleleFootprintType} system
     * property.
     *
     * @throws RuntimeException unless the system property is defined.
     */
    public static AlleleFootprintIndex global() {
        if (global == null)
            global = resolveGlobalIndex();

        return global;
    }

    private static AlleleFootprintIndex resolveGlobalIndex() {
        return resolveGlobalType().getAlleleFootprintIndex();
    }

    private static AlleleFootprintType resolveGlobalType() {
        return JamProperties.getRequiredEnum(TYPE_PROPERTY, AlleleFootprintType.class);
    }

    /**
     * Computes the neo-peptide footprint index for a single HLA
     * allele.
     *
     * @param pair the neo/self peptide pair.
     *
     * @param allele an HLA allele from a patient genotype.
     *
     * @return the neo-peptide footprint index for the peptides and
     * allele.
     */
    public abstract double compute(PeptidePair pair, Allele allele);

    /**
     * Returns the enumerated calculation type for this footprint.
     *
     * @return the enumerated calculation type for this footprint.
     */
    public abstract AlleleFootprintType getType();

    // -----------------------------------------------------------------

    private static final class LogAffinity extends AlleleFootprintIndex {
        @Override public double compute(PeptidePair pair, Allele allele) {
            //
            // Affinity is expressed as an IC50 concentration:
            // peptides with a lower IC50 bind more strongly, so we
            // invert the ratio relative to the stability model...
            //
            AffinityMethod method = AffinityMethod.NET_MHC_PAN;
            AffinityCache  cache  = AffinityCache.instance(method, allele);

            Peptide neoPeptide = pair.neo();
            Peptide selfPeptide = pair.self();

            double neoAffinity = cache.require(neoPeptide).getAffinity();
            double selfAffinity = cache.require(selfPeptide).getAffinity();

            return DoubleUtil.log2(selfAffinity / neoAffinity);
        }

        @Override public AlleleFootprintType getType() {
            return AlleleFootprintType.LOG_AFFINITY;
        }
    }

    // -----------------------------------------------------------------

    private static final class LogStability extends AlleleFootprintIndex {
        @Override public double compute(PeptidePair pair, Allele allele) {
            StabilityMethod method = StabilityMethod.NET_MHC_STAB_PAN;
            StabilityCache  cache  = StabilityCache.instance(method, allele);

            Peptide neoPeptide = pair.neo();
            Peptide selfPeptide = pair.self();

            double neoHalfLife = cache.require(neoPeptide).getHalfLife();
            double selfHalfLife = cache.require(selfPeptide).getHalfLife();

            return DoubleUtil.log2(neoHalfLife / selfHalfLife);
        }

        @Override public AlleleFootprintType getType() {
            return AlleleFootprintType.LOG_STABILITY;
        }
    }
}
