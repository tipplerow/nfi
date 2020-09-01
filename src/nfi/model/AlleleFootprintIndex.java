
package nfi.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jam.app.JamProperties;
import jam.math.DoubleUtil;

import jene.hla.Allele;
import jene.neo.PeptidePair;
import jene.peptide.Peptide;

import pepmhc.bind.BindPredictor;
import pepmhc.bind.BindRecord;
import pepmhc.bind.BindRecordMap;

import pepmhc.affy.net.NetMHCPanPredictor;
import pepmhc.stab.net.NetStabPredictor;

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
     * Computes neo-peptide footprint indexes for a single HLA
     * allele and neo/self peptide pair.
     *
     * @param allele the HLA allele of interest
     *
     * @param pair the neo/self peptide pair of interest.
     *
     * @return the footprint index record for the given allele and
     * neo/self peptide pair.
     */
    public AlleleFootprintRecord compute(Allele allele, PeptidePair pair) {
        return compute(allele, List.of(pair)).get(0);
    }

    /**
     * Computes neo-peptide footprint indexes for a single HLA
     * allele and a collection of neo/self peptide pairs.
     *
     * @param allele the HLA allele of interest
     *
     * @param pairs the neo/self peptide pairs of interest.
     *
     * @return a list containing the footprint index records for the
     * given allele and all neo/self peptide pairs.
     */
    public List<AlleleFootprintRecord> compute(Allele allele, Collection<PeptidePair> pairs) {
        //
        // It is more efficient to compute all binding records in a
        // single call to the underlying engine...
        //
        BindRecordMap bind = mapBinding(allele, pairs);
        AlleleFootprintType type = getFootprintType();

        List<AlleleFootprintRecord> records =
            new ArrayList<AlleleFootprintRecord>(pairs.size());

        for (PeptidePair pair : pairs) {
            double index =
                compute(allele, pair, bind);

            AlleleFootprintRecord record =
                AlleleFootprintRecord.create(allele, pair, type, index);

            records.add(record);
        }

        return records;
    }

    @SuppressWarnings("unchecked")
    private BindRecordMap mapBinding(Allele allele, Collection<PeptidePair> pairs) {
        return getBindPredictor().map(allele, PeptidePair.peptides(pairs));
    }

    /**
     * Computes neo-peptide footprint indexes for collections of HLA
     * alleles and neo/self peptide pairs.
     *
     * @param alleles the HLA alleles of interest.
     *
     * @param pairs the neo/self peptide pairs of interest.
     *
     * @return a list containing the footprint index records for all
     * allele-pair combinations.
     */
    public List<AlleleFootprintRecord> compute(Collection<Allele> alleles, Collection<PeptidePair> pairs) {
        int recordCount = alleles.size() * pairs.size();

        List<AlleleFootprintRecord> records =
            new ArrayList<AlleleFootprintRecord>(recordCount);

        for (Allele allele : alleles)
            records.addAll(compute(allele, pairs));

        return records;
    }

    /**
     * Computes the neo-peptide footprint index for a single HLA
     * allele and neo/self peptide pair.
     *
     * @param allele an HLA allele of interest.
     *
     * @param pair the neo/self peptide pair of interest.
     *
     * @param bind a map of binding records, which must contain
     * records for the neo- and self-peptides.
     *
     * @return the footprint index record for the given allele,
     * peptide pair, and binding records.
     *
     * @throws RuntimeException unless the binding map contains
     * records for the neo- and self-peptides.
     */
    public abstract double compute(Allele allele, PeptidePair pair, BindRecordMap bind);

    /**
     * Returns the enumerated calculation type for this footprint.
     *
     * @return the enumerated calculation type for this footprint.
     */
    public abstract AlleleFootprintType getFootprintType();

    /**
     * Returns the binding strength prediction method used by this
     * footprint.
     *
     * @return the binding strength prediction method used by this
     * footprint.
     */
    public abstract BindPredictor getBindPredictor();

    // -----------------------------------------------------------------

    private static final class LogAffinity extends AlleleFootprintIndex {
        @Override public double compute(Allele allele, PeptidePair pair, BindRecordMap bind) {
            //
            // Affinity is expressed as an IC50 concentration:
            // peptides with a lower IC50 bind more strongly, so we
            // invert the ratio relative to the stability model...
            //
            Peptide neoPeptide = pair.neo();
            Peptide selfPeptide = pair.self();

            double neoAffinity = bind.require(neoPeptide).getAffinity();
            double selfAffinity = bind.require(selfPeptide).getAffinity();

            return DoubleUtil.log2(selfAffinity / neoAffinity);
        }

        @Override public AlleleFootprintType getFootprintType() {
            return AlleleFootprintType.LOG_AFFINITY;
        }

        @Override public BindPredictor getBindPredictor() {
            return NetMHCPanPredictor.INSTANCE;
        }
    }

    // -----------------------------------------------------------------

    private static final class LogStability extends AlleleFootprintIndex {
        @Override public double compute(Allele allele, PeptidePair pair, BindRecordMap bind) {
            Peptide neoPeptide = pair.neo();
            Peptide selfPeptide = pair.self();

            double neoHalfLife = bind.require(neoPeptide).getHalfLife();
            double selfHalfLife = bind.require(selfPeptide).getHalfLife();

            return DoubleUtil.log2(neoHalfLife / selfHalfLife);
        }

        @Override public AlleleFootprintType getFootprintType() {
            return AlleleFootprintType.LOG_STABILITY;
        }

        @Override public BindPredictor getBindPredictor() {
            return NetStabPredictor.INSTANCE;
        }
    }
}
