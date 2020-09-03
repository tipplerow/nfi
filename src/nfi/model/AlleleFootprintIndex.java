
package nfi.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import jam.app.JamProperties;
import jam.math.DoubleUtil;

import jene.hla.Allele;
import jene.neo.PeptidePair;
import jene.neo.PeptidePairRecord;
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
     * @param pairRecord the neo/self peptide pair of interest.
     *
     * @return the footprint index record for the given allele and
     * neo/self peptide pair.
     */
    public AlleleFootprintRecord compute(Allele allele, PeptidePairRecord pairRecord) {
        return compute(allele, List.of(pairRecord)).get(0);
    }

    /**
     * Computes neo-peptide footprint indexes for a single HLA
     * allele and a collection of neo/self peptide pairs.
     *
     * @param allele the HLA allele of interest
     *
     * @param pairRecords the neo/self peptide pairs of interest.
     *
     * @return a list containing the footprint index records for the
     * given allele and all neo/self peptide pairs.
     */
    public List<AlleleFootprintRecord> compute(Allele allele, Collection<PeptidePairRecord> pairRecords) {
        //
        // It is more efficient to compute all binding records in a
        // single call to the underlying engine...
        //
        BindRecordMap bindingMap = mapBinding(allele, pairRecords);

        List<AlleleFootprintRecord> footprintRecords =
            new ArrayList<AlleleFootprintRecord>(pairRecords.size());

        for (PeptidePairRecord pairRecord : pairRecords)
            footprintRecords.add(compute(allele, pairRecord, bindingMap));

        return footprintRecords;
    }

    @SuppressWarnings("unchecked")
    private BindRecordMap mapBinding(Allele allele, Collection<PeptidePairRecord> pairRecords) {
        return getBindPredictor().map(allele, PeptidePairRecord.peptides(pairRecords));
    }

    private AlleleFootprintRecord compute(Allele patientAllele, PeptidePairRecord pairRecord, BindRecordMap bindingMap) {
        AlleleFootprintType footprintType = getFootprintType();

        Peptide neoPeptide = pairRecord.getNeoPeptide();
        Peptide selfPeptide = pairRecord.getSelfPeptide();

        BindRecord neoBindRecord = bindingMap.require(neoPeptide);
        BindRecord selfBindRecord = bindingMap.require(selfPeptide);

        double footprintIndex = compute(neoBindRecord, selfBindRecord);

        return AlleleFootprintRecord.create(pairRecord,
                                            patientAllele,
                                            footprintType,
                                            neoBindRecord,
                                            selfBindRecord,
                                            footprintIndex);
    }

    /**
     * Computes neo-peptide footprint indexes for collections of HLA
     * alleles and neo/self peptide pairs.
     *
     * @param alleles the HLA alleles of interest.
     *
     * @param pairRecords the neo/self peptide pairs of interest.
     *
     * @return a list containing the footprint index records for all
     * allele-pair combinations.
     */
    public List<AlleleFootprintRecord> compute(Collection<Allele> alleles, Collection<PeptidePairRecord> pairRecords) {
        int recordCount = alleles.size() * pairRecords.size();

        List<AlleleFootprintRecord> footprintRecords =
            new ArrayList<AlleleFootprintRecord>(recordCount);

        for (Allele allele : alleles)
            footprintRecords.addAll(compute(allele, pairRecords));

        return footprintRecords;
    }

    /**
     * Computes the footprint index for a neo/self peptide pair.
     *
     * @param neoBindRecord the neo-antigen binding record.
     *
     * @param selfBindRecord the self-antigen binding record.
     *
     * @return the footprint index for the given binding records.
     */
    public abstract double compute(BindRecord neoBindRecord, BindRecord selfBindRecord);

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
        @Override public double compute(BindRecord neoBindRecord, BindRecord selfBindRecord) {
            //
            // Affinity is expressed as an IC50 concentration:
            // peptides with a lower IC50 bind more strongly, so we
            // invert the ratio relative to the stability model...
            //
            double neoAffinity = neoBindRecord.getAffinity();
            double selfAffinity = selfBindRecord.getAffinity();

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
        @Override public double compute(BindRecord neoBindRecord, BindRecord selfBindRecord) {
            double neoHalfLife = neoBindRecord.getHalfLife();
            double selfHalfLife = selfBindRecord.getHalfLife();

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
