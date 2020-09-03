
package nfi.model;

import jam.io.Delimiter;
import jam.report.LineBuilder;

import jene.hla.Allele;
import jene.neo.PeptidePairRecord;

import pepmhc.bind.BindRecord;

/**
 * Encapsulates the results of a footprint index calculation for an
 * HLA allele and peptide pair.
 */
public final class AlleleFootprintRecord {
    private final PeptidePairRecord pairRecord;

    private final Allele patientAllele;
    private final AlleleFootprintType footprintType;

    private final double neoBindingQty;
    private final double neoBindingPct;
    private final double selfBindingQty;
    private final double selfBindingPct;
    private final double footprintIndex;

    private AlleleFootprintRecord(PeptidePairRecord   pairRecord,
                                  Allele              patientAllele,
                                  AlleleFootprintType footprintType,
                                  double              neoBindingQty,
                                  double              neoBindingPct,
                                  double              selfBindingQty,
                                  double              selfBindingPct,
                                  double              footprintIndex) {
        this.pairRecord = pairRecord;

        this.patientAllele = patientAllele;
        this.footprintType = footprintType;

        this.neoBindingQty = neoBindingQty;
        this.neoBindingPct = neoBindingPct;
        this.selfBindingQty = selfBindingQty;
        this.selfBindingPct = selfBindingPct;
        this.footprintIndex = footprintIndex;
    }

    /**
     * The standard delimiter for flat files containing peptide pair
     * records.
     */
    public static final Delimiter DELIM = Delimiter.TAB;

    /**
     * Creates a new footprint record with fixed attributes.
     *
     * @param pairRecord the target of the calculation.
     *
     * @param patientAllele the HLA allele in the calculation.
     *
     * @param footprintType the enumerated footprint calculation type.
     *
     * @param neoBindRecord the neo-antigen binding record.
     *
     * @param selfBindRecord the self-antigen binding record.
     *
     * @param footprintIndex the calculated neo-antigen footprint index.
     *
     * @return a new footprint record with the specified attributes.
     */
    public static AlleleFootprintRecord create(PeptidePairRecord   pairRecord,
                                               Allele              patientAllele,
                                               AlleleFootprintType footprintType,
                                               BindRecord          neoBindRecord,
                                               BindRecord          selfBindRecord,
                                               double              footprintIndex) {
        return create(pairRecord,
                      patientAllele,
                      footprintType,
                      neoBindRecord.getStrength(),
                      neoBindRecord.getPercentile(),
                      selfBindRecord.getStrength(),
                      selfBindRecord.getPercentile(),
                      footprintIndex);
    }

    /**
     * Creates a new footprint record with fixed attributes.
     *
     * @param pairRecord the target of the calculation.
     *
     * @param patientAllele the HLA allele in the calculation.
     *
     * @param footprintType the enumerated footprint calculation type.
     *
     * @param neoBindingQty the relevant neo-antigen binding strength
     * (affinity or stability).
     *
     * @param neoBindingPct the percentile ranking of the neo-antigen
     * binding strength.
     *
     * @param selfBindingQty the relevant self-antigen binding strength
     * (affinity or stability).
     *
     * @param selfBindingPct the percentile ranking of the self-antigen
     * binding strength.
     *
     * @param footprintIndex the calculated neo-antigen footprint index.
     *
     * @return a new footprint record with the specified attributes.
     */
    public static AlleleFootprintRecord create(PeptidePairRecord   pairRecord,
                                               Allele              patientAllele,
                                               AlleleFootprintType footprintType,
                                               double              neoBindingQty,
                                               double              neoBindingPct,
                                               double              selfBindingQty,
                                               double              selfBindingPct,
                                               double              footprintIndex) {
        return new AlleleFootprintRecord(pairRecord,
                                         patientAllele,
                                         footprintType,
                                         neoBindingQty,
                                         neoBindingPct,
                                         selfBindingQty,
                                         selfBindingPct,
                                         footprintIndex);
    }

    /**
     * Returns the header line for flat files containing allele
     * footprint records.
     *
     * @return the header line for flat files containing allele
     * footprint records.
     */
    public static String header() {
        LineBuilder builder = new LineBuilder(DELIM);

        builder.append(PeptidePairRecord.header(DELIM));
        builder.append("Patient_Allele");
        builder.append("Footprint_Type");
        builder.append("Neo_Binding_Qty");
        builder.append("Neo_Binding_Pct");
        builder.append("Self_Binding_Qty");
        builder.append("Self_Binding_Pct");
        builder.append("Footprint_Index");

        return builder.toString();
    }

    /**
     * Creates a new allele footprint record by parsing a delimited
     * line from a flat file.
     *
     * @param line the line to parse.
     *
     * @return the allele footprint record encoded in the specified
     * line.
     *
     * @throws RuntimeException unless the line contains a properly
     * formatted allele footprint record.
     */
    public static AlleleFootprintRecord parse(String line) {
        String[] fields = DELIM.split(line, 13);

        PeptidePairRecord pairRecord =
            PeptidePairRecord.parse(fields, 0);

        Allele patientAllele = Allele.instance(fields[6]);
        AlleleFootprintType footprintType = AlleleFootprintType.valueOf(fields[7]);

        double neoBindingQty = Double.parseDouble(fields[8]);
        double neoBindingPct = Double.parseDouble(fields[9]);
        double selfBindingQty = Double.parseDouble(fields[10]);
        double selfBindingPct = Double.parseDouble(fields[11]);
        double footprintIndex = Double.parseDouble(fields[12]);

        return create(pairRecord,
                      patientAllele,
                      footprintType,
                      neoBindingQty,
                      neoBindingPct,
                      selfBindingQty,
                      selfBindingPct,
                      footprintIndex);
    }

    /**
     * Formats this record for output to a delimited flat file.
     *
     * @return a string containing the formatted text.
     */
    public String format() {
        LineBuilder builder = new LineBuilder(DELIM);

        builder.append(pairRecord.format(DELIM));
        builder.append(patientAllele.shortKey());
        builder.append(footprintType.name());
        builder.append(neoBindingQty, "%.2f");
        builder.append(neoBindingPct, "%.2f");
        builder.append(selfBindingQty, "%.2f");
        builder.append(selfBindingPct, "%.2f");
        builder.append(footprintIndex, "%.4f");

        return builder.toString();
    }

    /**
     * Returns the computed neo-antigen footprint index.
     *
     * @return the computed neo-antigen footprint index.
     */
    public double getFootprintIndex() {
        return footprintIndex;
    }

    /**
     * Returns the enumerated footprint calculation type.
     *
     * @return the enumerated footprint calculation type.
     */
    public AlleleFootprintType getFootprintType() {
        return footprintType;
    }

    /**
     * Returns the percentile rank of the neo-antigen binding
     * strength.
     *
     * @return the percentile rank of the neo-antigen binding
     * strength.
     */
    public double getNeoBindingPct() {
        return neoBindingPct;
    }

    /**
     * Returns the relevant neo-antigen binding quantity (binding
     * affinity as an IC50 concentration or the half-life of the
     * peptide-MHC complex).
     *
     * @return the relevant neo-antigen binding quantity.
     */
    public double getNeoBindingQty() {
        return neoBindingQty;
    }

    /**
     * Returns the HLA allele used in the calculation.
     *
     * @return the HLA allele used in the calculation.
     */
    public Allele getPatientAllele() {
        return patientAllele;
    }

    /**
     * Returns the neo/self peptide target of the calculation.
     *
     * @return the neo/self peptide target of the calculation.
     */
    public PeptidePairRecord getPeptidePairRecord() {
        return pairRecord;
    }

    /**
     * Returns the percentile rank of the self-antigen binding
     * strength.
     *
     * @return the percentile rank of the self-antigen binding
     * strength.
     */
    public double getSelfBindingPct() {
        return selfBindingPct;
    }

    /**
     * Returns the relevant self-antigen binding quantity (binding
     * affinity as an IC50 concentration or the half-life of the
     * peptide-MHC complex).
     *
     * @return the relevant self-antigen binding quantity.
     */
    public double getSelfBindingQty() {
        return selfBindingQty;
    }
}
