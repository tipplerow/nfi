
package nfi.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import jam.app.JamApp;
import jam.app.JamLogger;
import jam.app.JamProperties;
import jam.io.IOUtil;
import jam.util.ListUtil;
import jam.util.StreamUtil;

import jene.hla.Allele;
import jene.hla.Genotype;
import jene.neo.PeptidePairRecord;
import jene.neo.PeptidePairTable;
import jene.tcga.TumorBarcode;
import jene.tcga.TumorGenotypeTable;

/**
 * Computes allele footprint index scores for a patient cohort.
 */
public final class AlleleFootprintDriver extends JamApp {
    private final String footprintFile;
    private final String peptidePairFile;
    private final String tumorPatientFile;
    private final String patientGenotypeFile;

    private final AlleleFootprintType footprintType;
    private final AlleleFootprintIndex footprintIndex;

    private PeptidePairTable peptidePairTable;
    private TumorGenotypeTable tumorGenotypeTable;
    private List<TumorBarcode> tumorBarcodes;
    private List<AlleleFootprintRecord> footprintRecords;

    private AlleleFootprintDriver(String... propFiles) {
        super(propFiles);

        this.footprintFile = resolveFootprintFile();
        this.peptidePairFile = resolvePeptidePairFile();
        this.tumorPatientFile = resolveTumorPatientFile();
        this.patientGenotypeFile = resolvePatientGenotypeFile();

        this.footprintType = resolveFootprintType();
        this.footprintIndex = footprintType.getAlleleFootprintIndex();
    }

    private static String resolveFootprintFile() {
        return JamProperties.getRequired(FOOTPRINT_FILE_PROPERTY);
    }

    private static AlleleFootprintType resolveFootprintType() {
        return JamProperties.getRequiredEnum(FOOTPRINT_TYPE_PROPERTY, AlleleFootprintType.class);
    }

    private static String resolvePeptidePairFile() {
        return JamProperties.getRequired(PEPTIDE_PAIR_FILE_PROPERTY);
    }

    private static String resolvePatientGenotypeFile() {
        return JamProperties.getRequired(PATIENT_GENOTYPE_FILE_PROPERTY);
    }

    private static String resolveTumorPatientFile() {
        return JamProperties.getRequired(TUMOR_PATIENT_FILE_PROPERTY);
    }

    /**
     * Name of the system property that specifies the full path name
     * of the output footprint file.
     */
    public static final String FOOTPRINT_FILE_PROPERTY = "AlleleFootprintDriver.footprintFile";

    /**
     * Name of the system property that specifies the allele footprint
     * calculation type to employ.
     */
    public static final String FOOTPRINT_TYPE_PROPERTY = "AlleleFootprintDriver.footprintType";

    /**
     * Name of the system property that specifies the full path name
     * of the input file mapping patients to their HLA genotypes.
     */
    public static final String PATIENT_GENOTYPE_FILE_PROPERTY = "AlleleFootprintDriver.patientGenotypeFile";

    /**
     * Name of the system property that specifies the full path name
     * of the input file containing neo/self-peptide pairs.
     */
    public static final String PEPTIDE_PAIR_FILE_PROPERTY = "AlleleFootprintDriver.peptidePairFile";

    /**
     * Name of the system property that specifies the full path name
     * of the input file mapping tumor barcodes to patient idenifiers.
     */
    public static final String TUMOR_PATIENT_FILE_PROPERTY = "AlleleFootprintDriver.tumorPatientFile";

    /**
     * Computes allele footprint index scores for a patient cohort.
     *
     * @param propFiles files containing the system properties that
     * define the runtime environment.
     *
     * @throws RuntimeException if any errors occur.
     */
    public static void run(String... propFiles) {
        AlleleFootprintDriver driver = new AlleleFootprintDriver(propFiles);
        driver.run();
    }

    private void run() {
        loadTables();
        sortBarcodes();
        processBarcodes();
        writeFootprints();

        JamLogger.info("DONE!");
    }

    private void loadTables() {
        peptidePairTable = PeptidePairTable.load(peptidePairFile);
        tumorGenotypeTable = TumorGenotypeTable.load(tumorPatientFile, patientGenotypeFile);
    }

    private void sortBarcodes() {
        tumorBarcodes = new ArrayList<TumorBarcode>(peptidePairTable.viewBarcodes());
        Collections.sort(tumorBarcodes);
    }

    private void processBarcodes() {
        List<List<AlleleFootprintRecord>> barcodeLists =
            StreamUtil.applyParallel(tumorBarcodes, barcode -> processBarcode(barcode));

        JamLogger.info("Concatenating footprint records...");
        footprintRecords = ListUtil.cat(barcodeLists);

        JamLogger.info("Sorting footprint records...");
        footprintRecords.sort(AlleleFootprintRecord.COMPARATOR);
    }

    private List<AlleleFootprintRecord> processBarcode(TumorBarcode barcode) {
        JamLogger.info("Processing [%s]...", barcode);

        try {
            Genotype patientGenotype = tumorGenotypeTable.require(barcode);
            Set<Allele> patientAlleles = patientGenotype.viewUniqueAlleles();
            List<PeptidePairRecord> peptidePairRecords = peptidePairTable.lookup(barcode);

            return footprintIndex.compute(patientAlleles, peptidePairRecords);
        }
        catch (RuntimeException ex) {
            JamLogger.warn(ex);
            return List.of();
        }
    }

    private void writeFootprints() {
        JamLogger.info("Writing [%s]...", footprintFile);
        IOUtil.writeLines(footprintFile, false, AlleleFootprintRecord.header());
        IOUtil.writeObjects(footprintFile, true, footprintRecords, record -> record.format());
    }

    private static void usage() {
        System.err.println("Usage: jam.neo.AlleleFootprintDriver PROP_FILE1 [PROP_FILE2 ...]");
        System.exit(1);
    }

    public static void main(String[] args) {
        if (args.length < 1)
            usage();

        run(args);
    }
}
