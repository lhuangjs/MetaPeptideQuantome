package phoenixcenter.metaproteomics;

import org.junit.Test;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Stream;


public class MPQPipelineTest {

    Mascot2XMLEnhancer mascot2XMLEnhancer = new Mascot2XMLEnhancer();

    PeptideProphetEnhancer peptideProphetEnhancer = new PeptideProphetEnhancer();

    private TaxAnalysis taxAnalysis = new TaxAnalysis();

    String decoyPrefix = "DECOY_";

    @Test
    public void dat2pepxml() throws IOException {
        String datFile = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/F008284.dat";
        String decoyPrefix = "DECOY_";
        String targetLibrary = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/biomass.fasta";
        String targetDecoyLibrary = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/biomass-td.fasta";
        String pepxmlFile = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/F008284.pep.xml";
//        mascot2XMLEnhancer.genLibraryWithDecoy(targetLibrary, decoyPrefix, targetDecoyLibrary);
        mascot2XMLEnhancer.convertAutoDecoyMode(datFile, decoyPrefix, targetDecoyLibrary, pepxmlFile);
    }

    @Test
    public void pepxml2pp() throws IOException, FileParsingException {
        String[] datIds = Stream.of("F008283", "F008284", "F008285", "F008286")
                .toArray(String[]::new);

        String targetLibrary = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/biomass.fasta";
        String targetDecoyLibrary = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/biomass-td.fasta";
        String[] ppPepxmlFiles = new String[datIds.length];
        Map<String, String> ppParams = new HashMap<>();
        ppParams.put("-PPM", null);
        for (int i = 0; i < datIds.length; i++) {
            String pepxmlFile = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/" + datIds[i] + ".pep.xml";
            mascot2XMLEnhancer.convertAutoDecoyMode("/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/" + datIds[i] + ".dat",
                    decoyPrefix, targetDecoyLibrary, pepxmlFile);
//            peptideProphetEnhancer.runPeptideProphet( "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/" + datIds[i] + ".pp.pep.xml",
//            ppParams, pepxmlFile);
        }

    }


    @Test
    public void pp2tsv() throws IOException, FileParsingException {
        String peptide2PSMCountFile = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/peptide2psmcount.tsv";
        String[] ppPepxmlFiles = Stream.of("F008283", "F008284", "F008285", "F008286")
                .map(id -> "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/pp-" + id + ".pep.xml")
                .toArray(String[]::new);
        peptideProphetEnhancer.pp2tsv(decoyPrefix, ppPepxmlFiles, 0.01, peptide2PSMCountFile);
    }

    @Test
    public void tsv2lca() throws IOException, FileParsingException {
        String peptide2PSMCountFile = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/peptide2psmcount_end_10000.tsv";
        String lcaFile = "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/peptide2psmcount_end_10000_lca.tsv";
        taxAnalysis.peptide2LCA(peptide2PSMCountFile, true, true, lcaFile);
    }

    @Test
    public void run() throws IOException, FileParsingException {

        // dat2pepxml
        String datFile = "/home/huangjs/Documents/mpq/F007790.dat";
        String decoyPrefix = "DECOY_";
        String targetLibrary = "/home/huangjs/Documents/mpq/biomass-nr.fasta";
        String targetDecoyLibrary = "/home/huangjs/Documents/mpq/biomass-nr.td.fasta";
        String pepxmlFile = "/home/huangjs/Documents/mpq/F007790.pep.xml";
//        mascot2XMLEnhancer.genLibraryWithDecoy(targetLibrary, decoyPrefix, targetDecoyLibrary);
//        mascot2XMLEnhancer.convertAutoDecoyMode(datFile, decoyPrefix, targetDecoyLibrary, pepxmlFile);
        // pepxml2pp
        PeptideProphetEnhancer peptideProphetEnhancer = new PeptideProphetEnhancer();
        Map<String, String> ppParams = new HashMap<>();
        ppParams.put("-PPM", null);
        String ppPepxmlFile = "/home/huangjs/Documents/mpq/filtered-F007790.pep.xml";
//        peptideProphetEnhancer.run(pepxmlFile, ppParams, ppPepxmlFile);
        // pp2peptide

        // peptide2lca

        // lca2quant

    }
}