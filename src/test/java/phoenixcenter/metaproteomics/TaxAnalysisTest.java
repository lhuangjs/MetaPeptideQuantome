package phoenixcenter.metaproteomics;

import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.dom4j.io.SAXReader;
import org.junit.Before;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class TaxAnalysisTest {
    TaxAnalysis taxAnalysis = new TaxAnalysis();

    @Test
    public void peptide2LCA() throws IOException {
        String peptideWithQuant = "/home/huangjs/Documents/mpq/sample_peptide.tsv";
        String lcaFile = "/home/huangjs/Documents/mpq/LCA_sample_peptide.tsv";
        taxAnalysis.peptide2LCA(peptideWithQuant, true, lcaFile);
    }

    @Test
    public void calPeptideTaxonDistribution() throws IOException {
        String lcaFile = "/home/huangjs/Documents/mpq/LCA_sample_peptide.tsv";
        taxAnalysis.calPeptideTaxonDistribution(lcaFile);
    }

    @Test
    public void dealLake() throws DocumentException, IOException {
        double peptideFDR = 0.01;
        String peptideFile = "/home/huangjs/Documents/mpq/lake/peptides.tsv";
        final SAXReader reader = new SAXReader();
        Map<String, int[]> peptide2PSMCount = new HashMap<>();
        String[] pepxmlFiles = {"/home/huangjs/Documents/mpq/lake/LCM-(01).pep.xml",
                "/home/huangjs/Documents/mpq/lake/GEM-(01).pep.xml"};
        for (int i = 0; i < pepxmlFiles.length; i++) {
            long totalPSM = 0L;
            long validPSM = 0L;
            Document document = reader.read(pepxmlFiles[i]);
            List<Element> queryElementList = document.getRootElement().element("msms_run_summary").elements();
            Iterator<Element> itr = queryElementList.iterator();
            Element queryElement = null;
            while (itr.hasNext()) {
                queryElement = itr.next();
                if (queryElement.getName().equals("spectrum_query")) {
                    break;
                }
            }
            while (queryElement != null) {
                totalPSM += 1L;
                Element searchHit = queryElement.element("search_result").element("search_hit");
                String peptide = searchHit.attributeValue("peptide");
                double qvalue = Double.parseDouble(searchHit.element("analysis_result")
                        .element("percolator_result").attribute("q-Value").getValue());
                if (qvalue < peptideFDR) {
                    validPSM += 1L;
                    int[] psmCount = peptide2PSMCount.get(peptide);
                    if (psmCount == null) {
                        psmCount = new int[pepxmlFiles.length];
                        peptide2PSMCount.put(peptide, psmCount);
                    }
                    psmCount[i] += 1;
                }
                if(itr.hasNext()){
                    queryElement = itr.next();
                }else{
                    queryElement = null;
                }
            }
            System.err.printf("Total PSM: %d, valid PSM: %d\n", totalPSM, validPSM);
        }

        BufferedWriter bw = Files.newBufferedWriter(Paths.get(peptideFile));
        bw.write(String.join("\t",
                "Peptide", "Lake1", "Lake2") + System.lineSeparator());
        for (Map.Entry<String, int[]> e : peptide2PSMCount.entrySet()) {
            bw.write(e.getKey() + "\t" +
                    Arrays.stream(e.getValue()).mapToObj(String::valueOf).collect(Collectors.joining("\t"))
                    + System.lineSeparator());
        }
        bw.close();
    }

    @Test
    public void runLakeLCA() throws Exception {
        String peptideFile = "/home/huangjs/Documents/mpq/lake/peptides.tsv";
        String lcaFile = "/home/huangjs/Documents/mpq/lake/lca.tsv";
        taxAnalysis.peptide2LCA(peptideFile, true, lcaFile);
    }

//    @Test
//    public void name() throws IOException {
//        List<String> peptideList = Files.lines(Paths.get( "/home/huangjs/Documents/mpq/lake/peptides.tsv"))
//                .skip(1)
//                .map(line -> line.split("\t", 2)[0])
//                .collect(Collectors.toList());
//        for (int i = 0; i < peptideList.size(); i++) {
//            String cmd = "unipept pept2lca --equate --all " + peptideList.subList(i, i + 50).stream().collect(Collectors.joining(" "));
//            CommandExecutor.exec(cmd);
//
//        }
//
//    }
}