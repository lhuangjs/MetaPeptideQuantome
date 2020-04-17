package phoenixcenter.metaproteomics;

import org.apache.logging.log4j.core.util.JsonUtils;
import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.dom4j.io.SAXReader;
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
        taxAnalysis.peptide2LCA(peptideWithQuant, true, true, lcaFile);
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
                if (itr.hasNext()) {
                    queryElement = itr.next();
                } else {
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
        taxAnalysis.peptide2LCA(peptideFile, true, true, lcaFile);
    }

    @Test
    public void calBiomassPeptideTaxonDistribution() throws IOException {
        taxAnalysis.calPeptideTaxonDistribution("/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/LCA/final_lca.tsv");
    }

    @Test
    public void lca2Quant() throws IOException {
        taxAnalysis.lca2Quant("/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/LCA/final_lca.tsv",
                "species",
                3,
                false,
                "/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/LCA/taxon2quant_species-3.tsv");
    }


    @Test
    public void compare() throws IOException {
        Map<String, Double> actualAbundance = new HashMap<String, Double>() {
            {
                put("Agrobacterium tumefaciens", 4.259);
                put("Alteromonas macleodii", 4.259);
                put("Bacillus subtilis", 4.259);
                put("Chlamydomonas reinhardtii", 4.259);
                put("Chromobacterium violaceum", 4.259);
                put("Cupriavidus metallidurans", 4.259);
                put("Enterobacteria phage ES18", 0.411);
                put("Escherichia coli", 4.259);
                put("Escherichia virus M13", 0.411);
                put("Escherichia virus MS2", 0.411);
                put("Nitrososphaera viennensis", 4.259);
                put("Paracoccus denitrificans", 4.259);
                put("Pseudomonas fluorescens", 4.259);
                put("Pseudomonas furukawaii", 4.259);
                put("Pseudomonas sp. ATCC 13867", 4.259);
                put("Rhizobium leguminosarum", 8.518);
                put("Roseobacter sp.", 4.259);
                put("Salmonella enterica", 12.777);
                put("Salmonella virus FelixO1", 0.411);
                put("Salmonella virus P22", 0.411);
                put("Staphylococcus aureus", 8.518);
                put("Stenotrophomonas maltophilia", 4.259);
                put("Thermus thermophilus", 4.259);
            }
        };
        double total = actualAbundance.values().stream().mapToDouble(val -> val).sum();
        int sampeSize = 4;
        Map<String, Double> estimateAbundance = Files.lines(Paths.get("/home/huangjs/Documents/mpq/biomass-orgin-Run1_P/LCA/taxon2quant_species-4.tsv"))
                .skip(4)
                .map(line -> {
                    String[] tmp = line.split("\t");
                    double quantVal = 0.0;
                    for (int i = 2; i < 2 + sampeSize; i++) {
                        quantVal += Double.parseDouble(tmp[i]);
                    }
                    return new AbstractMap.SimpleEntry<>(tmp[1], quantVal/ sampeSize);
                })
                .collect(Collectors.toMap(
                        e -> e.getKey(),
                        e -> e.getValue()
                ));
        /** cal MSE **/
        double mse = 0.0;

        System.out.println(String.join("\t",
                "Taxon name",
                "Actual value",
                "Estimate value",
                "Diff"
        ));
        for (Map.Entry<String, Double> e : actualAbundance.entrySet().stream()
                .sorted(Comparator.comparing(Map.Entry::getKey))
                .collect(Collectors.toList())
        ) {
            e.setValue(e.getValue() / total * 100);
            double estVal = 0.0;
            if (estimateAbundance.containsKey(e.getKey())) {
                estVal = estimateAbundance.get(e.getKey());
                System.out.println(e.getKey()
                        + "\t" + e.getValue()
                        + "\t" + estVal
                        + "\t" + (e.getValue() - estVal));
            } else {
                System.err.println(e.getKey()
                        + "\t" + e.getValue()
                        + "\t" + estVal
                        + "\t" + (e.getValue() - estVal));
            }
            mse += Math.pow(e.getValue() - estVal, 2.0);
        }
        mse /= actualAbundance.size();
        System.out.println(mse);
    }

}