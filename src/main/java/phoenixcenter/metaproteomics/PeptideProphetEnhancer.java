package phoenixcenter.metaproteomics;

import lombok.extern.log4j.Log4j2;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.PepXmlParser;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.PeptideprophetResult;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SearchHit;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SpectrumQuery;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Function:
 * 1. run PeptideProphet
 * 2. filter psm by FDR and convert pep.xml to tsv
 */
@Log4j2
public class PeptideProphetEnhancer {

    private final String xinteractLocation = GlobalConfig.getValue("xinteract");

    private final String qvalityLocation = GlobalConfig.getValue("qvality");

    public void runPeptideProphet(String pepxmlFile,
                                  Map<String, String> params,
                                  String ppPepxmlFile) throws IOException, FileParsingException {
        // run PeptideProphetEnhancer
        String command = String.join(" ", xinteractLocation,
                "-N" + ppPepxmlFile,
                params.entrySet().stream()
                        .map(e -> e.getKey() + (e.getValue() == null ? "" : e.getValue()))
                        .collect(Collectors.joining(" ")),
                pepxmlFile);
        CommandExecutor.exec(command);
        // delete unused files
        deleteUnusedFiles(Paths.get(ppPepxmlFile));
    }

    /**
     * Filter by FDR
     *
     * @param ppPepxmlPath  the pepxml generated by PeptideProphet
     * @param decoyPrefix
     * @param fdrThreshold
     * @param fppPepxmlPath the result pepxml file filtered by FDR
     * @throws IOException
     * @throws FileParsingException
     */
    public void filterByFDR(Path ppPepxmlPath, String decoyPrefix,
                            double fdrThreshold, Path fppPepxmlPath) throws IOException, FileParsingException {
        // copy file content
        Files.copy(ppPepxmlPath, fppPepxmlPath, StandardCopyOption.REPLACE_EXISTING);
        log.debug("finish copying {} to {}", fppPepxmlPath);
        double probThreshold = getProbThreshold(ppPepxmlPath, decoyPrefix, fdrThreshold);
        log.info("PeptideProphetEnhancer threshold: {} when FDR = {}", probThreshold, fdrThreshold);
        String command = String.join(" ", xinteractLocation,
                "-N" + fppPepxmlPath,
                "-nI -p" + probThreshold);
        CommandExecutor.exec(command);
        // delete unused files
        deleteUnusedFiles(fppPepxmlPath);
    }

    private double getProbThreshold(Path ppPepxmlPath, String decoyPrefix, double fdrThreshold) throws FileParsingException, IOException {
        // read pepxml file to get probability after run PeptideProphetEnhancer
        Path targetPath = Files.createTempFile(ppPepxmlPath.getParent(), "target", ".txt");
        Path decoyPath = Files.createTempFile(ppPepxmlPath.getParent(), "decoy", ".txt");
        Path qvalityRsPath = Files.createTempFile(ppPepxmlPath.getParent(), "qvality", ".txt");
        log.debug("create temp files for qvality: target file - {} ; decoy file - {} ; result - {}",
                targetPath, decoyPath, qvalityRsPath);
        BufferedWriter targetBW = Files.newBufferedWriter(targetPath);
        BufferedWriter decoyBW = Files.newBufferedWriter(decoyPath);

        List<SpectrumQuery> specQueryList = PepXmlParser.parse(ppPepxmlPath)
                .getMsmsRunSummary()
                .get(0)
                .getSpectrumQuery();
        long targetCount = 0L;
        long decoyCount = 0L;
        for (SpectrumQuery specQuery : specQueryList) {
            SearchHit searchHit = specQuery.getSearchResult().get(0).getSearchHit().get(0);
            // proteins
            Set<String> pidSet = new HashSet<>();
            pidSet.add(searchHit.getProtein());
            if (searchHit.getAlternativeProtein().size() > 0) {
                pidSet.addAll(
                        searchHit.getAlternativeProtein().stream()
                                .map(altProt -> altProt.getProtein())
                                .collect(Collectors.toSet())
                );
            }
            // probability
            double prob = ((PeptideprophetResult) (searchHit.getAnalysisResult().get(0).getAny().get(0)))
                    .getProbability();
            if (pidSet.stream().anyMatch(pid -> !pid.startsWith(decoyPrefix))) {
                // it is target
                targetCount++;
                targetBW.write(prob + System.lineSeparator());
            } else {
                decoyBW.write(prob + System.lineSeparator());
                decoyCount++;
            }
        }
        targetBW.close();
        decoyBW.close();
        log.debug("{}: target = {}, decoy = {}", ppPepxmlPath, targetCount, decoyCount);

        /** run qvality **/
        String command = String.join(" ", qvalityLocation,
                targetPath.toString(), decoyPath.toString(),
                "-d -Y",
                "-o", qvalityRsPath.toString()
        );
        CommandExecutor.exec(command);

        /** get min probability when FDR < fdrThreshold **/
        BufferedReader br = Files.newBufferedReader(qvalityRsPath);
        String line;
        double probThreshold = -1.0;
        br.readLine();
        while ((line = br.readLine()) != null) {
            String[] tmp = line.split("\t");
            if (Double.parseDouble(tmp[2]) >= fdrThreshold) {
                break;
            } else {
                probThreshold = Double.parseDouble(tmp[0]);
            }
        }
        br.close();
        // delete temporary files
        Files.delete(targetPath);
        Files.delete(decoyPath);
        Files.delete(qvalityRsPath);
        log.debug("delete temp files for qvality: target file - {} ; decoy file - {} ; result - {}",
                targetPath, decoyPath, qvalityRsPath);
        return probThreshold;
    }

    private void deleteUnusedFiles(Path pepxmlPath) throws IOException {
        String pepxmlName = pepxmlPath.getFileName().toString();
        pepxmlName = pepxmlName.substring(0, pepxmlName.lastIndexOf('.'));
        Path ppIndexPath = pepxmlPath.getParent().resolve(pepxmlName + ".xml.index");
        Path ppModelsPath = pepxmlPath.getParent().resolve(pepxmlName + "-MODELS.html");
        log.debug("{} will be deleted", ppIndexPath);
        log.debug("{} will be deleted", ppModelsPath);
        Files.delete(ppIndexPath);
        Files.delete(ppModelsPath);
    }

    public void pp2tsv(String decoyPrefix, String[] ppPepxmlFiles, double fdrThreshold, String peptide2PSMCountFile) throws IOException, FileParsingException {
        double[] probThresholds = new double[ppPepxmlFiles.length];
        for (int i = 0; i < ppPepxmlFiles.length; i++) {
            probThresholds[i] = getProbThreshold(Paths.get(ppPepxmlFiles[i]), decoyPrefix, fdrThreshold);
        }
        pp2tsv(decoyPrefix, ppPepxmlFiles, probThresholds, peptide2PSMCountFile);
    }

    public void pp2tsv(String decoyPrefix, String[] ppPepxmlFiles, String peptide2PSMCountFile) throws IOException, FileParsingException {
        double[] probThresholds = IntStream.range(0, ppPepxmlFiles.length)
                .mapToDouble(i -> 0.0)
                .toArray();
        pp2tsv(decoyPrefix, ppPepxmlFiles, probThresholds, peptide2PSMCountFile);
    }

    private void pp2tsv(String decoyPrefix, String[] ppPepxmlFiles,
                        double[] probThresholds, String peptide2PSMCountFile) throws FileParsingException, IOException {
        Map<String, int[]> peptide2PSMCounts = new HashMap<>();
        for (int i = 0; i < ppPepxmlFiles.length; i++) {
            Map<String, Integer> peptide2PSMCount = statPSMCount(decoyPrefix, ppPepxmlFiles[i], probThresholds[i]);
            final int repIdx = i;
            peptide2PSMCount.entrySet().stream()
                    .forEach(e -> {
                        int[] psmCounts = peptide2PSMCounts.get(e.getKey());
                        if (psmCounts == null) {
                            psmCounts = new int[ppPepxmlFiles.length];
                            peptide2PSMCounts.put(e.getKey(), psmCounts);
                        }
                        psmCounts[repIdx] = e.getValue();
                    });
        }
        BufferedWriter bw = Files.newBufferedWriter(Paths.get(peptide2PSMCountFile));
        String sampleNames = Arrays.stream(ppPepxmlFiles)
                .map(file -> Paths.get(file))
                .map(path -> {
                    String sampleName = path.getName(path.getNameCount() - 1).toString();
                    int suffixIdx;
                    if (sampleName.endsWith(".pep.xml")) {
                        return sampleName.substring(0, sampleName.length() - 8);
                    } else if ((suffixIdx = sampleName.lastIndexOf(".")) != -1) {
                        return sampleName.substring(0, suffixIdx);
                    }
                    return sampleName;
                })
                .collect(Collectors.joining("\t"));
        bw.write("Peptide\t" + sampleNames + System.lineSeparator());
        for (Map.Entry<String, int[]> e : peptide2PSMCounts.entrySet()) {
            bw.write(e.getKey() + "\t"
                    + Arrays.stream(e.getValue()).mapToObj(count -> String.valueOf(count)).collect(Collectors.joining("\t"))
                    + System.lineSeparator());
        }
        bw.close();
    }

    private Map<String, Integer> statPSMCount(String decoyPrefix, String ppPepxmlFile,
                                              double probThreshold) throws FileParsingException {
        List<SpectrumQuery> specQueryList = PepXmlParser.parse(Paths.get(ppPepxmlFile))
                .getMsmsRunSummary()
                .get(0)
                .getSpectrumQuery();
        Map<String, Integer> peptide2PSMCount = new HashMap<>();
        for (SpectrumQuery specQuery : specQueryList) {
            SearchHit searchHit = specQuery.getSearchResult().get(0).getSearchHit().get(0);
            // peptide
            String peptide = searchHit.getPeptide();
            // proteins
            Set<String> pidSet = new HashSet<>();
            pidSet.add(searchHit.getProtein());
            if (searchHit.getAlternativeProtein().size() > 0) {
                pidSet.addAll(
                        searchHit.getAlternativeProtein().stream()
                                .map(altProt -> altProt.getProtein())
                                .collect(Collectors.toSet())
                );
            }
            // target
            if (pidSet.stream().anyMatch(pid -> !pid.startsWith(decoyPrefix))) {
                double prob = ((PeptideprophetResult) (searchHit.getAnalysisResult().get(0).getAny().get(0)))
                        .getProbability();
                if (prob >= probThreshold) {
                    Integer psmCount = peptide2PSMCount.get(peptide);
                    if (psmCount == null) {
                        peptide2PSMCount.put(peptide, 1);
                    } else {
                        peptide2PSMCount.put(peptide, psmCount + 1);
                    }
                }
            }
        }
        return peptide2PSMCount;
    }
}
