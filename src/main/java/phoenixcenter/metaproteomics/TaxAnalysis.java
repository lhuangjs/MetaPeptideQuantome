package phoenixcenter.metaproteomics;

import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.extern.log4j.Log4j2;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.Executors;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

@Log4j2
public class TaxAnalysis {

    private CommandExecutor commandExecutor = new CommandExecutor(Executors.newFixedThreadPool(6));

    private final String python3 = GlobalConfig.getValue("python3");
    private final String peptideRankDistScript = GlobalConfig.getValue("peptide-rank-dist.script");
    private final String peptideTaxonDistScript = GlobalConfig.getValue("peptide-taxon-dist.script");

    private final ObjectMapper objectMapper = new ObjectMapper();

    private final int batchSize = GlobalConfig.getIntValue("unipept.batch.size");

    private final String unipeptHeader = "taxon_id	taxon_name	taxon_rank	superkingdom_id	superkingdom_name	kingdom_id	kingdom_name	subkingdom_id	subkingdom_name	superphylum_id	superphylum_name	phylum_id	phylum_name	subphylum_id	subphylum_name	superclass_id	superclass_name	class_id	class_name	subclass_id	subclass_name	infraclass_id	infraclass_name	superorder_id	superorder_name	order_id	order_name	suborder_id	suborder_name	infraorder_id	infraorder_name	parvorder_id	parvorder_name	superfamily_id	superfamily_name	family_id	family_name	subfamily_id	subfamily_name	tribe_id	tribe_name	subtribe_id	subtribe_name	genus_id	genus_name	subgenus_id	subgenus_name	species_group_id	species_group_name	species_subgroup_id	species_subgroup_name	species_id	species_name	subspecies_id	subspecies_name	varietas_id	varietas_name	forma_id	forma_name";

    private final String[] ranks = {"forma", "varietas", "subspecies", "species", "species_subgroup", "species_group", "subgenus", "genus", "subtribe",
            "tribe", "subfamily", "family", "superfamily", "parvorder", "infraorder", "suborder", "order", "superorder",
            "infraclass", "subclass", "class", "superclass", "subphylum", "phylum", "superphylum", "subkingdom", "kingdom",
            "superkingdom"
    };

    private final Map<String, Integer> rank2SortIdx = IntStream.range(0, ranks.length)
            .mapToObj(Integer::valueOf)
            .collect(Collectors.toMap(
                    i -> ranks[i],
                    i -> i
            ));

    public void peptide2LCA(String peptideFile,
                            boolean equalIL,
                            String lcaFile) throws IOException {
        BufferedReader br = Files.newBufferedReader(Paths.get(peptideFile));
        BufferedWriter bw = Files.newBufferedWriter(Paths.get(lcaFile));
        // peptide file contains quant data?
        String[] tmp = Arrays.stream(br.readLine().split("\t", 2)).map(e -> e.trim()).toArray(String[]::new);
        boolean containQuantData = tmp.length > 1;
        // write header
        if (containQuantData) {
            bw.write(String.join("\t", tmp[0], tmp[1], unipeptHeader) + System.lineSeparator());
        } else {
            bw.write(String.join("\t", tmp[0], unipeptHeader) + System.lineSeparator());
        }
        // read peptide file
        String line;
        long count = 0L;
        int sleepBatch = batchSize * 5;
        List<String[]> peptideInfoList = new ArrayList<>(batchSize);
        while ((line = br.readLine()) != null) {
            count += 1L;
            peptideInfoList.add(line.split("\t", 2));
            if (peptideInfoList.size() == batchSize) {
                processShard(peptideInfoList, containQuantData, equalIL, bw);
                bw.flush();
                peptideInfoList.clear();
                log.info("{} peptides have been processed", count);
            }
            if (count % sleepBatch == 0) {
                try {
                    Thread.sleep(60_000);
                } catch (InterruptedException e) {
                    throw new RuntimeException(e);
                }
            }
        }
        if (peptideInfoList.size() > 0) {
            processShard(peptideInfoList, containQuantData, equalIL, bw);
        }
        bw.close();
        br.close();
        log.info("total {} peptides have been processed", count);
        Path parentDirPath = Paths.get(lcaFile).getParent();
        String rank2PeptideCountJsonFile = Files.createTempFile("rank2PeptideCount", ".json").toString();
        String peptideCount2TaxonCountJsonFile = Files.createTempFile("peptideCount2TaxonCount", ".json").toString();
        calPeptideTaxonDistribution(lcaFile);
        log.debug("rank2PeptideCountJsonFile: {}, peptideCount2TaxonCountJsonFile: {} ",
                rank2PeptideCountJsonFile, peptideCount2TaxonCountJsonFile);
    }

    private void processShard(List<String[]> peptideInfoList,
                              boolean containQuantData,
                              boolean equalIL,
                              BufferedWriter bw) throws IOException {
        Map<String, List<String>> peptide2Info = peptideInfoList.stream()
                .collect(Collectors.groupingBy(peptInfo -> peptInfo[0],
                        Collectors.mapping((String[] peptInfo) -> String.join("\t", peptInfo),
                                Collectors.toList())
                ));
        Map<String, String> peptide2LCA = new HashMap<>();
        peptide2LCA(peptide2Info.keySet(), equalIL, lcaInfo -> peptide2LCA.put(lcaInfo[0], lcaInfo[1]));
        for (Map.Entry<String, List<String>> e : peptide2Info.entrySet()) {
            String peptide = e.getKey();
            for (String peptideInfo : e.getValue()) {
                String lcaInfo = peptide2LCA.get(peptide);
                if (lcaInfo != null) {
                    bw.write(peptideInfo + "\t" + lcaInfo + System.lineSeparator());
                } else {
                    bw.write(peptideInfo + System.lineSeparator());
                }
            }
        }
//        List<String> peptideList = peptideInfoList.stream()
//                .map(entry -> entry[0])
//                .collect(Collectors.toList());
//        List<String[]> lcaInfoList = new ArrayList<>(peptideList.size());
//        peptide2LCA(peptideList, equalIL, lcaInfo -> lcaInfoList.add(lcaInfo));
//        Iterator<String[]> peptideInfoItr = peptideInfoList.iterator();
//        Iterator<String[]> lcaInfoItr = lcaInfoList.iterator();
//        while (lcaInfoItr.hasNext()) {
//            String[] lcaInfo = lcaInfoItr.next();
//            while (peptideInfoItr.hasNext()) {
//                StringJoiner joiner = new StringJoiner("\t");
//                String[] peptideInfo = peptideInfoItr.next();
//                joiner.add(peptideInfo[0]);
//                if (containQuantData) {
//                    joiner.add(peptideInfo[1]);
//                }
//                if (lcaInfo[0].equals(peptideInfo[0])) {
//                    joiner.add(lcaInfo[1].replace(",", "\t"));
//                    bw.write(joiner.toString() + System.lineSeparator());
//                    break;
//                }
//            }
//        }
//        while (peptideInfoItr.hasNext()) {
//            if (containQuantData) {
//                bw.write(String.join("\t", peptideInfoItr.next()) + System.lineSeparator());
//            } else {
//                bw.write(peptideInfoItr.next() + System.lineSeparator());
//            }
//        }
    }

    public void peptide2LCA(Collection<String> peptides,
                            boolean equalIL,
                            Consumer<String[]> consumer) throws IOException {
//        Path tempFilePath = Files.createTempFile("unipept", ".tsv");
        String command = String.join(" ",
                "unipept pept2lca",
//                "--output", tempFilePath.toString(),
                equalIL ? "--equate" : "",
                "--all",
                peptides.stream().collect(Collectors.joining(" "))
        );
        commandExecutor.exec(command, (String line) -> {
            if (!line.startsWith("peptide")) {
                consumer.accept(line.split(",", 2));
            }
        });
//        Files.lines(tempFilePath)
//                .skip(1)
//                .forEach(line -> consumer.accept(line.split(",", 2)));
//        Files.delete(tempFilePath);
    }

    public void calPeptideTaxonDistribution(String lcaFile) throws IOException {
        Path lcaFilePath = Paths.get(lcaFile);
        BufferedReader br = Files.newBufferedReader(lcaFilePath);
        String[] tmp = br.readLine().split("\t");
        Map<String, Integer> rank2Index = Arrays.stream(ranks)
                .collect(Collectors.toMap(
                        r -> r,
                        r -> -1
                ));
        int tidIdx = -1;
        int trankIdx = -1;
        for (int i = 0; i < tmp.length; i++) {
            if (tmp[i].equals("taxon_id")) {
                tidIdx = i;
            }
            if (tmp[i].endsWith("_name")) {
                String rank = tmp[i].substring(0, tmp[i].length() - 5);
                if (rank2Index.containsKey(rank)) {
                    rank2Index.put(rank, i);
                }
            } else if (tmp[i].equals("taxon_id")) {
                tidIdx = i;
            } else if (tmp[i].equals("taxon_rank")) {
                trankIdx = i;
            }
        }
        // calculate distribution
        String line;
        List<String[]> lineageList = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            tmp = line.split("\t", rank2Index.get(ranks[0]) + 1);
            if (tidIdx < tmp.length) {
                String[] lineage = new String[ranks.length];
                lineageList.add(lineage);
                for (int i = 0; i < ranks.length; i++) {
                    int rankIdx = rank2Index.get(ranks[i]);
                    String tname = tmp[rankIdx];
                    if (tname.length() > 0) {
                        lineage[i] = tname;
                    } else {
                        // 如果当前rank是低于LCA的rank，那么说明说明当前rank不存在taxon，否则，当前rank暂时没有分配名字。
                        // 等级越高数字越大
                        if (rank2SortIdx.get(ranks[i]) < rank2SortIdx.get(tmp[trankIdx])) {
                            lineage[i] = null;
                        } else {
                            lineage[i] = ranks[i] + "_" + tmp[tidIdx];
                        }
                    }
                }
            }
        }
        br.close();
        // cal distribution
        Map<String, PeptideCount> rank2PeptideCount = new HashMap<>();
        Map<String, Map<Integer, Integer>> peptideCount2TaxonCountInRank = new HashMap<>();
        int peptideCountForSubrank = 0;
        for (int i = 0; i < ranks.length; i++) {
            final int rankIdx = i;
            Map<String, Integer> taxon2PeptideCount = lineageList.parallelStream()
                    .map(lineage -> lineage[rankIdx])
                    .filter(tname -> tname != null)
                    .collect(Collectors.groupingBy(tname -> tname, Collectors.reducing(0, e -> 1, Integer::sum)));
            int peptideCount = taxon2PeptideCount.values().stream().mapToInt(Integer::intValue).sum();
            rank2PeptideCount.put(ranks[i], new PeptideCount(peptideCount - peptideCountForSubrank, peptideCountForSubrank));
            System.out.println(String.join("\t", ranks[i], (peptideCount - peptideCountForSubrank) + "", peptideCountForSubrank + ""));
            peptideCountForSubrank = peptideCount;
            peptideCount2TaxonCountInRank.put(ranks[i],
                    taxon2PeptideCount.entrySet().stream()
                            .collect(Collectors.groupingBy((Map.Entry<String, Integer> e) -> e.getValue(),
                                    Collectors.reducing(0, e -> 1, Integer::sum)))
            );
        }
        // write
        Path parentPath = lcaFilePath.getParent();
        String lcaFileName = lcaFilePath.getName(lcaFilePath.getNameCount() - 1).toString();
        int suffixIdx = lcaFileName.lastIndexOf(".");
        lcaFileName = suffixIdx == -1 ? lcaFileName : lcaFileName.substring(0, suffixIdx);
        Path peptideRankDistChartPath = parentPath.resolve(lcaFileName + "_peptide_distribution.png");
        Path peptideTaxonDistChartPath = parentPath.resolve(lcaFileName + "_peptide_taxon_distribution.png");
        Path rank2PeptideCountJsonPath = Files.createTempFile("rank2PeptideCount", ".json");
        Path peptideCount2TaxonCountJsonPath = Files.createTempFile("peptideCount2TaxonCount", ".json");
        objectMapper.writeValue(rank2PeptideCountJsonPath.toFile(), rank2PeptideCount);
        objectMapper.writeValue(peptideCount2TaxonCountJsonPath.toFile(), peptideCount2TaxonCountInRank);
        // plot
        String peptideRankDistChartCmd = String.join(" ",
                python3,
                peptideRankDistScript,
                rank2PeptideCountJsonPath.toString(),
                peptideRankDistChartPath.toString()
        );
        commandExecutor.exec(peptideRankDistChartCmd);
        String peptideTaxonDistChartCmd = String.join(" ",
                python3,
                peptideTaxonDistScript,
                peptideCount2TaxonCountJsonPath.toString(),
                peptideTaxonDistChartPath.toString()
        );
        commandExecutor.exec(peptideTaxonDistChartCmd);
        Files.delete(peptideRankDistChartPath);
        Files.delete(peptideTaxonDistChartPath);
    }

    public void lca2Quant(String lcaFile, String peptide2QuantFile, boolean equateIL) {


    }

    @Data
    @NoArgsConstructor
    @AllArgsConstructor
    class PeptideCount {
        private int peptideCountForRank;

        private int peptideCountForSubrank;
    }
}