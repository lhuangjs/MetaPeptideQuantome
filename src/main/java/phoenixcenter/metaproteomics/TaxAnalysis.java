package phoenixcenter.metaproteomics;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.extern.log4j.Log4j2;
import org.apache.http.HttpHeaders;
import org.apache.http.NoHttpResponseException;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.StringEntity;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;
import phoenixcenter.metaproteomics.entity.QuantPeptide;
import phoenixcenter.metaproteomics.entity.UnipeptLCA;
import phoenixcenter.metaproteomics.entity.UnipeptTaxon;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

@Log4j2
public class TaxAnalysis {

    private final String python3 = GlobalConfig.getValue("python3");

    private final String peptideRankDistScript = GlobalConfig.getValue("peptide-rank-dist.script");

    private final String peptideTaxonDistScript = GlobalConfig.getValue("peptide-taxon-dist.script");

    private final String taxonQuantScript = GlobalConfig.getValue("taxon-quant.script");

    private final int pept2dataSleepSecond = GlobalConfig.getIntValue("unipept.pept2data.sleep.second");

    private final int pept2dataRetrySecond = GlobalConfig.getIntValue("unipept.pept2data.retry.second");

    private final CloseableHttpClient httpClient = HttpClients.createDefault();

    private final Map<Integer, UnipeptTaxon> tid2Taxon = new ConcurrentHashMap<>(3000);

    private final ObjectMapper objectMapper = new ObjectMapper();

    private final double log2Val = Math.log(2);

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
                            boolean missedCleavage,
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
        List<String[]> peptideInfoList = new ArrayList<>(batchSize);
        while ((line = br.readLine()) != null) {
            count += 1L;
            peptideInfoList.add(line.split("\t", 2));
            if (peptideInfoList.size() == batchSize) {
                processShard(peptideInfoList, equalIL, missedCleavage, bw);
                bw.flush();
                peptideInfoList.clear();
                try {
                    Thread.sleep(pept2dataSleepSecond * 1_000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                log.info("{} peptides have been processed", count);
            }
        }
        if (peptideInfoList.size() > 0) {
            processShard(peptideInfoList, equalIL, missedCleavage, bw);
        }
        bw.close();
        br.close();
        log.info("total {} peptides have been processed", count);
        calPeptideTaxonDistribution(lcaFile);
    }

    private void processShard(List<String[]> peptideInfoList,
                              boolean equalIL,
                              boolean missedCleavage,
                              BufferedWriter bw) throws IOException {
        Map<String, List<String>> peptide2Info = peptideInfoList.stream()
                .collect(Collectors.groupingBy(peptInfo -> peptInfo[0],
                        Collectors.mapping((String[] peptInfo) -> String.join("\t", peptInfo),
                                Collectors.toList())
                ));
        Map<String, String> peptide2LCA = new HashMap<>();
        peptide2LCA(peptide2Info.keySet(), equalIL, missedCleavage, unipeptLCA -> {
            List<String> lineageNames = unipeptLCA.getLineageNames();
            List<Integer> lineageIds = unipeptLCA.getLineageIds();
            peptide2LCA.put(unipeptLCA.getSequence(),
                    String.join("\t", unipeptLCA.getLcaId().toString(),
                            unipeptLCA.getLcaName(), unipeptLCA.getLcaRank())
                            + "\t"
                            + IntStream.range(0, ranks.length)
                            .mapToObj(i -> lineageIds.get(i) == null ? "\t" : lineageIds.get(i) + "\t" + lineageNames.get(i))
                            .collect(Collectors.joining("\t"))
            );
        });
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
    }

    public void peptide2LCA(Collection<String> peptides,
                            boolean equalIL,
                            boolean missedCleavage,
                            Consumer<UnipeptLCA> consumer) throws IOException {

        /** LCA search **/
        HttpPost lcaRequest = new HttpPost("https://unipept.ugent.be/mpa/pept2data");
        lcaRequest.addHeader(HttpHeaders.ACCEPT, "application/json");
        lcaRequest.addHeader(HttpHeaders.CONTENT_TYPE, "application/json");
        Map<String, Object> lcaParams = new HashMap<>();
        lcaParams.put("peptides", peptides);
        lcaParams.put("equate_il", equalIL);
        lcaParams.put("missed", missedCleavage);
        lcaRequest.setEntity(new StringEntity(objectMapper.writeValueAsString(lcaParams), StandardCharsets.UTF_8));
        // LCA response
        CloseableHttpResponse lcaResponse = null;
        int retryCount = 1;
        do {
            try {
                lcaResponse = httpClient.execute(lcaRequest);
            } catch (NoHttpResponseException ntre) {
                try {
                    Thread.sleep(pept2dataRetrySecond * retryCount);
                    log.warn("network error when execute LCA api, retry {} time", retryCount++);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                continue;
            }
            if (lcaResponse.getStatusLine().getStatusCode() != 200) {
                try {
                    Thread.sleep(pept2dataRetrySecond * retryCount);
                    log.warn("network error when execute LCA api, retry {} time", retryCount++);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            } else {
                break;
            }
        } while (true);
        List<UnipeptLCA> unipeptLCAList = new ArrayList<>();
        Set<Integer> unsearchTidSet = new HashSet<>();
        String lcaInfos = EntityUtils.toString(lcaResponse.getEntity());
        Iterator<JsonNode> peptItr = objectMapper.readTree(lcaInfos).get("peptides").elements();
        while (peptItr.hasNext()) {
            JsonNode peptNode = peptItr.next();
            UnipeptLCA unipeptLCA = new UnipeptLCA();
            unipeptLCA.setSequence(peptNode.get("sequence").asText());
            unipeptLCA.setLcaId(peptNode.get("lca").asInt());
            List<Integer> lineageIds = new ArrayList<>();
            Iterator<JsonNode> idsItr = peptNode.get("lineage").elements();
            while (idsItr.hasNext()) {
                Number tid = idsItr.next().numberValue();
                lineageIds.add(tid == null ? null : tid.intValue());
            }
            unipeptLCA.setLineageIds(lineageIds);
            unipeptLCAList.add(unipeptLCA);
            // add unsearch tid
            if (!tid2Taxon.containsKey(unipeptLCA.getLcaId())) {
                unsearchTidSet.add(unipeptLCA.getLcaId());
            }
            unipeptLCA.getLineageIds().stream()
                    .filter(tid -> tid != null && !tid2Taxon.containsKey(tid))
                    .forEach(tid -> unsearchTidSet.add(tid));
        }
        /** taxon search **/
        if (unsearchTidSet.size() > 0) {
            HttpPost taxaRequest = new HttpPost("https://unipept.ugent.be/private_api/taxa");
            taxaRequest.addHeader(HttpHeaders.ACCEPT, "application/json");
            taxaRequest.addHeader(HttpHeaders.CONTENT_TYPE, "application/json");
            Map<String, Object> taxaParams = new HashMap<>();
            taxaParams.put("taxids", unsearchTidSet);
            taxaRequest.setEntity(new StringEntity(objectMapper.writeValueAsString(taxaParams), StandardCharsets.UTF_8));
            CloseableHttpResponse taxaResponse = null;
            retryCount = 1;
            do {
                try {
                    taxaResponse = httpClient.execute(taxaRequest);
                } catch (NoHttpResponseException nre) {
                    try {
                        Thread.sleep(pept2dataRetrySecond * retryCount);
                        log.warn("network error when execute taxa api, retry {} time", retryCount++);
                    } catch (InterruptedException ie) {
                        throw new RuntimeException(ie);
                    }
                    continue;
                }
                if (taxaResponse.getStatusLine().getStatusCode() != 200) {
                    try {
                        Thread.sleep(pept2dataRetrySecond * retryCount);
                        log.warn("network error when execute taxa api, retry {} time", retryCount++);
                    } catch (InterruptedException ie) {
                        throw new RuntimeException(ie);
                    }
                } else {
                    break;
                }
            } while (true);
            String taxonInfo = EntityUtils.toString(taxaResponse.getEntity());
            List<UnipeptTaxon> taxonList = objectMapper.readValue(taxonInfo, new TypeReference<List<UnipeptTaxon>>() {
            });
            taxonList.stream()
                    .forEach(taxon -> tid2Taxon.put(taxon.getId(), taxon));
        }
        unipeptLCAList.stream()
                .forEach(unipeptLCA -> {
                    // set tname
                    UnipeptTaxon taxon = tid2Taxon.get(unipeptLCA.getLcaId());
                    unipeptLCA.setLcaName(taxon.getName());
                    unipeptLCA.setLcaRank(taxon.getRank());
                    List<Integer> lineageIds = unipeptLCA.getLineageIds();
                    List<String> lineageNames = new ArrayList<>(lineageIds.size());
                    unipeptLCA.setLineageNames(lineageNames);
                    for (int i = 0; i < lineageIds.size(); i++) {
                        lineageNames.add(null);
                        if (lineageIds.get(i) != null) {
                            lineageNames.set(i, tid2Taxon.get(lineageIds.get(i)).getName());
                        }
                    }
                    // consume
                    consumer.accept(unipeptLCA);
                });
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
        int lcaIdIdx = -1;
        int lcaRankIdx = -1;
        for (int i = 0; i < tmp.length; i++) {
            if (tmp[i].equals("taxon_id")) {
                lcaIdIdx = i;
            }
            if (tmp[i].endsWith("_name")) {
                String rank = tmp[i].substring(0, tmp[i].length() - 5);
                if (rank2Index.containsKey(rank)) {
                    rank2Index.put(rank, i);
                }
            } else if (tmp[i].equals("taxon_id")) {
                lcaIdIdx = i;
            } else if (tmp[i].equals("taxon_rank")) {
                lcaRankIdx = i;
            }
        }
        // calculate distribution
        String line;
        List<String[]> lineageList = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            tmp = line.split("\t", rank2Index.get(ranks[0]) + 1);
            if (lcaIdIdx < tmp.length) {
                // skip when lca is root
                if (tmp[lcaIdIdx].equals("1")) {
                    continue;
                }
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
                        if (rank2SortIdx.get(ranks[i]) < rank2SortIdx.get(tmp[lcaRankIdx].replace(" ", "_"))) {
                            lineage[i] = null;
                        } else {
                            lineage[i] = ranks[i] + "_" + tmp[lcaIdIdx];
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
        CommandExecutor.exec(peptideRankDistChartCmd);
        String peptideTaxonDistChartCmd = String.join(" ",
                python3,
                peptideTaxonDistScript,
                peptideCount2TaxonCountJsonPath.toString(),
                peptideTaxonDistChartPath.toString()
        );
        CommandExecutor.exec(peptideTaxonDistChartCmd);
//        Files.delete(rank2PeptideCountJsonPath);
//        Files.delete(peptideCount2TaxonCountJsonPath);
    }

    public void lca2Quant(String lcaFile,
                          String rank,
                          int minPeptideForTaxon,
                          boolean log2,
                          String taxonQuantFile) throws IOException {
        @Data
        class TaxonQuantChartData {
            String rank;
            boolean log2;
            List<String> samples;
            Map<String, List<Double>> taxon2Quants;
        }

        class TaxonQuant {
            UnipeptTaxon unipeptTaxon;
            List<QuantPeptide> quantPeptideList;
            List<Double> quantValList;

            public TaxonQuant(UnipeptTaxon unipeptTaxon, List<QuantPeptide> quantPeptideList) {
                this.unipeptTaxon = unipeptTaxon;
                this.quantPeptideList = quantPeptideList;
            }
        }
        /** filter **/
        List<TaxonQuant> taxonQuantList = lca2Quant(lcaFile, rank)
                .entrySet().stream()
                .filter(e -> e.getValue().size() >= minPeptideForTaxon)
                .map(e -> new TaxonQuant(e.getKey(), e.getValue()))
                .collect(Collectors.toList());
        /** write result file **/
        String headers = Files.lines(Paths.get(lcaFile)).limit(1).findFirst().get();
        int sampleStartIdx = headers.indexOf("\t") + 1;
        int sampleEndIdx = headers.indexOf("taxon_id") - 1;
        String samplesStr = headers.substring(sampleStartIdx, sampleEndIdx);
        // for python plot
        TaxonQuantChartData taxonQuantChartData = new TaxonQuantChartData();
        taxonQuantChartData.rank = rank;
        taxonQuantChartData.log2 = log2;
        taxonQuantChartData.samples = Arrays.stream(samplesStr.split("\t")).collect(Collectors.toList());
        int sampleSize = taxonQuantChartData.samples.size();
        taxonQuantChartData.taxon2Quants = new HashMap<>();
        // tsv result file
        Path taxonQuantFilePath = Paths.get(taxonQuantFile);
        BufferedWriter bw = Files.newBufferedWriter(taxonQuantFilePath);
        bw.write(String.join(System.lineSeparator(),
                "#Rank: " + rank,
                "#Min peptide count for taxon : " + minPeptideForTaxon,
                log2 ? "# Enable Log2 transform" : "# Disable Log2 transform")
                + System.lineSeparator());
        bw.write(String.join("\t", "Taxon Id", "Taxon name",
                samplesStr, "Peptides") + System.lineSeparator());
        List<Double> totalTaxonQuantVals = IntStream.range(0, sampleSize).mapToObj(i -> 0.0)
                .collect(Collectors.toList());
        for (TaxonQuant taxonQuant : taxonQuantList) {
            List<Double> taxonQuantVals = IntStream.range(0, sampleSize).mapToObj(i -> 0.0)
                    .collect(Collectors.toList());
            for (QuantPeptide quantPeptide : taxonQuant.quantPeptideList) {
                List<Double> peptQunatVals = quantPeptide.getQuantValList();
                for (int j = 0; j < peptQunatVals.size(); j++) {
                    taxonQuantVals.set(j, taxonQuantVals.get(j) + peptQunatVals.get(j));
                    totalTaxonQuantVals.set(j, totalTaxonQuantVals.get(j) + peptQunatVals.get(j));
                }
            }
            taxonQuantChartData.taxon2Quants.put(taxonQuant.unipeptTaxon.getName(), taxonQuantVals);
            taxonQuant.quantValList = taxonQuantVals;
        }

        for (TaxonQuant taxonQuant : taxonQuantList) {
            List<Double> quantValList = IntStream.range(0, taxonQuant.quantValList.size())
                    .mapToObj(i -> log2
                            ? Math.log(taxonQuant.quantValList.get(i) / totalTaxonQuantVals.get(i)) / log2Val
                            : taxonQuant.quantValList.get(i) / totalTaxonQuantVals.get(i) * 100
                    )
                    .collect(Collectors.toList());
            UnipeptTaxon taxon = taxonQuant.unipeptTaxon;
            bw.write(String.join("\t", taxon.getId().toString(), taxon.getName(),
                    quantValList.stream()
                            .map(quant -> quant.toString())
                            .collect(Collectors.joining("\t")),
                    taxonQuant.quantPeptideList.stream()
                            .map(pept -> pept.getSequence())
                            .collect(Collectors.joining("\t"))
            ) + System.lineSeparator());
        }
        bw.close();
        // run python to plot
        Path taxonQuantJsonFilePath = Files.createTempFile("taxon-quant", ".json");
        objectMapper.writeValue(taxonQuantJsonFilePath.toFile(), taxonQuantChartData);
        Path taxonQuantChartPath = taxonQuantFilePath.getParent().resolve("taxon-quant-" + rank + ".png");
        String taxonQuantChartCmd = String.join(" ",
                python3,
                taxonQuantScript,
                taxonQuantJsonFilePath.toString(),
                taxonQuantChartPath.toString()
        );
        CommandExecutor.exec(taxonQuantChartCmd);
        Files.delete(taxonQuantJsonFilePath);
    }

    private Map<UnipeptTaxon, List<QuantPeptide>> lca2Quant(String lcaFile,
                                                           String rank) throws IOException {
        BufferedReader br = Files.newBufferedReader(Paths.get(lcaFile));
        // peptide file contains quant data?
        List<String> headerList = Arrays.stream(br.readLine().split("\t")).collect(Collectors.toList());
        int sampleStartIdx = 1;
        int sampleEndIdx = headerList.indexOf("taxon_id") - 1;
        int rankIdIdx = headerList.indexOf(rank + "_id");
        int rankNameIdx = headerList.indexOf(rank + "_name");
        boolean containQuant = sampleEndIdx == sampleEndIdx;

        /** stat quantitative **/
        Map<UnipeptTaxon, List<QuantPeptide>> taxon2QuantPeptides = new HashMap<>();
        String line;
        while ((line = br.readLine()) != null) {
            String[] tmp = Arrays.stream(line.split("\t", rankNameIdx + 2))
                    .map(String::trim)
                    .toArray(String[]::new);
            String rankId;
            if (tmp.length > rankIdIdx && (rankId = tmp[rankIdIdx]).length() > 0) {
                UnipeptTaxon taxon = UnipeptTaxon.builder()
                        .id(Integer.parseInt(rankId))
                        .name(tmp[rankNameIdx])
                        .build();
                List<QuantPeptide> quantPeptideList = taxon2QuantPeptides.get(taxon);
                if (quantPeptideList == null) {
                    quantPeptideList = new ArrayList<>();
                    taxon2QuantPeptides.put(taxon, quantPeptideList);
                }
                quantPeptideList.add(
                        QuantPeptide.builder()
                                .sequence(tmp[0])
                                .quantValList(containQuant ? IntStream.rangeClosed(sampleStartIdx, sampleEndIdx)
                                        .mapToObj(i -> tmp[i].length() == 0 ? 0 : Double.valueOf(tmp[i]))
                                        .collect(Collectors.toList())
                                        : null)
                                .build()
                );
            }
        }
        br.close();
        return taxon2QuantPeptides;
    }

    @Data
    @NoArgsConstructor
    @AllArgsConstructor
    class PeptideCount {
        private int peptideCountForRank;

        private int peptideCountForSubrank;
    }
}