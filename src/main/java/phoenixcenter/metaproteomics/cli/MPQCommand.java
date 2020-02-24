package phoenixcenter.metaproteomics.cli;


import phoenixcenter.metaproteomics.Mascot2XMLEnhancer;
import phoenixcenter.metaproteomics.PeptideProphetEnhancer;
import phoenixcenter.metaproteomics.TaxAnalysis;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;

@Command(name = "MPQ", mixinStandardHelpOptions = true, version = "MPQ 1.0",
        description = "MetaPeptideQuantome")
public class MPQCommand implements Callable<Integer> {

    private Mascot2XMLEnhancer mascot2XMLEnhancer = new Mascot2XMLEnhancer();

    private PeptideProphetEnhancer peptideProphetEnhancer = new PeptideProphetEnhancer();

    private TaxAnalysis taxAnalysis = new TaxAnalysis();

    @Command(name = "target2decoy", description = "Generate a new library contains target sequences and reversed decoy sequences")
    public void genTargetDecoyLibrary(
            @Option(names = "-i", description = "Path of the target library file", required = true) String targetLibrary,
            @Option(names = "-P", description = "The prefix of decoy protein sequence", required = true) String decoyPrefix,
            @Option(names = "-o", description = "Path of the output target-decoy database file", required = true) String tdLibrary
    ) throws IOException {
        mascot2XMLEnhancer.genLibraryWithDecoy(targetLibrary, decoyPrefix, tdLibrary);
    }

    @Command
    public void dat2xml(
            @Option(names = "-i", description = "Path of the input DAT file", required = true) String datFile,
            @Option(names = "-D", description = "Sequence database path", required = true) String targetDecoyLibrary,
            @Option(names = "-E", description = "Enzyme", required = true) String enzyme,
            @Option(names = "-M", defaultValue = "false",
                    description = "Whether it is a Dat file generated in an manual decoy database search mode, default is false") boolean isManual,
            @Option(names = "-P", description = "The prefix of decoy protein sequence, only required in in automatic decoy database search mode") String decoyPrefix,
            @Option(names = "-o", description = "Path of the output pepxm file", required = true) String pepxmlFile) throws IOException {
        if (isManual) {
            mascot2XMLEnhancer.convertManualDecoyMode(datFile, enzyme, targetDecoyLibrary, pepxmlFile);
        } else {
            mascot2XMLEnhancer.convertAutoDecoyMode(datFile, decoyPrefix, enzyme, targetDecoyLibrary, pepxmlFile);
        }
    }

    @Command(name = "xml2pp", description = "Run PeptideProphet")
    public void xml2pp(
            @Option(names = "-i", description = "Path of the input pepXML file", required = true) String pepxmlFile,
            @Option(names = "-FDR", description = "The peptide-level FDR") Double fdr,
            @Option(names = "-P", description = "The prefix of decoy protein sequence, required when filtering by FDR") String decoyPrefix,
            @Option(names = "-D", description = "The sequence database", required = true) String libraryFile,
            @Option(names = "-PPM ", defaultValue = "true",
                    description = "Use PPM instead of daltons in Accurate Mass Model, default true") boolean ppm,
            @Option(names = "-o", description = "Path of the output pepXML file", required = true) String ppPepxmlFile
    ) throws IOException, FileParsingException {
        Map<String, String> params = new HashMap<>();
        params.put("-D", libraryFile);
        if (ppm) {
            params.put("-PPM", null);
        }
        if (fdr != null && fdr < 1) {
            Path tmpPath = Files.createTempFile("pp", ".pep.xml");
            peptideProphetEnhancer.runPeptideProphet(pepxmlFile, params, tmpPath.toString());
            peptideProphetEnhancer.filterByFDR(tmpPath, decoyPrefix, fdr, Paths.get(ppPepxmlFile));
            Files.delete(tmpPath);
        } else {
            peptideProphetEnhancer.runPeptideProphet(pepxmlFile, params, ppPepxmlFile);
        }
    }

    @Command(name = "xml2tsv", description = "Extract PSM count from pepXML files")
    public void xml2tsv(
            @Option(names = "-i", description = "Path of the input pepXML files, separated by commas", required = true) String pepxmlFiles,
            @Option(names = "-FDR", description = "The peptide-level FDR") Double fdr,
            @Option(names = "-P", description = "The prefix of decoy protein sequence") String decoyPrefix,
            @Option(names = "-o", description = "Path of the output tsv file", required = true) String tsvFile
    ) throws IOException, FileParsingException {
        if (fdr == null) {
            peptideProphetEnhancer.pp2tsv(decoyPrefix, pepxmlFiles.split(","), tsvFile);
        } else {
            peptideProphetEnhancer.pp2tsv(decoyPrefix, pepxmlFiles.split(","),
                    fdr, tsvFile);
        }
    }

    @Command(name = "pept2lca", description = "Retrieve LCA information from Unipept")
    public void pept2lca(
            @Option(names = "-i", description = "Peptide list file path", required = true) String peptideFile,
            @Option(names = "-il", defaultValue = "true",
                    description = "Equate I and L, default true") boolean equateIL,
            @Option(names = "-m", defaultValue = "false",
                    description = "Advanced missing cleavage handling, default false") boolean missedCleavage,
            @Option(names = "-o", description = "The output LCA file path", required = true) String lcaFile
    ) throws IOException {
        taxAnalysis.peptide2LCA(peptideFile, equateIL, missedCleavage, lcaFile);
    }

    @Command(name = "lca2quant", description = "Taxonomic quantative analysis")
    public void lca2quant(
            @Option(names = "-i", description = "The LCA file path", required = true) String lcaFile,
            @Option(names = "-rank", description = "forma, varietas, subspecies, species, species_subgroup, species_group, subgenus, genus, subtribe, " +
                    "tribe, subfamily, family, superfamily, parvorder, infraorder, suborder, order, superorder, " +
                    "infraclass, subclass, class, superclass, subphylum, phylum, superphylum, subkingdom, kingdom," +
                    "superkingdom") String rank,
            @Option(names = "-log2", defaultValue = "false",
                    description = "Log2 transform, default false") boolean log2,
            @Option(names = "-min", description = "Min peptide count for each taxon") int minPeptideForTaxon,
            @Option(names = "-o", description = "The output path of taxonomic quantative analysis file", required = true) String taxonQuantFile
    ) throws IOException {
        taxAnalysis.lca2Quant(lcaFile, rank, minPeptideForTaxon, log2, taxonQuantFile);
    }


    public static void main(String[] args) {
        int exitCode = new CommandLine(new MPQCommand()).execute(args);
        System.exit(exitCode);
    }

    @Override
    public Integer call() throws Exception {

        return 0;
    }
}
