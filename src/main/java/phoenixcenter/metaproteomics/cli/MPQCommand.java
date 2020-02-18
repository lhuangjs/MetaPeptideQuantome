package phoenixcenter.metaproteomics.cli;


import phoenixcenter.metaproteomics.TaxAnalysis;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.util.concurrent.Callable;

@Command(name = "MPQ", mixinStandardHelpOptions = true, version = "MPQ 1.0",
        description = "MetaPeptideQuantome")
public class MPQCommand implements Callable<Integer> {

    @Option(names = {"-i", "--input"}, required = true)
    private String peptideFile;

    @Option(names = {"--il"}, defaultValue = "true")
    private boolean equalIL;

    @Option(names = {"--mc"}, defaultValue = "false")
    private boolean missedCleavage;

    @Option(names = {"-o", "--output"})
    private String lcaFile;


    public static void main(String[] args) {
        int exitCode = new CommandLine(new MPQCommand()).execute(args);
        System.exit(exitCode);
    }

    @Override
    public Integer call() throws Exception {
        TaxAnalysis taxAnalysis = new TaxAnalysis();
        taxAnalysis.peptide2LCA(peptideFile, equalIL, missedCleavage, lcaFile);
        return 0;
    }
}
