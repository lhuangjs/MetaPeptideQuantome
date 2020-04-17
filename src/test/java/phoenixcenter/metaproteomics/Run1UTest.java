package phoenixcenter.metaproteomics;

import org.junit.Test;

import java.io.IOException;

public class Run1UTest {

    TaxAnalysis taxAnalysis = new TaxAnalysis();

    @Test
    public void peptide2LCA() throws IOException {
        String peptideWithQuant = "peptide.tsv";
        String lcaFile = "peptide.tsv";
        taxAnalysis.peptide2LCA(peptideWithQuant, true, true, lcaFile);
    }
}
