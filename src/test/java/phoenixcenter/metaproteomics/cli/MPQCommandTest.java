package phoenixcenter.metaproteomics.cli;

import org.junit.Test;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.util.Arrays;
import java.util.stream.Collectors;

public class MPQCommandTest {

    MPQCommand mpqCommand = new MPQCommand();

    @Test
    public void genTargetDecoyLibrary() throws IOException {
        /*
         MPQ target2decoy -i /home/huangjs/Documents/mpq/test/biomass.fasta \
         -P DECOY_ \
         -o /home/huangjs/Documents/mpq/test/biomass.td.fasta
         */
        String targetLibrary = "/home/huangjs/Documents/mpq/test/biomass.fasta";
        String decoyPrefix = "DECOY_";
        String tdLibrary = "/home/huangjs/Documents/mpq/test/biomass.td.fasta";
        mpqCommand.genTargetDecoyLibrary(targetLibrary, decoyPrefix, tdLibrary);
    }

    @Test
    public void dat2xml() throws IOException {
        /*
         MPQ dat2xml -i /home/huangjs/Documents/mpq/test/F008283.dat \
         -P DECOY_ \
         -D /home/huangjs/Documents/mpq/test/biomass.td.fasta \
         -E trypsin \
         -P DECOY_ \
         -o /home/huangjs/Documents/mpq/test/F008283.pep.xml
         */
        String datFile = "/home/huangjs/Documents/mpq/test/F008283.dat";
        String targetDecoyLibrary = "/home/huangjs/Documents/mpq/test/biomass.td.fasta";
        String enzyme = "trypsin";
        String decoyPrefix = "DECOY_";
        String pepxmlFile = "/home/huangjs/Documents/mpq/test/F008283.pep.xml";
        mpqCommand.dat2xml(datFile, targetDecoyLibrary, enzyme, false, decoyPrefix, pepxmlFile);
    }

    @Test
    public void xml2pp() throws IOException, FileParsingException {
         /*
         MPQ xml2pp -i /home/huangjs/Documents/mpq/test/F008283.pep.xml \
         -P DECOY_ \
         -D /home/huangjs/Documents/mpq/test/biomass.fasta \
         -o /home/huangjs/Documents/mpq/test/F008283.pp.pep.xml
         */
        String pepxmlFile = "/home/huangjs/Documents/mpq/test/F008283.pep.xml";
        String decoyPrefix = "DECOY_";
        String libraryFile = "/home/huangjs/Documents/mpq/test/biomass.fasta";
        String ppPepXML = "/home/huangjs/Documents/mpq/test/F008283.pp.pep.xml";
        mpqCommand.xml2pp(pepxmlFile, null, decoyPrefix, libraryFile, true, ppPepXML);
    }

    @Test
    public void xml2tsv() throws IOException, FileParsingException {
        /*
         MPQ xml2tsv -i /home/huangjs/Documents/mpq/test/F008283.pp.pep.xml,\
         /home/huangjs/Documents/mpq/test/F008284.pp.pep.xml,\
         /home/huangjs/Documents/mpq/test/F008285.pp.pep.xml,\
         /home/huangjs/Documents/mpq/test/F008286.pp.pep.xml \
         -fdr 0.01 \
         -P DECOY_ \
         -o /home/huangjs/Documents/mpq/test/peptide.tsv
         */
        String[] pepxmlFiles = {
                "/home/huangjs/Documents/mpq/test/F008283.pp.pep.xml",
                "/home/huangjs/Documents/mpq/test/F008284.pp.pep.xml",
                "/home/huangjs/Documents/mpq/test/F008285.pp.pep.xml",
                "/home/huangjs/Documents/mpq/test/F008286.pp.pep.xml"
        };
        Double fdr = 0.01;
        String decoyPrefix = "DECOY_";
        String tsvFile = "/home/huangjs/Documents/mpq/test/peptide.tsv";
        mpqCommand.xml2tsv(Arrays.stream(pepxmlFiles).collect(Collectors.joining(",")), fdr, decoyPrefix, tsvFile);
    }

    @Test
    public void pept2lca() throws IOException {
        /*
         MPQ pept2lca -i /home/huangjs/Documents/mpq/test/peptide.tsv \
         -m \
         -o /home/huangjs/Documents/mpq/test/lca.tsv
         */
        String peptideFile = "/home/huangjs/Documents/mpq/test/peptide.tsv";
        boolean equateIL = true;
        boolean missedCleavage = true;
        String lcaFile = "/home/huangjs/Documents/mpq/test/lca.tsv";
        mpqCommand.pept2lca(peptideFile, equateIL, missedCleavage, lcaFile);
    }

    @Test
    public void lca2quant() throws IOException {
        /*
         MPQ lca2quant -i /home/huangjs/Documents/mpq/test/lca.tsv \
         -rank species \
         -E trypsin \
         -min 3 \
         -o /home/huangjs/Documents/mpq/test/taxon-quant-species.tsv
         */
        String lcaFile = "/home/huangjs/Documents/mpq/test/lca.tsv";
        String rank = "species";
        int minPeptideForTaxon = 3;
        boolean log2 = false;
        String taxonQuantFile = "/home/huangjs/Documents/mpq/test/taxon-quant-species.tsv";
        mpqCommand.lca2quant(lcaFile, rank, log2, minPeptideForTaxon, taxonQuantFile);
    }
}