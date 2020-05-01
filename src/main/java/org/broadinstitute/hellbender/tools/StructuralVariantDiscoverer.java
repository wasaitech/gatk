package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.BasicReference;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments.AlignmentSignatureBasicType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleChimera;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.SimpleNovelAdjacencyAndChimericAlignmentEvidence;
import org.broadinstitute.hellbender.tools.spark.sv.utils.CNVInputReader;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

public class StructuralVariantDiscoverer extends ReadWalker {
    @Argument(doc = "External CNV calls file. Should be single sample VCF, and contain only confident autosomal non-reference CNV calls (for now).",
            fullName = "cnv-calls", optional = true)
    private String cnvCallsFile;

    @Argument(doc = "file containing non-canonical chromosome names (e.g chrUn_KI270588v1) in the reference, human reference (hg19 or hg38) assumed when omitted",
            shortName = "alt-tigs",
            fullName = "non-canonical-contig-names-file", optional = true)
    private String nonCanonicalChromosomeNamesFile;

    private static final double SCORE_DIFF_TOLERANCE = 0.;

    private String sampleId;
    private SAMSequenceDictionary refDict;
    private BasicReference reference;
    private SVIntervalTree<VariantContext> cnvCalls;
    private Set<String> canonicalChromosomes;

    private String currentContigName = null;
    private final List<GATKRead> readsForCurrentContig = new ArrayList<>();
    private final Map<NovelAdjacencyAndAltHaplotype, SimpleNovelAdjacencyAndChimericAlignmentEvidence> simpleMap = new HashMap<>();

    @Override public boolean requiresReads() { return true; }
    @Override public boolean requiresReference() { return true; }

    @Override public List<ReadFilter> getDefaultReadFilters() {
        return Arrays.asList(ReadFilterLibrary.MAPPED, ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
    }

    @Override public void onTraversalStart() {
        final SAMFileHeader header = getHeaderForReads();
        if ( header.getSortOrder() != SAMFileHeader.SortOrder.queryname ) {
            throw new UserException("This tool requires a queryname-sorted source of reads.");
        }
        sampleId = SVUtils.getSampleId(header);
        refDict = header.getSequenceDictionary();
        cnvCalls = cnvCallsFile == null ? null : CNVInputReader.loadCNVCalls(cnvCallsFile, header);
        canonicalChromosomes = SVUtils.getCanonicalChromosomes(nonCanonicalChromosomeNamesFile, refDict);
    }

    @Override public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        reference = referenceContext;
        final String contigName = read.getName();
        if ( !contigName.equals(currentContigName) ) {
            if ( !readsForCurrentContig.isEmpty() ) {
                processContigAlignments(readsForCurrentContig);
                readsForCurrentContig.clear();
            }
            currentContigName = contigName;
        }
        readsForCurrentContig.add(read);
    }

    @Override public Object onTraversalSuccess() {
        final Object result = super.onTraversalSuccess();
        final List<VariantContext> variants = new ArrayList<>(2 * simpleMap.size());
        for ( final SimpleNovelAdjacencyAndChimericAlignmentEvidence novelAdjacencyAndEvidence : simpleMap.values() ) {
            final List<SvType> svTypes =
                    novelAdjacencyAndEvidence.getNovelAdjacencyReferenceLocations().toSimpleOrBNDTypes(reference);
            variants.addAll(novelAdjacencyAndEvidence.turnIntoVariantContexts(svTypes, sampleId, refDict, cnvCalls));
        }
        return result;
    }

    private void processContigAlignments( final List<GATKRead> contigAlignments ) {
        final List<AlignmentInterval> alignmentIntervals = new ArrayList<>(contigAlignments.size());
        String contigName = null;
        byte[] contigSequence = null;
        for ( final GATKRead read : contigAlignments ) {
            contigName = read.getName();
            if ( !read.isSupplementaryAlignment() ) contigSequence = read.getBasesNoCopy();
            alignmentIntervals.add(new AlignmentInterval(read));
        }
        if ( contigSequence == null ) {
            throw new UserException("No primary line for " + contigName);
        }
        final AlignedContig alignedContig = new AlignedContig(contigName, contigSequence, alignmentIntervals);
        if ( !alignedContig.notDiscardForBadMQ() ) return;

        final List<AssemblyContigWithFineTunedAlignments> fineTunedAlignmentsList =
            alignedContig.reConstructContigFromPickedConfiguration(canonicalChromosomes, SCORE_DIFF_TOLERANCE);
        for ( final AssemblyContigWithFineTunedAlignments fineTunedAlignment : fineTunedAlignmentsList ) {
            if ( fineTunedAlignment.getAlignmentSignatureBasicType() == AlignmentSignatureBasicType.SIMPLE_CHIMERA ) {
                if ( SimpleChimera.splitPairStrongEnoughEvidenceForCA(fineTunedAlignment.getHeadAlignment(),
                                                                      fineTunedAlignment.getTailAlignment()) ) {
                    final SimpleChimera simpleChimera = fineTunedAlignment.extractSimpleChimera(refDict);
                    final NovelAdjacencyAndAltHaplotype novelAdjacency =
                            new NovelAdjacencyAndAltHaplotype(simpleChimera, fineTunedAlignment.getContigSequence(), refDict);
                    SimpleNovelAdjacencyAndChimericAlignmentEvidence val = simpleMap.get(novelAdjacency);
                    if ( val != null ) {
                        val.getAlignmentEvidence().add(simpleChimera);
                    } else {
                        final SimpleNovelAdjacencyAndChimericAlignmentEvidence novelAdjacencyEvidence =
                                new SimpleNovelAdjacencyAndChimericAlignmentEvidence(novelAdjacency,
                                                                Collections.singletonList(simpleChimera));
                        simpleMap.put(novelAdjacency, novelAdjacencyEvidence);
                    }
                }
            }
        }
    }
}
