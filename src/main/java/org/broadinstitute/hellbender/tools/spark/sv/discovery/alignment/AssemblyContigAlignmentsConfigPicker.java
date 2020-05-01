package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.Serializable;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;

/**
 * A simple heuristic optimizer based on extensive manual review of alignments
 * produced by the aligner (currently "bwa mem -x intractg") with the aim for
 * picking a configuration that provides "optimal coverage" for the input
 * assembly contig.
 * Future improvements are definitely welcome and could benefit the whole pipeline.
 */
public class AssemblyContigAlignmentsConfigPicker {

    /**
     * Filters input alignments of single-ended long reads, e.g. local assembly contigs,
     * with the objective of
     * choosing a set of alignments that provide "optimal coverage" of the assembly contig.
     *
     * Currently "optimality" is defined based on an heuristic scoring scheme
     * {@link AlignedContig#computeScoreOfConfiguration(List, Set, int)}.
     *
     * <p>
     *     It goes through four major steps:
     *     <ul>
     *         <li>
     *             first group the alignments of each local assembly contig together
     *         </li>
     *         <li>
     *             then exhaustively iterate through all possible combinations (named configuration) of the alignments,
     *             score them, and pick the best ones with the "optimal coverage"; note that sometimes
     *             there are multiple "optimal" configurations, and when this happens, all of them are returned,
     *             with each of them having the field {@link AssemblyContigWithFineTunedAlignments#hasEquallyGoodAlnConfigurations}
     *             set to true.
     *         </li>
     *         <li>
     *             alignments containing large gaps (i.e. insertions and deletions)
     *             are finally split at gap start and end locations.
     *             (note, see warning in {@link AlignedContig#splitGaps(AlignedContig.GoodAndBadMappings, boolean)})
     *         </li>
     *         <li>
     *             A final round of pruning, in order to further remove uninformative alignments.
     *         </li>
     *     </ul>
     * </p>
     *
     * @param assemblyAlignments    long read alignments
     * @param header                header for the long reads
     * @param canonicalChromosomes  a set of chromosome names that are defined as canonical, e.g. for Human, chr1-chr22, and chrX and chrY
     * @param scoreDiffTolerance    a tolerance where if two configurations' scores differ by less than or equal to this amount, they are considered equally good
     * @param toolLogger            logger for outputting summary and debugging
     *
     * @return              contigs with alignments filtered and custom formatted as {@link AlignmentInterval}
     */
    public static JavaRDD<AssemblyContigWithFineTunedAlignments> createOptimalCoverageAlignmentSetsForContigs(final JavaRDD<GATKRead> assemblyAlignments,
                                                                                                              final SAMFileHeader header,
                                                                                                              final Set<String> canonicalChromosomes,
                                                                                                              final double scoreDiffTolerance,
                                                                                                              final Logger toolLogger) {
        final JavaRDD<AlignedContig> parsedContigAlignments =
                convertRawAlignmentsToAlignedContigAndFilterByQuality(assemblyAlignments, header, toolLogger);

        return parsedContigAlignments
                .flatMap(tig ->
                        tig.reConstructContigFromPickedConfiguration(canonicalChromosomes, scoreDiffTolerance).iterator());
    }

    // step 1: parse and primitive filter ==============================================================================

    /**
     * Parses input alignments into custom {@link AlignmentInterval} format, and
     * performs a primitive filtering implemented in
     * {@link AlignedContig#notDiscardForBadMQ()} that
     * gets rid of contigs with no good alignments.
     *
     * It's important to remember that this step doesn't select alignments,
     * but only parses alignments and either keeps the whole contig or drops it completely.
     */
    private static JavaRDD<AlignedContig> convertRawAlignmentsToAlignedContigAndFilterByQuality(final JavaRDD<GATKRead> assemblyAlignments,
                                                                                                final SAMFileHeader header,
                                                                                                final Logger toolLogger) {
        assemblyAlignments.cache();
        toolLogger.info( "Processing " + assemblyAlignments.count() + " raw alignments from " +
                         assemblyAlignments.map(GATKRead::getName).distinct().count() + " contigs.");

        final JavaRDD<AlignedContig> parsedContigAlignments =
                new SAMFormattedContigAlignmentParser(assemblyAlignments, header, false)
                        .getAlignedContigs()
                        .filter(AlignedContig::notDiscardForBadMQ).cache();
        assemblyAlignments.unpersist();
        toolLogger.info( "Filtering on MQ left " + parsedContigAlignments.count() + " contigs.");
        return parsedContigAlignments;
    }

    public static final class SAMFormattedContigAlignmentParser extends AlignedContigGenerator implements Serializable {
        private static final long serialVersionUID = 1L;

        private final JavaRDD<GATKRead> unfilteredContigAlignments;
        private final SAMFileHeader header;
        private final boolean splitGapped;

        public SAMFormattedContigAlignmentParser(final JavaRDD<GATKRead> unfilteredContigAlignments,
                                                 final SAMFileHeader header, final boolean splitGapped) {
            this.unfilteredContigAlignments = unfilteredContigAlignments;
            this.header = header;
            this.splitGapped = splitGapped;
        }

        @Override
        public JavaRDD<AlignedContig> getAlignedContigs() {
            return unfilteredContigAlignments
                    .filter(r -> !r.isSecondaryAlignment())
                    .groupBy(GATKRead::getName)
                    .map(Tuple2::_2)
                    .map(gatkReads ->
                            parseReadsAndOptionallySplitGappedAlignments(
                                    Utils.stream(gatkReads).map(r->r.convertToSAMRecord(header)).collect(Collectors.toList()),
                                    GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, splitGapped));
        }

        /**
         * Iterates through the input {@code noSecondaryAlignments}, which are assumed to contain no secondary alignment (i.e. records with "XA" tag),
         * converts to custom {@link AlignmentInterval} format and
         * split the records when the gap in the alignment reaches the specified {@code sensitivity}.
         * The size of the returned iterable of {@link AlignmentInterval}'s is guaranteed to be no lower than that of the input iterable.
         */
        @VisibleForTesting
        public static AlignedContig parseReadsAndOptionallySplitGappedAlignments(final Iterable<SAMRecord> noSecondaryAlignments,
                                                                                 final int gapSplitSensitivity,
                                                                                 final boolean splitGapped) {

            Utils.validateArg(noSecondaryAlignments.iterator().hasNext(), "input collection of GATK reads is empty");

            final SAMRecord primaryAlignment
                    = Utils.stream(noSecondaryAlignments).filter(sam -> !sam.getSupplementaryAlignmentFlag())
                    .findFirst()
                    .orElseThrow(() -> new GATKException("no primary alignment for read " + noSecondaryAlignments.iterator().next().getReadName()));

            Utils.validate(!primaryAlignment.getCigar().containsOperator(CigarOperator.H),
                    "assumption that primary alignment does not contain hard clipping is invalid for read: " + primaryAlignment.toString());

            final byte[] contigSequence = primaryAlignment.getReadBases().clone();
            final List<AlignmentInterval> parsedAlignments;
            if ( primaryAlignment.getReadUnmappedFlag() ) { // the Cigar
                parsedAlignments = Collections.emptyList();
            } else {
                if (primaryAlignment.getReadNegativeStrandFlag()) {
                    SequenceUtil.reverseComplement(contigSequence);
                }

                final Stream<AlignmentInterval> unSplitAIList = Utils.stream(noSecondaryAlignments).map(AlignmentInterval::new);
                if (splitGapped) {
                    final int unClippedContigLength = primaryAlignment.getReadLength();
                    parsedAlignments = unSplitAIList.map(ar ->
                            ContigAlignmentsModifier.splitGappedAlignment(ar, gapSplitSensitivity, unClippedContigLength))
                            .flatMap(Utils::stream).collect(Collectors.toList());
                } else {
                    parsedAlignments = unSplitAIList.collect(Collectors.toList());
                }
            }
            return new AlignedContig(primaryAlignment.getReadName(), contigSequence, parsedAlignments);
        }
    }
}
