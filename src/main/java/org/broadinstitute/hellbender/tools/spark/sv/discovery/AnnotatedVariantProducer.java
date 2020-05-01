package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.Serializable;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigAlignmentsSparkArgumentCollection;

/**
 * Given identified pair of breakpoints for a simple SV and its supportive evidence, i.e. chimeric alignments,
 * produce an annotated {@link VariantContext}.
 */
public class AnnotatedVariantProducer implements Serializable {
    private static final long serialVersionUID = 1L;


    public static VariantContext produceAnnotatedVcFromEvidenceTargetLink(final EvidenceTargetLink evidenceTargetLink,
                                                                          final SvType svType) {
        final PairedStrandedIntervals pairedStrandedIntervals = evidenceTargetLink.getPairedStrandedIntervals();
        final StrandedInterval strandedIntervalLeft = pairedStrandedIntervals.getLeft();
        final StrandedInterval strandedIntervalRight = pairedStrandedIntervals.getRight();
        final int start = strandedIntervalLeft.getInterval().midpoint();
        final int end = strandedIntervalRight.getInterval().midpoint();
        final VariantContextBuilder builder = svType
                .getBasicInformation()
                .attribute(GATKSVVCFConstants.CIPOS, produceCIInterval(start, strandedIntervalLeft.getInterval()))
                .attribute(GATKSVVCFConstants.CIEND, produceCIInterval(end, strandedIntervalRight.getInterval()))
                .attribute(GATKSVVCFConstants.READ_PAIR_SUPPORT, evidenceTargetLink.getReadPairs())
                .attribute(GATKSVVCFConstants.SPLIT_READ_SUPPORT, evidenceTargetLink.getSplitReads());
        return builder.make();
    }

    public static List<VariantContext> annotateBreakpointBasedCallsWithImpreciseEvidenceLinks(final List<VariantContext> assemblyDiscoveredVariants,
                                                                                              final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                                                              final ReadMetadata metadata,
                                                                                              final ReferenceMultiSparkSource reference,
                                                                                              final DiscoverVariantsFromContigAlignmentsSparkArgumentCollection parameters,
                                                                                              final Logger localLogger) {

        final int originalEvidenceLinkSize = evidenceTargetLinks.size();
        final List<VariantContext> result = assemblyDiscoveredVariants
                .stream()
                .map(variant -> annotateWithImpreciseEvidenceLinks(
                        variant,
                        evidenceTargetLinks,
                        reference.getReferenceSequenceDictionary(null),
                        metadata, parameters.assemblyImpreciseEvidenceOverlapUncertainty))
                .collect(Collectors.toList());
        localLogger.info("Used " + (originalEvidenceLinkSize - evidenceTargetLinks.size()) + " evidence target links to annotate assembled breakpoints");
        return result;
    }

    //==================================================================================================================

    private static VariantContext annotateWithImpreciseEvidenceLinks(final VariantContext variant,
                                                                     final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                                                     final SAMSequenceDictionary referenceSequenceDictionary,
                                                                     final ReadMetadata metadata,
                                                                     final int defaultUncertainty) {
        if (variant.getStructuralVariantType() == StructuralVariantType.DEL) {
            SVContext svc = SVContext.of(variant);
            final int padding = (metadata == null) ? defaultUncertainty : (metadata.getMaxMedianFragmentSize() / 2);
            PairedStrandedIntervals svcIntervals = svc.getPairedStrandedIntervals(metadata, referenceSequenceDictionary, padding);

            final Iterator<Tuple2<PairedStrandedIntervals, EvidenceTargetLink>> overlappers = evidenceTargetLinks.overlappers(svcIntervals);
            int readPairs = 0;
            int splitReads = 0;
            while (overlappers.hasNext()) {
                final Tuple2<PairedStrandedIntervals, EvidenceTargetLink> next = overlappers.next();
                readPairs += next._2.getReadPairs();
                splitReads += next._2.getSplitReads();
                overlappers.remove();
            }
            final VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variant);
            if (readPairs > 0) {
                variantContextBuilder.attribute(GATKSVVCFConstants.READ_PAIR_SUPPORT, readPairs);
            }
            if (splitReads > 0) {
                variantContextBuilder.attribute(GATKSVVCFConstants.SPLIT_READ_SUPPORT, splitReads);
            }

            return variantContextBuilder.make();
        } else {
            return variant;
        }
    }

    /**
     * Produces the string representation of a VCF 4.2-style SV CI interval centered around 'point'.
     */
    @VisibleForTesting
    static String produceCIInterval(final int point, final SVInterval ciInterval) {
        Utils.validate(ciInterval.getStart() <= point && ciInterval.getEnd() >= point, "Interval must contain point");
        return String.join(",",
                String.valueOf(ciInterval.getStart() - point),
                String.valueOf(ciInterval.getEnd() - point));
    }
}
