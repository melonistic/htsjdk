package htsjdk.samtools;

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * If you would like to consider the library as the major sort key, please specify the headers when constructing this comparator.
 *
 * @author nhomer
 */
public class SAMRecordDuplicateComparator implements SAMRecordComparator {

    //private static final byte F = 0, R = 1, FF = 2, FR = 3, RR = 4, RF = 5;
    private static final byte FF = 0, FR = 1, F = 2, RF = 3, RR = 4, R = 5;

    private final Map<String, Short> libraryIds = new HashMap<String, Short>(); // from library string to library id
    private short nextLibraryId = 1;
    
    private ScoringStrategy scoringStrategy = ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH;
    
    public SAMRecordDuplicateComparator() {}

    public SAMRecordDuplicateComparator(final List<SAMFileHeader> headers) {
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.duplicate, headers, false);
        
        // pre-populate the library names
        for (final SAMFileHeader header : headers) {
            for (final SAMReadGroupRecord readGroup : header.getReadGroups()) {
                final String libraryName = readGroup.getLibrary();
                if (null != libraryName) {
                    final short libraryId = this.nextLibraryId++;
                    this.libraryIds.put(libraryName, libraryId);
                }
            }
        }
    }
    
    public void setScoringStrategy(final ScoringStrategy scoringStrategy) {
        this.scoringStrategy = scoringStrategy;
    }


    /**
     * Gets the library name from the header for the record. If the RG tag is not present on
     * the record, or the library isn't denoted on the read group, a constant string is
     * returned.
     */
    private static String getLibraryName(final SAMRecord rec) {
        final String readGroupId = (String) rec.getAttribute("RG");

        if (readGroupId != null) {
            final SAMReadGroupRecord rg = rec.getHeader().getReadGroup(readGroupId);
            if (rg != null) {
                final String libraryName = rg.getLibrary();
                if (null != libraryName) return libraryName;
            }
        }

        return "Unknown Library";
    }

    /** Get the library ID for the given SAM record. */
    private short getLibraryId(final SAMRecord rec) {
        final String library = getLibraryName(rec);
        Short libraryId = this.libraryIds.get(library);

        if (libraryId == null) {
            libraryId = this.nextLibraryId++;
            this.libraryIds.put(library, libraryId);
        }

        return libraryId;
    }
    
    private int compareOrientationByteCollapseOrientation(final int orientation1, final int orientation2) {
        // F == FR, F == FF
        // R == RF, R == RR
        if (F == orientation1 || R == orientation1) { // first orientation is fragment
            /**
             * We want 
             * F == FR, F == FF
             * R == RF, R == RR
             */
            if (F == orientation1) {
                if (F == orientation2 || FR == orientation2 || FF == orientation2) {
                    return 0;
                }
            }
            else { // R == orientation1
                if (R == orientation2 || RF == orientation2 || RR == orientation2) {
                    return 0;
                }
            }
        }
        else if (F == orientation2 || R == orientation2) { // first orientation is paired, second is fragment
            return -compareOrientationByteCollapseOrientation(orientation2, orientation1);
        }

        return orientation1 - orientation2;
    }
    
    /**
     * Returns a single byte that encodes the orientation of the two reads in a pair.
     */
    private static byte getPairedOrientationByte(final boolean read1NegativeStrand, final boolean read2NegativeStrand) {
        if (read1NegativeStrand) {
            if (read2NegativeStrand) return SAMRecordDuplicateComparator.RR;
            else return SAMRecordDuplicateComparator.RF;
        } else {
            if (read2NegativeStrand) return SAMRecordDuplicateComparator.FR;
            else return SAMRecordDuplicateComparator.FF;
        }
    }
    
    private int getFragmentOrientation(final SAMRecord record) {
         return record.getReadNegativeStrandFlag() ? SAMRecordDuplicateComparator.R : SAMRecordDuplicateComparator.F;
    }

    private int getPairedOrientation(final SAMRecord record) {
        if (record.getReadPairedFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) {
            return getPairedOrientationByte(record.getReadNegativeStrandFlag(), record.getMateNegativeStrandFlag());
        } else {
            return getFragmentOrientation(record);
        }
    }

    private int getMateReferenceIndex(final SAMRecord record) {
        if (record.getReadPairedFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) {
            return record.getMateReferenceIndex();
        } else {
            return -1;
        }
    }

    private int getMateCoordinate(final SAMRecord record) {
        if (record.getReadPairedFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) {
            return record.getMateNegativeStrandFlag() ? SAMUtils.getMateUnclippedEnd(record) : SAMUtils.getMateUnclippedStart(record);
        } else {
            return -1;
        }
    }
    
    private boolean hasUnmappedEnd(final SAMRecord record) {
        // HERE
        //System.err.println("hasUnmappedEnd: (" + record.getReadUnmappedFlag() + " || " + (record.getReadPairedFlag() && record.getMateUnmappedFlag()) + ")");
        //System.err.print(record.getSAMString());
        return (record.getReadUnmappedFlag() || (record.getReadPairedFlag() && record.getMateUnmappedFlag()));
    }

    private boolean hasMappedEnd(final SAMRecord record) {
        return (!record.getReadUnmappedFlag() || (record.getReadPairedFlag() && !record.getMateUnmappedFlag()));
    }
    
    private boolean pairedEndAndBothMapped(final SAMRecord record) {
        return (record.getReadPairedFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag());
        
    }
    
    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        int cmp;

        // temporary variables for comparisons
        int samRecord1Value, samRecord2Value;

        cmp = fileOrderCompare(samRecord1, samRecord2);
        // the duplicate scoring strategy
        if (cmp == 0) {
            cmp = DuplicateScoringStrategy.compare(samRecord1, samRecord2, this.scoringStrategy, true);
        }
        // the read name
        if (cmp == 0) {
            cmp = samRecord1.getReadName().compareTo(samRecord2.getReadName());
        }
        /**
         * If both reads are paired and both ends mapped, always prefer the first end over the second end.  This is needed to
         * properly choose the first end for optical duplicate identification when both ends are mapped to the same position etc.
         */
        if (cmp == 0) {
            if (samRecord1.getReadPairedFlag() && samRecord2.getReadPairedFlag()) {
                samRecord1Value = samRecord1.getFirstOfPairFlag() ? 0 : 1;
                samRecord2Value = samRecord2.getFirstOfPairFlag() ? 0 : 1;
                cmp = samRecord1Value - samRecord2Value;
            }
        }

        return cmp;
    }

    private int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2, final boolean collapseOrientation, final boolean collapseSingleEndMapped) {
        int cmp;

        // temporary variables for comparisons
        int samRecord1Value, samRecord2Value;

        // library identifier
        {
            samRecord1Value = getLibraryId(samRecord1);
            samRecord2Value = getLibraryId(samRecord2);
            cmp = samRecord1Value - samRecord2Value;
        }
        // reference index
        if (cmp == 0) {
            samRecord1Value = samRecord1.getReferenceIndex();
            samRecord2Value = samRecord2.getReferenceIndex();
            // NB: this accounts for unmapped reads to be placed at the ends of the file
            if (samRecord1Value == -1) {
                cmp = (samRecord2Value == -1) ? 0 : 1;
            }
            else if (samRecord2Value == -1) {
                cmp = -1;
            }
            else {
                cmp = samRecord1Value - samRecord2Value;
            }
        }
        // read coordinate
        if (cmp == 0) {
            samRecord1Value = samRecord1.getReadNegativeStrandFlag() ? samRecord1.getUnclippedEnd() : samRecord1.getUnclippedStart();
            samRecord2Value = samRecord2.getReadNegativeStrandFlag() ? samRecord2.getUnclippedEnd() : samRecord2.getUnclippedStart();
            cmp = samRecord1Value - samRecord2Value;
        }
        // orientation
        if (cmp == 0) {
            samRecord1Value = getPairedOrientation(samRecord1);
            samRecord2Value = getPairedOrientation(samRecord2);
            if (collapseOrientation) {
                cmp = compareOrientationByteCollapseOrientation(samRecord1Value, samRecord2Value);
            }
            else {
                cmp = samRecord1Value - samRecord2Value;
            }
        }
        // both ends need to be mapped
        if (pairedEndAndBothMapped(samRecord1) && pairedEndAndBothMapped(samRecord2)) {
            // mate's reference index
            if (cmp == 0) {
                samRecord1Value = getMateReferenceIndex(samRecord1);
                samRecord2Value = getMateReferenceIndex(samRecord2);
                cmp = samRecord1Value - samRecord2Value;
            }
            // mate's coordinate
            if (cmp == 0) {
                samRecord1Value = getMateCoordinate(samRecord1);
                samRecord2Value = getMateCoordinate(samRecord2);
                cmp = samRecord1Value - samRecord2Value;
            }
        }
        if (cmp == 0) {
            samRecord1Value = hasMappedEnd(samRecord1) ? 0 : 1;
            samRecord2Value = hasMappedEnd(samRecord2) ? 0 : 1;
            cmp = samRecord1Value - samRecord2Value;
        }
        // if both paired or both unpaired, then check if one of the two ends (or single end) is unmapped
        // else prefer the one that is paired end
        if (cmp == 0 && !collapseSingleEndMapped) {
            if (samRecord1.getReadPairedFlag() == samRecord2.getReadPairedFlag()) {
                // Is this unmapped or its mate?
                samRecord1Value = hasUnmappedEnd(samRecord1) ? 1 : 0;
                samRecord2Value = hasUnmappedEnd(samRecord2) ? 1 : 0;
                cmp = samRecord1Value - samRecord2Value;
            }
            else { // if we care if one is paired and the other is not
                cmp = samRecord1.getReadPairedFlag() ? -1 : 1;
            }
        }

        return cmp;
    }

    /**
     * Compare for equivalence when creating a set of duplicates, not discriminating between records that are duplicates of each other.
     * 
     * Major difference between this and fileOrderCompare is how we compare the orientation byte.  Here we want:
     *   F == FR, F == FF
     *   R == RF, R == RR
     */
    public int duplicateSetCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        return fileOrderCompare(samRecord1, samRecord2, true, true);
    }

    /**
     *
     * @param samRecord1
     * @param samRecord2
     * @return
     */
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        return fileOrderCompare(samRecord1, samRecord2, false, false);
    }
}
