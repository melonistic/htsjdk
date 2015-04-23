package htsjdk.samtools;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;
import java.util.Arrays;

/**
 * Make singletons for unmapped and secondary or supplementary.
 *
 * @author nhomer
 */
public class DuplicateSetIterator implements CloseableIterator<DuplicateSet> {
    
    private final CloseableIterator<SAMRecord> wrappedIterator;
    
    private DuplicateSet duplicateSet = null;

    // TODO: expose these
    private final int maxRecordsInRam = SAMFileWriterImpl.getDefaultMaxRecordsInRam();
    private final File tmpDir = new File(System.getProperty("java.io.tmpdir"));

    private final SAMRecordDuplicateComparator comparator;

    public DuplicateSetIterator(final CloseableIterator<SAMRecord> iterator, final SAMFileHeader header) {
        this(iterator, header, false);
    }

    public DuplicateSetIterator(final CloseableIterator<SAMRecord> iterator, final SAMFileHeader header, final boolean preSorted) {
        this.comparator = new SAMRecordDuplicateComparator(Arrays.asList(header));

        if (preSorted) {
            this.wrappedIterator = iterator;
        }
        else {
            // Sort it!
            final SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
                    new BAMRecordCodec(header), comparator,
                    maxRecordsInRam, tmpDir);

            while (iterator.hasNext()) {
                final SAMRecord record = iterator.next();
                alignmentSorter.add(record);
            }
            iterator.close();
            
            this.wrappedIterator = alignmentSorter.iterator();
        }
        
        this.duplicateSet = new DuplicateSet();
        
        
        if (hasNext()) {
            final SAMRecord record = this.wrappedIterator.next();
            //System.err.print("DuplicateSet::DuplicateSetIterator::RECORD: " + record.getSAMString());
            //this.duplicateSet.add(this.wrappedIterator.next());
            this.duplicateSet.add(record);
        }
        
    }

    public void setScoringStrategy(final DuplicateScoringStrategy.ScoringStrategy scoringStrategy) {
        this.comparator.setScoringStrategy(scoringStrategy);
    }

    public DuplicateSet next() {
        DuplicateSet duplicateSet = null;
        
        int cmp = 0;

        while (0 == cmp) {
            if (!wrappedIterator.hasNext()) { // no more!
                //System.err.println("No more!");
                duplicateSet = this.duplicateSet;
                this.duplicateSet = new DuplicateSet();
                break;
            }
            else {
                // get another one
                final SAMRecord record = this.wrappedIterator.next();
                //System.err.print("DuplicateSet::next::RECORD: " + record.getSAMString());

                // assumes that the duplicate set always has at least one record inside!
                final SAMRecord representative = this.duplicateSet.getRepresentative();
                
                if (representative.getReadUnmappedFlag() || representative.isSecondaryOrSupplementary()) {
                    duplicateSet = this.duplicateSet;
                    this.duplicateSet = new DuplicateSet();
                    this.duplicateSet.add(record);
                    break; // exits the 0 == cmp loop
                }
                else {
                    // compare against the representative for set membership, not ordering
                    cmp = this.comparator.duplicateSetCompare(representative, record);
                    /*
                    System.err.println("*******");
                    System.err.println("DuplicateSet::next::cmp: " + cmp);
                    System.err.print("DuplicateSet::next::REPRE: " + representative.getSAMString());
                    System.err.print("DuplicateSet::next::OTHER: " + record.getSAMString());
                    */

                    if (0 < cmp) {
                        throw new SAMException("The input records were not sorted in duplicate order:\n" +
                        representative.getSAMString() + record.getSAMString());
                    } else if (0 == cmp) {
                        this.duplicateSet.add(record);
                        // loop around again!
                        //System.err.print("NEWRE: " + this.duplicateSet.getRepresentative().getSAMString());
                    } else {
                        duplicateSet = this.duplicateSet;
                        this.duplicateSet = new DuplicateSet();
                        this.duplicateSet.add(record);
                    }
                    //System.err.println("*******");
                }
            }
        }
        
        return duplicateSet;
    }

    public void close() { wrappedIterator.close(); }

    public boolean hasNext() { 
        return (!duplicateSet.isEmpty() || wrappedIterator.hasNext());
    }

    // Does nothing!
    public void remove() { }
}
