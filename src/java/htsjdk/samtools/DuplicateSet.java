package htsjdk.samtools;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Stores a set of records that are duplicates of each other.  The first records in the list of records is 
 * considered the representative of the duplicate, and typically does not have it's duplicate flag set.  
 * The records' duplicate flag will be set appropriately as records are added.
 * 
 * At this time, this does not track optical duplicates.
 */
public class DuplicateSet {
    
    private final List<SAMRecord> records;
    
    private static final SAMRecordDuplicateComparator comparator = new SAMRecordDuplicateComparator();
    
    private boolean needsSorting = false;

    public DuplicateSet() {
        records = new ArrayList<SAMRecord>(10);
    }
    
    public void add(final SAMRecord record) {
        this.records.add(record);
        needsSorting = true;
    }
    
    private void sort() {
        if (!records.isEmpty()) {
            Collections.sort(records, this.comparator);
            
            final SAMRecord representative = records.get(0);
            
            // reset duplicate flags
            for (final SAMRecord record : records) {
                if (!record.getReadUnmappedFlag() && !record.isSecondaryOrSupplementary() && !record.getReadName().equals(representative.getReadName())) {
                    record.setDuplicateReadFlag(true);
                }
            }
            records.get(0).setDuplicateReadFlag(false);
        }
        needsSorting = false; // this could be in the if above if you think hard about it
    }
    
    public List<SAMRecord> getRecords() {
        if (needsSorting) {
            sort();
        }
        
        return this.records;
    }
    
    public SAMRecord getRepresentative() {
        if (needsSorting) {
            sort();
        }
        
        return records.get(0);
    }
    
    public int size() {
        return this.records.size();
    }
    
    public boolean isEmpty() {
        return this.records.isEmpty();
    }
}
