import glob
import sys
import pysam

def check_read(r):
    if r.is_secondary: return False
    if r.is_unmapped: return False
    if r.is_supplementary: return False
    return True


file_pairs = [12066571, 11557403, 11876270, 12888782, 11275417, 7374023, 12136999, 3018723, 8376929]
NPAIRS = sum(file_pairs)

if __name__ == "__main__":

    workdir = sys.argv[1]
    counts = []
    counts_unique = []
    n = 0
    n2 = 0
    
    for filename, pairs in zip(sorted(glob.glob(workdir + "/*_sorted.bam")), file_pairs):
        print filename
        with pysam.AlignmentFile(filename, 'rb') as bamfile:
            filereads = 0
            cond_counts = []
            n += sum([x[1] for x in bamfile.get_index_statistics()]) 
            print "mapped + unmapped = ", bamfile.mapped, bamfile.unmapped, bamfile.mapped+ bamfile.unmapped
            
            reads = {}
            reads_u = {}
            for read in bamfile.fetch():
                if not (read.is_unmapped):
                    name = read.query_name
                    if read.is_read1: name += "/1"
                    else: name += "/2"
                    if name in reads:
                        reads[name] += 1
                    else:
                        reads[name] = 1
                    if read.mapping_quality > 3: #uniquely mapped:
                        reads_u[name] = 1
                    
            
            counts.append(len(reads))
            #counts_unique.append(reads.values().count(1))
            counts_unique.append(reads_u.values().count(1))
            #print bamfile.count(read_callback=check_read), len(reads), reads.values().count(1), max(reads.values())
            print "mapped:", len(reads), "mapped unique:", len(reads_u), "max times mapped:", max(reads.values())
            
    mean = sum([c/(x*2.) for c, x in zip(counts, file_pairs)])/len(counts)
    mean2 = sum([c/(x*2.) for c, x in zip(counts_unique, file_pairs)])/len(counts)
    
    print "mean file coverage (primary)", mean
    print "mean file coverage (unique)", mean2
    
