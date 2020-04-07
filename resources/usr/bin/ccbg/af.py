import sys

sys.path.append("/software/dev/lib//python2.7/site-packages/")
import pysam

def af_vcf(chrom, pos, ref, alt, filename):
    """ pulls info field information for a variant from a vcf

    If nothing is found returns None

    Args:
      chrom(str): chromsome name
      pos (int): position
      ref (str): reference allele
      alt (str): alt allele

    Returns:
      info( VariantRecordInfo ): the info field from the file, None if 

    Raises:
      Raises an IOError on file not found
    """

    handle = pysam.VariantFile( filename )

    for record in handle.fetch(chrom, pos - 1, pos + 1 ):
        if ( record.chrom == chrom and 
             record.pos   == pos   and
             record.ref   == ref   and
             alt in record.alts):
            return record.info['AF']


    return None


def latest_file( directory ):
    """ Find the files, sorted by name and version, in a directory

    Args:
      directory(str): directory to look in

    Returns:
      dict of lists of filenames and their versions

    """
    pass


if __name__ == "__main__":
    print( af_vcf("1", 906272, "A", "C", "/data/gemini/gemini_freq/gemini_180117.vcf.gz"))
    print( af_vcf("1", 906273, "A", "C", "/data/gemini/gemini_freq/gemini_180117.vcf.gz"))
    print( af_vcf("1", 139011,  "TTGAG", "T", "/data/projects/Population_databases/gnomad/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz"))
    print( af_vcf("1", 139011,  "TTGAG", "TT", "/data/projects/Population_databases/gnomad/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz"))
    print( af_vcf("1", 721694,  "A", "G", "/data/projects/matt/brca/data/ExAC.r0.3.1.sites.vep.vcf.gz"))
