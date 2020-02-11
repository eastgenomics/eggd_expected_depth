
import sys
import re
import os 
import pprint as pp

import pysam

import ccbg.toolbox as toolbox



# ============== Generate depth functions =========================

INTERVALS = [20, 10, 5, 1, 0]
binned_depth = False

def make_region( bamfile, chrom=None, start=None, end=None, fh=None):
    """Finds the depth for a region in a bamfile defined by chrom, start and end
    
    
    Args:
        bamfile: pysam bamfile handle
        chrom (str): chromosome
        start (int): start position (included)
        end (int): end position (included)
        fh(file-handle): filehandle to write to, otherwise stdout

    Returns:
        print depths to stdout, but returns nothing

    
    Raises:
        No exceptions are caught by the function
    """


    start = int( start )
    end   = int( end   ) 

    first_entry = True

    prev_end = None

    for rec in bamfile.pileup(chrom, start, end):
        # if the first position is after the wanted start position, pad with 0's
        if first_entry and rec.pos + 1 > start:
            for i in range( start, rec.pos + 1):
                toolbox.print_stdout_or_file( "\t".join(map(str, [bamfile.getrname(rec.tid), i, 0]))+"\n", fh)
        
        first_entry = False

        # old code, that might be useful later
        # If we are doing binned depths, iterate through the values to we find the correc bin.
        if ( binned_depth ):
            for bin_depth in sorted( bins, reverse=True ):
                if ( rec.n >= bin_depth ):
                    rec.n = bin_depth
                    break



        # check positions as we get reads spanning the region and thus would get depths outside the wanted region.
        if rec.pos + 1 < start:
            continue
        if rec.pos + 1 > end:
            return

        # there is an area in the wanted region that does not have any reads, fill in with 0's
        if ( prev_end is not None and prev_end < rec.pos + 1 ):
            for i in range( prev_end+1, rec.pos + 1):
                toolbox.print_stdout_or_file( "\t".join(map(str, [bamfile.getrname(rec.tid), i, 0]))+"\n", fh)

        toolbox.print_stdout_or_file( "\t".join(map(str, [bamfile.getrname(rec.tid), rec.pos + 1, rec.n]))+"\n", fh)
        prev_end = rec.pos + 1

    # No reads in the region, so make a lot of 0's
    if ( prev_end is None):
        for i in range( start, end +1):
            toolbox.print_stdout_or_file( "\t".join(map(str, [chrom, i, 0]))+"\n", fh)

    # if the past position is before the wanted endpos, pad with 0's
    elif ( prev_end < end ):
        for i in range( prev_end+1, end +1):
            toolbox.print_stdout_or_file( "\t".join(map(str, [bamfile.getrname(rec.tid), i, 0]))+"\n", fh)



def make_region_in_blocks(bamfile, chrom=None, start=None, end=None, fh=None):
    """Print the depth, in blocks, for a region in a bamfile defined by chrom, start and end
    
    Args:
        bamfile: pysam bamfile handle
        chrom (str): chromosome
        start (int): start position (included)
        end (int): end position (included)
        fh(file-handle): filehandle to write to, otherwise stdout
    
    Returns:
        print depths to stdout, but returns nothing

    
    """


    start = int( start )
    end   = int( end   ) 

    ( prev_chrom, prev_start, prev_ref, prev_end, prev_depth) = ( None, None, None, None, None)
    prev_printed = None
    for rec in bamfile.pileup(chrom, start, end):
        # if the first position is after the wanted start position, pad with 0's
        if prev_start is  None and rec.pos + 1 > start:
            toolbox.print_stdout_or_file( "\t".join(map(str, [bamfile.getrname(rec.tid), start, rec.pos, 0]))+"\n", fh)

        # check positions as we get reads spanning the region and thus would get depths outside the wanted region.
        if rec.pos + 1 < start:
            continue
        if rec.pos + 1 > end:
            break

        # old code, that might be useful later
        # If we are doing binned depths, iterate through the values to we find the correc bin.
        if ( binned_depth ):
            for bin_depth in sorted( bins, reverse=True ):
                if ( rec.n >= bin_depth ):
                    rec.n = bin_depth
                    break


        # there is an area in the wanted region that does not have any reads.
        # Print prev region and print a new 0 reads block
        if ( prev_end is not None and prev_end + 1 < rec.pos ):
            toolbox.print_stdout_or_file( "\t".join(map(str, [prev_chrom, prev_start, prev_end, prev_depth]))+"\n", fh)
            toolbox.print_stdout_or_file( "\t".join(map(str, [prev_chrom, prev_end+1, rec.pos , 0]))+"\n", fh)
            ( prev_chrom, prev_start, prev_end, prev_depth) = ( bamfile.getrname(rec.tid), rec.pos+1, rec.pos+1, rec.n)
        
        # first data set the counters
        if ( prev_depth is None):
            ( prev_chrom, prev_start, prev_end, prev_depth) = ( bamfile.getrname(rec.tid), rec.pos+1, rec.pos+1, rec.n)
        # Same depth as the prev base increase prev_endpos with one
        elif (  prev_depth == rec.n ):
            prev_end = rec.pos + 1
        else:
            toolbox.print_stdout_or_file( "\t".join(map(str, [prev_chrom, prev_start, prev_end, prev_depth]))+"\n", fh)
            prev_printed =  rec.pos 
            ( prev_chrom, prev_start, prev_end, prev_depth) = ( bamfile.getrname(rec.tid), rec.pos+1, rec.pos+1, rec.n)

    # print the last block seen
    if ( prev_chrom is not None):
        toolbox.print_stdout_or_file( "\t".join(map(str, [prev_chrom, prev_start, prev_end, prev_depth]))+"\n", fh)

    if ( prev_end is None):
        toolbox.print_stdout_or_file( "\t".join(map(str, [chrom, start, end, 0]))+"\n", fh)


    # if the past position is before the wanted endpos, pad with a 0 block
    if ( prev_end is not None and prev_end < end ):
        toolbox.print_stdout_or_file( "\t".join(map(str, [prev_chrom, prev_end+1, end, 0]))+"\n", fh)

def make_regions_from_bedfile( bamfile, bedfile, block = False, fh=None):
    """Calculates depths based on entries in a bedfile
    
    
    Args:
        bamfile: pysam bamfile handle
        bedfile(str): name of the bedfile to read regions from
        block(bool): report as blocks or single base resolution, default false
        fh(file-handle): filehandle to write to, otherwise stdout

    Returns:
        print depths to stdout, but returns nothing

    
    Raises:
        No exceptions are caught by the function
    """

    bed_fh = open( bedfile, 'r')
    for line in bed_fh.readlines():
        line = line.strip("\n")

        fields = line.split("\t")

        chrom, start, end =  fields[:3]
        start = int( start)
        end   = int( end  )

        if block and block is not None:
            make_region_in_blocks(bamfile, chrom, start + 1, end , fh)
        else:
            make_region(bamfile, chrom, start + 1, end, fh)




def compress_depth_file( filename, outfile=None, delete_file=True):
    """compresses a tab file.

    Args:
        filename (str): file to compress
        outfile ( str): name of compressed file, default is [filename].gz
        delete_file (bool): delete original file after compression, default is True

    Returns:
        filename (str): filename of compressed file

    
    """

    if ( outfile is None):
        outfile = "{}.gz".format( filename )

    pysam.tabix_compress( filename, outfile , force=True )

    if ( delete_file ):
        os.unlink( filename )

    return outfile



def index_depth_file( filename, single_base=False ):
    """indexes a compressed depth tab file
    
    prev existing index files will be overwritten

    Args:
        filename (str): file to index
        single_base (bool): the file contains single base resolution data

    Returns:
        None

    
    """

    if ( single_base ):
        pysam.tabix_index(filename, force=True , seq_col=0,start_col=1, end_col=1)
    else:
        pysam.tabix_index(filename, force=True , seq_col=0,start_col=1, end_col=2)

    return None





# =============== reporting functions ====================


def view_region( depth_file, chrom=None, start=None, end=None):
    """reports the depths for, optionally a region, in a depth-file. Region is defined by chrom, start and end
    
    
    Args:
        depth_file: name of file to read from
        chrom (str): chromosome
        start (int): start position (included)
        end (int): end position (included)

    Returns:
        print depths to stdout, but returns nothing

    
    Raises:
        No exceptions are caught by the function
    """

    depth_in = pysam.TabixFile( depth_file )
   

    if ( chrom is None):
        iterator = depth_in.fetch()
    else:
        try:
            iterator = depth_in.fetch(chrom, int(start)-1, int(end)-1)
        except:
            return

    for row in iterator:
        print( row )


def view_region_from_bedfile( depth_file, bed_file):
    """reports the depths for, optionally a region, in a depth-file. Region is defined by chrom, start and end
    
    
    Args:
        depth_file: name of file to read from
        bedfile(str): name of the bedfile to read regions from

    Returns:
        print depths to stdout, but returns nothing

    
    Raises:
        No exceptions are caught by the function
    """

    bed_fh = open( bed_file, 'r')
    for line in bed_fh.readlines():
        line = line.strip("\n")

        fields = line.split("\t")

        chrom, start, end =  fields[:3]

        start = int( start )
        end   = int( end   )

        view_region( depth_file, chrom, start + 1, end)

def view_regions( depth_file, regions):
    """reports the depths for, regions
    
    
    Args:
        depth_file: name of file to read from
        regions(list): list of list of regions consiting of chrom, start, end

    Returns:
        print depths to stdout, but returns nothing

    
    Raises:
        No exceptions are caught by the function
    """

    for region in regions:
        chrom, start, end = region

        view_region( depth_file, chrom, int(start), int(end))




def coverage_region( depth_file, chrom=None, start=None, end=None, min_coverage=1):
    """reports a full depths/coverage report for a depth-file, limited by region defined by chrom, start, and end
        
    Args:
        depth_file: name of file to read from
        chrom (str): chromosome
        start (int): start position (included)
        end (int): end position (included)
        min_coverage(int): depth cutoff for percent coverage, default 1

    Returns:
        dict of values. 

    
    Raises:
        No exceptions are caught by the function
    """

    depth_in = pysam.TabixFile( depth_file )

    start        = int( start )
    end          = int( end )
    min_coverage = int( min_coverage )

    coverage = {}

    coverage[ 'bases_above_min' ] = 0
    coverage[ 'min_depth' ]       = None
    coverage[ 'mean_depth' ]      = 0
    coverage[ 'max_depth' ]       = None
    coverage[ 'summed_depth' ]    = 0
    coverage[ 'percent']          = 0
    coverage[ 'length' ]          = end-start +1
    
    coverage[ '20x' ]         = 0
    coverage[ '10-19x' ]      = 0
    coverage[ '6-9x' ]        = 0
    coverage[ '1-5x' ]        = 0
    coverage[ '0x' ]          = 0
    coverage[ 'N/A' ]         = 0

    coverage['blocks'] = {}
    coverage['blocks'][ '10-19x' ]      = []
    coverage['blocks'][ '6-9x' ]        = []
    coverage['blocks'][ '1-5x' ]        = []
    coverage['blocks'][ '0x' ]          = []
    coverage['blocks'][ 'N/A' ]         = []

    first_row = True

    try:
        iterator = depth_in.fetch(chrom, int(start) - 1, int(end) )
    except:
        return None


    for row in iterator:

#        print( row )
        
        # read in line, split and cast to int where appropriate
        block_chrom, block_start, block_end, block_depth = row.split("\t")
        block_start = int( block_start )
        block_end   = int( block_end   )
        block_depth = int( block_depth )

        # Asked for a defined region, and the returned block is outside the requested area, trim them back
        if ( chrom is not None):
            if ( block_start < start ):
                block_start = start

            if ( block_end > end ):
                block_end = end

            # The region asked for is not available
            if (first_row and start < block_start ):
#                print( "{} < {}".format( start, d_start ))
                coverage[ 'N/A' ] += block_end - block_start +1
                coverage['blocks'][ 'N/A' ].append( "{}:{}-{}".format(block_chrom, start, block_start - 1))

            first_row = False
                

            if ( coverage[ 'max_depth' ] is None or 
                 block_depth > coverage[ 'max_depth' ]):
                coverage[ 'max_depth' ] =  block_depth 
                      
            if ( coverage[ 'min_depth' ] is None or 
                 block_depth < coverage[ 'min_depth' ]):
                coverage[ 'min_depth' ] = block_depth 

            if ( block_depth >= min_coverage):
                coverage[ 'bases_above_min' ] += block_end  -  block_start  + 1
#                bases_above_min += int( block_end ) - int( block_start ) + 1

            coverage[ 'summed_depth' ] +=  block_depth *(block_end - block_start + 1)

            if ( block_depth >= 20):
                coverage[ '20x' ] += block_end  -  block_start  + 1
            elif( block_depth >= 10):
                coverage[ '10-19x' ] += block_end -  block_start +1
                coverage['blocks'][ '10-19x' ].append("{}:{}-{}".format(chrom, block_start, block_end))
            elif( block_depth >= 6):
                coverage[ '6-9x' ] += block_end -  block_start +1
                coverage['blocks'][ '6-9x' ].append("{}:{}-{}".format(chrom, block_start, block_end))
            elif( block_depth >= 1):
                coverage[ '1-5x' ] += block_end -  block_start +1
                coverage['blocks'][ '1-5x' ].append("{}:{}-{}".format(chrom, block_start, block_end))
            elif (block_depth == 0 ):
#                print ("{}:{}-{}".format(chrom, d_start, d_end))
                coverage[ '0x' ] += block_end -  block_start +1
                coverage['blocks'][ '0x' ].append("{}:{}-{}".format(chrom, block_start, block_end))
            else:
                print( "Unknown depth: {}".format( depth ))
                exit()

    coverage['blocks'][ '10-19x' ]      = _merge_regions( coverage['blocks'][ '10-19x' ])
    coverage['blocks'][ '6-9x' ]        = _merge_regions( coverage['blocks'][ '6-9x' ]  )
    coverage['blocks'][ '1-5x' ]        = _merge_regions( coverage['blocks'][ '1-5x' ]  )
    coverage['blocks'][ '0x' ]          = _merge_regions( coverage['blocks'][ '0x' ]    )
    coverage['blocks'][ 'N/A' ]         = _merge_regions( coverage['blocks'][ 'N/A' ]   )

    if ( block_start is None and block_end is None):
        coverage[ 'N/A' ] += end - start + 1
        coverage['blocks'][ 'N/A' ].append( "{}:{}-{}".format(chrom, start, end))

    elif ( block_end is not None and end > block_end ):
        coverage[ 'N/A' ] += end - block_end
        coverage['blocks'][ 'N/A' ].append( "{}:{}-{}".format(block_chrom, block_end + 1, end))

    coverage['mean_depth'] = 1.0*coverage['summed_depth']/coverage['length']
    coverage['percent']    = 1.0*coverage['bases_above_min']/coverage['length']


#    pp.pprint( coverage )

    return coverage

def coverage_regions_from_bedfile( depth_file, bed_file, min_coverage=20):
    """reports a full depths/coverage report for a depth-file, limited by regions in the bedfile
    
    
    Args:
        depth_file: name of file to read from
        chrom (str): chromosome
        start (int): start position (included)
        end (int): end position (included)
        min_coverage(int): depth cutoff for coverage, default 20

    Returns:
        dict of values. 

    
    Raises:
        No exceptions are caught by the function
    """

    coverages = {}

    total_length          = 0
    total_bases_above_min = 0
    total_summed_depth    = 0
    total_min_depth       = None
    total_max_depth       = None


    bed_fh = open( bed_file, 'r')
    for line in bed_fh.readlines():
        line = line.strip("\n")

        fields = line.split("\t")

        chrom, start, end =  fields[:3]
        start = int( start )
        end   = int( end   )
        id = "_".join(fields[3:])
        if id is None or id == "":
            id = "{}:{}-{}".format( chrom, start + 1, end)

        coverage = coverage_region( depth_file, chrom, start + 1, end, min_coverage)
        
        total_length += end - start # as 0 indexed and non inclusive end, we dont have to add one this time!
        if ( coverage is not None ):
            total_bases_above_min += coverage[ 'bases_above_min' ]
            total_summed_depth += coverage[ 'summed_depth' ]

        
            if ( total_min_depth is None or total_min_depth > coverage['min_depth']):
                total_min_depth = coverage['min_depth']

                if ( total_max_depth is None or total_max_depth < coverage['max_depth']):
                    total_max_depth = coverage['max_depth']
        


        coverages[ id ] = coverage

#    return { 'percent': 1.0*bases_above_min/total_length, 'length': total_length, 'bases_above_min': bases_above_min }

    coverages[ 'stats'] = {}
    coverages[ 'stats'][ 'length']          = total_length
    coverages[ 'stats'][ 'bases_above_min'] = total_bases_above_min
    coverages[ 'stats'][ 'min_depth']       = total_min_depth
    coverages[ 'stats'][ 'mean_depth']      = 1.0*total_summed_depth/total_length
    coverages[ 'stats'][ 'max_depth']       = total_max_depth


    return coverages



def coverage_regions( depth_file, regions, min_coverage=20):
    """reports a full depths/coverage report for a depth-file, limited by regions
    
    
    Args:
        depth_file: name of file to read from
        regions(list): list of list of regions consiting of chrom, start, end
        min_coverage(int): depth cutoff for coverage, default 20

    Returns:
        dict of values. 

    
    Raises:
        No exceptions are caught by the function
    """

    coverages = {}

    total_length          = 0
    total_bases_above_min = 0
    total_summed_depth    = 0
    total_min_depth       = None
    total_max_depth       = None


    for region in regions:
        chrom, start, end =  region[:3]
        start = int( start )
        end   = int( end )
        id = "_".join(region[3:])
        if id is None or id == "":
            id = "{}:{}-{}".format( chrom, start, end)

        coverage = coverage_region( depth_file, chrom, start, end, min_coverage)
        
        total_length += end - start + 1
        if coverage is not None:
            total_bases_above_min += coverage[ 'bases_above_min' ]
            total_summed_depth += coverage[ 'summed_depth' ]
        
            if ( total_min_depth is None or total_min_depth > coverage['min_depth']):
                total_min_depth = coverage['min_depth']
                
                if ( total_max_depth is None or total_max_depth < coverage['max_depth']):
                    total_max_depth = coverage['max_depth']
                    
        coverages[ id ] = coverage

#    return { 'percent': 1.0*bases_above_min/total_length, 'length': total_length, 'bases_above_min': bases_above_min }

    coverages[ 'stats' ] = {}
    coverages[ 'stats' ][ 'length']          = total_length
    coverages[ 'stats' ][ 'bases_above_min'] = total_bases_above_min
    coverages[ 'stats' ][ 'min_depth']       = total_min_depth
    coverages[ 'stats' ][ 'mean_depth']      = 1.0*total_summed_depth/total_length
    coverages[ 'stats' ][ 'max_depth']       = total_max_depth


    return coverages



def _merge_regions( regions ):
    """ merges neighbouring regions 

    Args:
      regions(list): list of region strings

    Returns:
        list of merged region strings

    Raises:
        No exceptions are caught by the function
    """

    merged_list = []

    prev_chrom, prev_start, prev_end = None, None, None
    for region in regions:
        region_chrom, region_start, region_end = re.split(r'[:-]', region)
        region_start = int(region_start)
        region_end   = int(region_end  )

        # This is the first entry, so set the prev-vars for tracking
        if ( prev_chrom is None):
            prev_chrom, prev_start, prev_end = region_chrom, region_start, region_end

        # This is part of the block we are currently growing, important check is the chromosome
        elif ( prev_end + 1  == region_start and prev_chrom == region_chrom ):
            prev_end = region_end

        else:
            # Append block, and start the next one
            merged_list.append("{}:{}-{}".format( prev_chrom, prev_start, prev_end))
            prev_chrom, prev_start, prev_end = region_chrom, region_start, region_end

    
    # append the last entry, if prev_chrom is not none. This only happens if the input is an empty list
    if ( prev_chrom is not None ):
        merged_list.append("{}:{}-{}".format( prev_chrom, prev_start, prev_end))
                               

    return merged_list
