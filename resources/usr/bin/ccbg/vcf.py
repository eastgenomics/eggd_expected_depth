
import sys
import re
import os 
import pprint as pp

def get_minimal_representation(pos, ref, alt): 
    """ produce the minimal representation of a variant

    Get the minimal representation of a variant, based on the ref + alt alleles in a VCF
    This is used to make sure that multiallelic variants in different datasets, 
    with different combinations of alternate alleles, can always be matched directly. 

    adopted from: http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/

    Args: 
        pos (int): genomic position in a chromosome
        ref (str): ref allele string
        alt (str): alt allele string

    Returns: 
        tuple: (pos, ref, alt) of remapped coordinate

    """

    # Ensure both the ref and alt bases are upper case.
    ref = ref.upper()
    alt = alt.upper()

    assert len(ref) >= 1 and len(alt) >= 1, "One or both alleles are empty"


    pos = int(pos)
    # If it's a simple SNV, don't remap anything
    if len(ref) == 1 and len(alt) == 1: 
        return pos, ref, alt
    else:
        # strip off identical suffixes
        while(alt[-1] == ref[-1] and min(len(alt),len(ref)) > 1):
            alt = alt[:-1]
            ref = ref[:-1]
        # strip off identical prefixes and increment position
        while(alt[0] == ref[0] and min(len(alt),len(ref)) > 1):
            alt = alt[1:]
            ref = ref[1:]
            pos += 1
        return pos, ref, alt 


def dephase_alleles( ref_allele, alleles ):
    """ dephase NMV alleles

    The alleles being dephased have to have the same length.

    Args:
      ref_mnv (str): the reference mnv
      mnvs (list of str): list of mnvs as strings, eg: ['ATC', 'CTG']

    Returns:
      ??

    """

    assert len(alleles) <=2, "Can only dephase a maximum of 2 alt alleles"


    # get rid of any alleles that is identical to the reference allele, and alleles with data in them
    non_ref_alleles = []
    for allele in alleles:
        assert len( allele ) >=1, "Cannot dephase an empty allele"
        if ref_allele != allele:
            non_ref_alleles.append( allele )

    alleles = non_ref_alleles

    # ensure all alleles are the same length
    allele_lengths =  [len(x) for x in alleles ]
    assert len( ref_allele ) == min( allele_lengths ) and len( ref_allele ) == max( allele_lengths ), "Alleles have different lengths"
        
    decomp_vars = {}
    # first gt is the reference
    for allele_index, allele in enumerate(alleles):
        
        for pos in range(0, len( ref_allele)):
            # the reference and the alt differs at this position
            if ( ref_allele[ pos ] != allele[ pos ]):

                # Extract the base
                base = allele[ pos ]
                # Never seen the pos and or base before, so add them to the dictionary
                if pos not in decomp_vars:
                    decomp_vars[ pos ] = {}

                if base not in decomp_vars[ pos ]:
                    decomp_vars[ pos ][base] = []

                # A little trick here. I like the single non-phased
                # non-ref variant genotype to be 0/1. Why? I dont
                # know, but by setting the allele nr to 1 in the case
                # of testing a single allle this makes things come out
                # as I want.
                if ( len( alleles ) == 1):
                    decomp_vars[ pos ][ base ].append( 1  )
                else:
                    decomp_vars[ pos ][ base ].append( allele_index  )

    res = []
    for pos in decomp_vars:

        # there is only one non-ref allele at this position
        if len(decomp_vars[ pos ]) == 1:

            for base in decomp_vars[ pos ] :
                # Both are non-ref alleles (count of two)
                if len ( decomp_vars[ pos ][ base ] ) == 2:
                    res.append({ 'offset': pos,"ref": ref_allele[pos], 'alts': [base], 'GT':(1,1)})
                # Allele 1 is non ref allele
                elif decomp_vars[ pos ][ base ][ 0 ] == 0:
                    res.append({ 'offset': pos,"ref": ref_allele[pos], 'alts': [base], 'GT':(1,0)})
                # Allele 2 is non ref allele
                else:
                    res.append({ 'offset': pos,"ref": ref_allele[pos], 'alts': [base], 'GT':(0,1)})

        else:

            # there are two non-ref bases at this position, so this is
            # the best way I can figure out how to do this for now
            bases = {}

            for base in decomp_vars[ pos ]:
                bases[ decomp_vars[ pos ][ base ][0] ] = base

            res.append({ 'offset': pos, "ref": ref_allele[pos], 'alts': [bases[0], bases[1]], 'GT':(1,2)})

    return res



def smith_waterman_alignment( seq1, seq2 ):


    # These scores are taken from Wikipedia.
    # en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
    match    = 2
    mismatch = -1
    gap      = -1


    score_matrix = [[0 for col in range(len(seq1) + 1 )] for row in range(len(seq2)+1)]
    

    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            score = calc_score(score_matrix, i, j)
            if score > max_score:
                max_score = score
                max_pos   = (i, j)

            score_matrix[i][j] = score


    pp.pprint( score_matrix )
    return;
