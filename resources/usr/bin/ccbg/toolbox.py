#!/usr/bin/python
#
# collection of misc useful python functions 
#
#
#
# Kim Brugger (21 Sep 2017), contact: kbr@brugger.dk


from __future__ import print_function, unicode_literals
import sys
import os
import re
import gzip




def open_file( filename ):
    """ universally open a file for reading

     if the filename contains a .gz, the file will be opened with gzip
     instead of the standard python open

    Args:
      filename (str): filename to read from

    Returns:
      filehandle

    """

    fh = None
    if ".gz" in filename:
        fh = gzip.open(filename, 'rb')
    else:
        fh = open(filename, 'rb')

    return fh




def print_stdout_or_file(line, fh=None):
    ''' wrapping function giving the option to either print to a file handle or stdout

    Args:
      line (str): string to print, no newline is added!
      fh (filehandle): if set, function will print to this instead

    Returns:
      None

    Raises:
        No exceptions are caught by the function
      
    '''


    if ( fh is not None):
        fh.write( line )
    else:
        sys.stdout.write(line)


    return None


def get_runfolder( name ):
    """Tries to identify the runfolder from a name/directory


    Args:
        name (str) : name to extract runfolder from

    Returns:
        runfolder (str)
        
    Raises:
        None

    ToDo: 
        should raise an error if providing a non-existing file
    """
    
    path = os.path.realpath( name )

    if ( not os.path.exists( name )):
#        print( "{} does not exist".format( name )) 
        return None

    # Remove the name, if it is a file
    if ( os.path.isfile(  path )):
        path = re.sub( r'(.*\/).*', r'\1', path)

    path = re.sub( 'vcfs', '', path)
    path = re.sub( 'stats', '', path)
    path = re.sub( 'bams', '', path)
    path = re.sub( 'logs', '', path)
    path = re.sub( 'tmp', '', path)
    path = re.sub( r'\/+', r'/', path)
    path = re.sub( r'\/+$', r'', path)
    dirs = path.split("/")
    return dirs[-1]

def get_sample_name( filename ):
    """extract the sample-name from a filename


    Args:
        name (str) : file to extract sample-name from

    Returns:
        sample-name (str)
        
    Raises:
        None

    ToDo: 
        should raise an error if providing a directory, or a
        non-existing file

    """

    if ( not os.path.isfile( filename ) ):
#        print( "{} is not a file".format( filename ))
        return None
    
    name = re.sub( r'.*\/(.*)', r'\1', filename)
    name = re.sub( r'(.*?)(\..*)', r'\1', name)

    return name



def get_sample_and_runfolder( filename ):
    """extract the sample and runfolder information from a filename


    Args:
        name (str) : file to extract sample-name from

    Returns:
        sample-name (str)
        runfolder-name (str)
        
    Raises:        
        None

    """
    sample_name = get_sample_name( filename )
    runfolder   = get_runfolder( filename )

    return sample_name, runfolder
