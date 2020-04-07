#
#
#
#
#
# Kim Brugger (21 Sep 2017), contact: kbr@brugger.dk


from __future__ import print_function, unicode_literals
import inspect

import pprint as pp
import sys
import os 
import re
import subprocess
import tempfile


sys.path.append('/software/dev/lib/python2.7/site-packages/')
import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.ext.automap import automap_base

from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound

import pysam



def init_database( dbase, db_host, db_user, db_passwd, db_model='mysql') :

    '''
    initiate the connection to the database.

    The user can override the connection details if using the non-production database

    dbase = database name
    db_host = the server name of the database 
    db_user = who to connect as
    db_passwd = the password for the user, if applicable 
    

    '''
    Session, Base= None, None

    try:
        database_string = "{model}://{user}:{passwd}@{host}/{db}".format(model=db_model, user=db_user, passwd=db_passwd, host=db_host, db=dbase)
        engine = create_engine(
            database_string,
            isolation_level="READ UNCOMMITTED"
            )
        
        Base = automap_base()
        Base.prepare(engine, reflect=True)
        
        session = sessionmaker(bind=engine)
        Session = session()
    except Exception as e: 
        print( e )
        raise RuntimeWarning("Could not connect to database %s" % database_string)


    return Session, Base

def database_connect( db_host, db_user, db_passwd, db_model='mysql') :

    '''
    initiate the connection to the database-server. This is for
    lowlevel sql functions only such as creating new databases on the
    fly etc

    The user can override the connection details if using the non-production database

    db_host = the server name of the database 
    db_user = who to connect as
    db_passwd = the password for the user, if applicable 
    

    '''
    Session, Base= None, None

    try:
        database_string = "{model}://{user}:{passwd}@{host}".format(model=db_model, user=db_user, passwd=db_passwd, host=db_host)
        engine = create_engine(
            database_string,
            isolation_level="READ UNCOMMITTED"
            )

        connection = engine.connect()
    except Exception as e: 
        print( e )
        raise RuntimeWarning("Could not connect to database %s" % database_string)


    return connection


def to_dict( result ):

    if  isinstance(result, sqlalchemy.orm.query.Query):
        result = results
        res = []
        for result in results.all():
            r = result.__dict__
            if '_sa_instance_state' in r:
                del (r['_sa_instance_state'])
            res.append( r )

        if ( len( res ) == 1):
            return res[0]
            
        return res

    else:
        results = {}
        for key in result.__dict__:
            if key.startswith("_"):
                continue
            results[ key ] = result.__dict__[ key ]
            
        return results




def get_or_create(session, 
                  model,
                  create_method='',
                  create_method_kwargs=None,
                  **kwargs):
    '''
    Checks if a database entry with the wanted data already
    exists, if it does return it, otherwise create entry and return
    that.
    
    Found on the internet at: http://skien.cc/blog/2014/01/15/sqlalchemy-and-race-conditions-implementing/
    
    tweked slightly
    '''

    try:
        return session.query(model).filter_by(**kwargs).one(), False
    except NoResultFound:
        kwargs.update(create_method_kwargs or {})
        created = getattr(model, create_method, model)(**kwargs)
        try:
            session.add(created)
            return created, True
        except IntegrityError:
            session.rollback()
            return Session.query(model).filter_by(**kwargs).one(), False




def sry_present( binned_file ):

    regions = regions_for_gene( 'SRY', ref='NM_003140.2' )
    
    overall_coverage, coverage = depths_for_regions( binned_file, regions)


    if ( overall_coverage <  100):
        return False

    return True


def write_regions_to_file( regions, filename ):

    fh = open( filename, "w")
    for region in regions:
        os.write(fh,  "\t".join(map(str, region))+"\n")
    os.close( fh )


def write_regions_to_tempfile( regions ):

    fh, filename = tempfile.mkstemp(dir="/tmp/")
    for region in regions:
        os.write( fh,  "\t".join(map(str, region))+"\n")
    os.close( fh )

    return filename





