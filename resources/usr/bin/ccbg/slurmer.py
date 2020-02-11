import subprocess
import shlex
import time
import os
import re
import pprint as pp


job_ids    = []
job_status = []

RUNNING    = 1
COMPLETED  = 2
FAILED     = 4
UNKNOWN    = 8


def submit( cmd, limits="" ):
    

    SLURM_cmd  = " sbatch -J {} ".format( "CMD" )
    if ( limits is not None):
        SLURM_cmd += " {} ".format(limits)


    p = subprocess.Popen(shlex.split(SLURM_cmd), shell=False, 
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE)

    (stdout, stderr) = p.communicate("#!/bin/bash \n{}".format( cmd));
    status = p.wait()

    if status == 0:
#        print stdout
        job_id = re.match('Submitted batch job (\d+)', stdout)
        job_id = job_id.group( 1 )
        job_ids.append( job_id )
        return True
    else:
        print "Failed submitting {} job".formt(cmd )
        return False

def status( q_job_id ):

    p = subprocess.Popen(shlex.split("sacct --format JobIDRaw,ExitCode,State,Elapsed,CPUTime,MaxRss -np -j {}".format(q_job_id)), 
                         shell=False,
                         stdout=subprocess.PIPE)

    status = p.wait()
    (stdout, stderr) = p.communicate()


    if status is None:
        # This should not happen!
        job.status = manager.Job_status.UNKNOWN
    elif status == 0:

#        print stdout

        # If the job has finised it does a batch line as well, so
        # take this into account and return the second last line,
        # The last one is empty!

        if ( stdout.count("\n") > 1):
            lines = stdout.split("\n")
            stdout = lines[-2]

        fields = stdout.split("|")
        (job_id, exit_code, status, elapsed, cputime, max_mem, undef) = fields

        if ( status == 'FAILED'    or 
             status == "NODE_FAIL" or
             status == "CANCELLED" or
             status == "TIMEOUT"   or
             status == "PREEMPTED" ):
            return FAILED
        elif (status == 'RUNNING' or 
              status == 'PENDING' or
              status == 'SUSPENDED'):
            return RUNNING
        elif (status == 'COMPLETED' ):
            if ( os.path.exists("slurm-{}.out".format(q_job_id))):
                os.unlink("slurm-{}.out".format(q_job_id))
            return COMPLETED
                
        else:
            return UNKNOWN


def check_jobs():
    global job_status
    job_status = [-1]*len(job_ids )
    for i, job_id in enumerate( job_ids ):
        job_status[ i] = status( job_id )





def wait(sleep_time = 30):

    check_jobs()

    if len(job_status) == 0:
        return True

    while( job_status.count( RUNNING )):
        time.sleep( sleep_time )
        check_jobs()
        print("D:{} R:{} F:{} U:{}".format(job_status.count( COMPLETED),
                                           job_status.count( RUNNING),
                                           job_status.count( FAILED ),
                                           job_status.count( UNKNOWN )))


    if ( job_status.count( FAILED)):
        return False

    print "All jobs done..."
    return True


def reset():
    """ Resets the tracking of old commands
    Args:
      None

    Returns:
      None

    Raises:
      None

    """
    
    global job_ids
    job_ids    = []
    global job_status
    job_status = []
