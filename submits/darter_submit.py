import sys, os, math
import shlex, subprocess, shutil

def run_subproc(args):
    """ Run subprocess given args as a command line string"""
    print(args)
    args = shlex.split(args)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pcomm = p.communicate()
    poutput = pcomm[0] #p.stdout.read()
    perr = pcomm[1] #p.stderr.read()
    print(poutput)
    print(perr)
    preturncode = p.wait()
    if(preturncode != 0):
        raise Exception("{0} failed!".format(args[0]))
    return pcomm

def submit_job(paramfile, prev_jobid, s, es):
    """ Pass None to prev_jobid if using param_file.in as paramfile """

    # Call sed to change starting step number to s
    args = "sed -i '5s/.*/{0} {1}\t\t# starting step, ending step/' {2}".format(s, es, paramfile)
    pcomm = run_subproc(args)

    # Do the submitting
    if(prev_jobid == None):
        args = "sed -i '11s/.*/aprun \-n $PBS_NNODES .\/hrmc $PBS_JOBID {0}/' {1}".format(paramfile.replace('/', '\/'), 'submits/pbs.sh')
        pcomm = run_subproc(args)
        #args = "qsub submits/pbs.sh {0}".format(paramfile)
        args = "qsub submits/pbs.sh"
    else:
        # Call sed to change modelfile based on previous job id
        args = "sed -i '2s/.*/model_final_{0}.ocoee.nics.utk.edu.xyz/' {1}".format(prev_jobid,paramfile)
        pcomm = run_subproc(args)
        args = "sed -i '11s/.*/aprun \-n $PBS_NNODES .\/hrmc $PBS_JOBID {0}/' {1}".format(paramfile.replace('/', '\/'), 'submits/pbs.sh')
        pcomm = run_subproc(args)
        #args = "qsub --depend=afterok:{0} submits/pbs.sh {1}".format(prev_jobid,paramfile)
        args = "qsub -W depend=afterok:{0} submits/pbs.sh".format(prev_jobid)
    pcomm = run_subproc(args)
    jobid = int( pcomm[0].strip().split('\n').pop().split('.')[0]) #'595863.ocoee.nics.utk.edu'
    #jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])
    return jobid

def main():
    """ If you want to resume a job, all you have to do is change ss """
    ts = 4000000 # total number of steps in the simulation
    js = 1500000 # number of steps per job
    ss = 0 # starting step, this number gets changed in the loop ~ incremented by js on each submit
    starting_param_filename = 'param_file.in' # It is expected that the starting parameters in here are correct

    for s in range(ss, ts, js):
        if(s > ss):
            # Copy the starting parameter file to another filename which includes the step being started
            src = os.getcwd() + '/' + starting_param_filename
            dst = os.getcwd() + '/param_file.{0}.in'.format(s)
            shutil.copyfile(src, dst)
            # Submit the job using that parameter file, but change it first
            jobid = submit_job(dst, jobid, s, min(s+js, ts))
        else:
            # On the first run, just submit the starting parameter file
            dst = starting_param_filename
            jobid = submit_job(dst, None, s, min(s+js, ts))
        print("Jobid: {0}".format(jobid))


if __name__ == '__main__':
    main()
