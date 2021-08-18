# PBS TORQUE plug-in

import os
import numpy as np
from jobdefs import Job
from identity import *

SUB = "qsub"

_BATCHSCRIPT = ("#!/bin/bash -l                                                  \n"+
              "#PBS -l nodes=%d:ppn=%d                                         \n"+
              "#PBS -q %s                                                      \n"+
              "#PBS -m %s                                                      \n"+
              "#PBS -r n                                                        \n"+
              "#PBS -l walltime=48:00:00                                        \n"+
              "#PBS -N %s                                                       \n"
              "# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE"+
              " nodes,ppn,walltime and my_job_name VALUES                       \n"+
              "cd $PBS_O_WORKDIR                                                \n"+
              "module load "+GCCMOD              +"                            \n"+
              "module load "+PYTHONMOD           +"                            \n"+
              "module load "+INTELMOD            +"                              \n"+
              "module load "+MPIMOD              +"                              \n")

def BATCHSCRIPT(job,notify):
    return _BATCHSCRIPT%(1,job.ncores,job.queue,notify,job.name)
    

MODELS = {"plasim":1,                #tasks per node (1 workq node on Sunnyvale has 8 threads)
          "sbdart":8,           #Here we use 'task' to mean a Sunnyvale job, as opposed to the
          "sbdart_earth":8,     #HPC convention of a task being a thread or process. This way our
          "sbdart_locked":8,    #code is MPI/OpenMP-agnostic.
          "postprocess":8,      
          "postprocess_earth":8,
          "postprocess_locked":8,
          "lmdz":8,
          "mitgcm":6}             

def getjobs():
    print("Checking jobs")
    os.system("qstat -u "+USER+" > cjobs.tmp")
    cjf = open("cjobs.tmp","r")
    joblist = cjf.read().split('\n')[5:-1]
    cjf.close()
    os.system("rm cjobs.tmp")
    resources={}
    for m in list(MODELS.keys()):
        resources[m] = np.zeros(256)
    tags = []
    #This part may need changing depending on how job tags are handled.
    for j in joblist:
        job = j.split()
        #if job[3][5:]!="lmdz-":
            #tags.append(job[0])
        tags.append(job[0])
    for t in tags:
        print("Looking up "+t)
        os.system("qstat -f "+t+" > jinfo.tmp")
        jf = open("jinfo.tmp","r")
        jinfo = jf.read().split('\n')[1:-2]
        while '' in jinfo:
            jinfo.remove('')
        jf.close()
        os.system("rm jinfo.tmp")
        ncpus = 1
        for l in jinfo:
            if len(l.split())>0:
                if l.split()[0]=="init_work_dir":
                    workdir = l.split()[2]
                #if l.split()[0]=="Resource_List.ncpus":
                    #ncpus = int(l.split()[2])
                if l.split()[0]=="Resource_List.nodes":
                    ncpus = int(l.split()[2].split("=")[1])
        ourjob=True
        try:
            job = np.load(workdir+"/job.npy").item()
        except:
            for nl in range(0,len(jinfo)):
                l = jinfo[nl]
                if len(l.split())>0:
                    if l.split()[0]=="init_work_dir":
                        try:
                            workdir = l.split()[2] + jinfo[nl+1].split()[0]
                        except:
                            workdir = l.split()[2]
            try:
                job = np.load(workdir+"/job.npy").item()
            except:
                ourjob=False
        if ourjob:
            jid = job.home
            if jid>=len(resources[job.model]):
                tmp = np.zeros(jid+100)
                tmp[:len(resources[job.model])] = resources[job.model][:]
                resources[job.model] = tmp
            resources[job.model][jid] = float(ncpus)/8.0#MODELS[job.model]
    
    return resources
