# SLURM plug-in

import os
import numpy as np
from crawldefs import Job
from identity import *

SUB = "sbatch"

USER = "t-98b023"

def BATCHSCRIPT(job,notify):
    return _BATCHSCRIPT%(job.name,job.name,job.name,job.ncores,16,job.queue,
                         notify,job.top+"plasim/job"+str(job.home))

_BATCHSCRIPT = ("#!/bin/bash                                                  \n"+
                "#SBATCH --job-name=%s                                        \n"+
                "#SBATCH --output=%j_%s.out                                   \n"+
                "#SBATCH --error=%j_%s.err                                    \n"+
                "#SBATCH --ntasks=%d                                          \n"+
                "##SBATCH --ntasks-per-node=%d                                \n"+
                "#SBATCH --mem-per-cpu=2000M                                  \n"+
                "#SBATCH --account=%s                                         \n"%ACCOUNT+
                "#SBATCH --partition=%s                                       \n"+
                "#SBATCH --time=36:00:00                                      \n"+
                "#SBATCH --mail-type=%s                                       \n"+
                "#SBATCH --mail-user=%s                                       \n"%EMAIL+
                "#make sure newFitModel.m code is in the following location   \n"+
                "cd %s                                                        \n")

#job.queue could for example be 'broadw1'


MODELS = {"rossby-block":8,     #tasks per node (1 workq node on Sunnyvale has 8 threads)
          "plasim":1}           #Here we use 'task' to mean a Sunnyvale job, as opposed to the
                                #HPC convention of a task being a thread or process. This way our
                                #code is MPI/OpenMP-agnostic.

def getjobs():
    print("Checking jobs")
    os.system("squeue -u "+USER+" > cjobs.tmp")
    cjf = open("cjobs.tmp","r")
    joblist = cjf.read().split('\n')[1:-1]
    cjf.close()
    os.system("rm cjobs.tmp")
    resources={}
    for m in list(MODELS.keys()):
        resources[m] = np.zeros(256)
    tags = []
    for j in joblist:
        job = j.split()
        tags.append(job[0])
    for t in tags:
        print("Looking up job "+t)
        os.system("scontrol show job "+t+" > jinfo.tmp")
        jf = open("jinfo.tmp","r")
        jinfo = jf.read().split('\n')[:-2]
        while '' in jinfo:
            jinfo.remove('')
        jf.close()
        os.system("rm jinfo.tmp")
        ncpus = 1
        for l in jinfo:
            if len(l.split('='))>0:
                if l.split('=')[0]=="WorkDir":
                    workdir = l.split('=')[1]
                #if l.split()[0]=="Resource_List.ncpus":
                    #ncpus = int(l.split()[2])
                if l.split()[0].split('=')[0]=="NumNodes":
                    ncpus = int(l.split()[1].split("=")[1])
        ourjob=True
        try:
            job = np.load(workdir+"/job.npy").item()
        except:
            #for nl in range(0,len(jinfo)):
                #l = jinfo[nl]
                #if len(l.split())>0:
                    #if l.split()[0]=="init_work_dir":
                        #try:
                            #workdir = l.split()[2] + jinfo[nl+1].split()[0]
                        #except:
                            #workdir = l.split()[2]
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
