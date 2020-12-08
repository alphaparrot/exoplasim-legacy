import os
import glob
import launchset
import sys
import numpy as np
from jobdefs import *
from batch_system import *
from sets import Set


if __name__=="__main__":    
  #Check which resources are in use and compare to our max allotment
  
  #Create auxiliary folders if they don't exist
  os.system("mkdir hopper")
  #os.system("mkdir waitlist")
  
  try:
    f=open(".home","r")
    top = f.read().split('\n')[0]
    f.close()
  except:
    os.system("echo $(pwd)>.home")
    f=open(".home","r")
    top = f.read().split('\n')[0]
    f.close()
    
  
  if "nosub" in sys.argv[:]:
      dryrun=True
  else:
      dryrun=False
  
  #f=open("nnodes.crwl","r")
  #nnodes=float(f.read().split('\n')[0])+1.0 #Chances are this is being called by one of the current jobs!
  #f.close()
  #resources={}
  #for m in MODELS.keys():
    #rf=open("running_"+m+".crwl","r")
    #resources[m] = rf.read().split('\n')[0]
    #resources[m] = resources[m].split()
    ##print resources[m]
    #rf.close()

  njobs=0
  #print "Moving on to normal tasks"
  #print running,nnodes
  capacityflag = False
  openjobs=True
  while openjobs: #We are using less than our full allocation, and the priority list is empty.
    
    #Get next task
    
    f=open("tasks.crwl","r")
    tasks=f.read().split('\n')
    f.close()
    ready=False
    queued = False
    print("Getting job....")
    while not ready:       #Search for a job to run 
      queued = False
      for i in range(0,len(tasks)-1):
        if tasks[i]!='':
          if tasks[i][0]!="#":
            task = tasks[i].split()
            if int(task[1])==0:
              queued=True
              taskname = task[0]
              mark=i                                     #so a plasim job might put us over the limit.
              f=open("tasklog.crwl","a")
              f.write("\nFound job "+' '.join(task))
              f.close()
              break
      if not queued:
          openjobs=False
          ready=False
          print("We have run out of jobs and ARE TRYING TO BREAK")
          break
      ready=False
      f=open("tasklog.crwl","a")
      f.write("\nSearching for header...")
      f.close()
      header=False
      for i in range(mark,-1,-1): #Find header for the task
        if tasks[i][0]=="#":
          header = tasks[i]
          if (len(header.split())-1)==(len(task)): #first header with the right number of args
            break
          else:
            header=False
      if not header:                #None match!
        f=open("tasklog.crwl","a")
        f.write("\nTask "+str(queued)+" header mismatch; skipping")
        f.close()
      else:                         #We found the header. Onward!
        ready=True
        f=open("tasklog.crwl","a")
        f.write("\nTask "+str(queued)+" header found")
        f.close()          
      

    #Engage next task
    if ready:
      
      print("Creating a Job with header\n",header,"\n and arguments \n",' '.join(task)) 
      newjob = Job(header,' '.join(task),-1)       #Collect and organize the job parameters
      
      newjob.home = taskname
      

      task[1] = '1'
      task = ' '.join(task)
      tasks[mark] = task
      tasks='\n'.join(tasks)
      if not Set(task).issubset(Set(" \n")): #Make sure we have something to write!
          f=open("tasks.crwl","w")
          f.write(tasks)
          f.close()
      
      launchset.newtask(newjob,dryrun=dryrun)            #Set up the job and submit it
      np.save(newjob.home+'/job.npy',newjob)
      if not dryrun:
          newjob.getID()
      else:
          newjob.tag='xxxxx.doug'
      newjob.write()
      
      njobs+=1
      
    
  #folders = []
  ##print MODELS.keys()
  #for m in MODELS.keys():
    #fr=open("running_"+m+".crwl","w")
    ##print resources[m]
    ##print '--'
    #fr.write(' '.join(resources[m])+'\n')
    #fr.close()
    #folders+=sorted(glob.glob(m+"/job*/"))
  
  ##Keep track of the current jobs
  
  #htext = "Current job array:\n"
  #for f in folders:
    #jf = open(f+"/job.crwl")
    #jc = jf.read().split('\n')
    #jf.close()
    #htext+=f+" "+jc[3]+"\n"
  #cjb = open("currentjobs.crwl","w")
  #cjb.write(htext)
  #cjb.close()
  #os.system("echo 'Relinquishing control'>>tasklog.crwl")
  #f=open(top+"/inuse.crwl","w")        #Release ownership of crawler.py
  #f.write('0')
  #f.close()
  if dryrun:
    print("%d jobs have been created!"%njobs)
  else:
    print("%d jobs have been created and submitted!"%njobs)