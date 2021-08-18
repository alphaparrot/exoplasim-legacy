import os
import numpy as np
from identity import *

class Job:
  def __init__(self,header,args,resource):
    self.args = args.split()
    self.name  = self.args[0]
    self.stat  = self.args[1]
    self.args = self.args[2:]
    self.ncores = int(self.args[0])
    self.queue = self.args[1]
    self.fields = header.split()[3:]
    self.top = os.getcwd()
    self.model = "plasim"
    
    self.parameters = {}
    n=0
    for n in range(0,len(self.fields)):
      self.parameters[self.fields[n]] = self.args[n]
      #print self.fields[n],":",self.args[n]
      
    self.home = resource
    self.jobname = self.name+".cl"
    
  def getID(self):
    os.system("qstat -u "+USER+" > cjobs.tmp")
    jf = open("cjobs.tmp","r")
    jlist = jf.read().split('\n')[5:-1]
    jf.close()
    os.system("rm cjobs.tmp")
    tag = None
    if len(jlist)>0:
      for j in jlist:
        job = j.split()
        name = job[3]
        if name==self.jobname:
          tag = job[0]
          break
    self.tag = tag
    return tag

  def write(self):
    jf = open(self.home+"/job.crwl","w")
    try:
        jt = ' '.join(self.fields)+'\n'+' '.join(self.args)+'\n'+self.name+'\n'+self.tag
    except:
        jt = ' '.join(self.fields)+'\n'+' '.join(self.args)+'\n'+self.name+'\nFAILED'
    jf.write(jt)
    jf.close()

  def kill(self):
    tag = self.getID()
    if tag:
      os.system("qdel "+tag)
      return 1
    else:
      return 0
 
 