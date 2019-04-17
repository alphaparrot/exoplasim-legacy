import sys
import getpass
import os

'''
Usage: python buildhpc.py option1=arg1 option2=arg2

Options:
    email=<you@mail>
    account=<allocation>
    
Both of these are optional if your job submission system doesn't require them.

'''

if __name__=="__main__":
  USER = getpass.getuser()
  EMAIL = ""
  ACCOUNT = ""
  GCCMOD = "gcc/4.9.1"             
  PYTHOMOD = "python/2.7.9"          
  INTELMOD = "intel/intel-17"        
  MPIMOD = "openmpi/2.0.1-intel-17"
  for k in sys.argv[:]:
      if k.split('=')[0]=="email":
          EMAIL=k.split('=')[1]
      if k.split('=')[0]=="account":
          ACCOUNT=k.split('=')[1]
      if k.split('=')[0]=="gcc":
          GCCMOD = "gcc/"+k.split('=')[1]
      if k.split('=')[0]=="python":
          PYTHONMOD = "python/"+k.split('=')[1]
      if k.split('=')[0]=="intel":
          INTELMOD = "intel/"+k.split('=')[1]
      if k.split('=')[0]=="openmpi":
          MPIMOD = "openmpi/"+k.split('=')[1]
      if k=="automod":
          modules = os.environ['LOADEDMODULES'].split(os.pathsep)
          for m in modules:
              if m.split('/')[0]=="gcc":
                  GCCMOD = m
              elif m.split('/')[0]=='python':
                  PYTHONMOD = m
              elif m.split('/')[0]=="intel":
                  INTELMOD = m
              elif m.split('/')[0]=="openmpi":
                  MPIMOD = m
          
          
          
          
  identitypy = ("        \n"+
                'USER = "'+USER+'"        \n'+
                "        \n"+
                'EMAIL = "'+EMAIL+'"        \n'+
                "        \n"+
                'ACCOUNT = "'+ACCOUNT+'"        \n'+
                "        \n"+
                'GCCMOD = "'+GCCMOD+'"          \n'+
                "        \n"+
                'PYTHONMOD = "'+PYTHONMOD+'"          \n'+
                "        \n"+
                'INTELMOD = "'+INTELMOD+'"          \n'+
                "        \n"+
                'MPIMOD = "'+MPIMOD+'"          \n'+
                "        \n")
  
  with open("identity.py","w") as pyfile:
      pyfile.write(identitypy)
      
    
                
                