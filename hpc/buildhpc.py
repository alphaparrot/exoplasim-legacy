import sys
import getpass

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
  for k in sys.argv[:]:
      if k.split('=')[0]=="email":
          EMAIL=k.split('=')[1]
      if k.split('=')[0]=="account":
          ACCOUNT=k.split('=')[1]
          
          
  identitypy = ("        \n"+
                'USER = "'+USER+'"        \n'+
                "        \n"+
                'EMAIL = "'+EMAIL+'"        \n'+
                "        \n"+
                'ACCOUNT = "'+ACCOUNT+'"        \n')
  
  with open("identity.py","w") as pyfile:
      pyfile.write(identitypy)
      
    
                
                