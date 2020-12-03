import sys

# Usage: python convertf902dat.py input.txt output.dat

if __name__=="__main__":
  output = []
  with open(sys.argv[1],"r") as infile:
    inp = infile.read().split('\n')
  while "real" not in inp[0].split():
    inp.pop(0)
  while inp[-1]=="":
    inp.pop()
  inp[0] = inp[0].split('/')[1].split()[0].split(',')[:-1]
  inp[-1] = inp[-1].split()[1].split(',')
  for n in range(1,len(inp)-1):
    inp[n] = inp[n].split()[1].split(',')[:-1]
  dat = []
  for line in inp:
    for k in line:
      dat.append(str(float(k)*0.01))
  with open("wvref.txt","r") as wfile:
    wdat = wfile.read().split("\n")
  while wdat[-1]=="":
    wdat.pop()
  for n in range(len(dat)):
    output.append([wdat[n].split()[0],dat[n]])
  outtext = ""
  for o in output:
    outtext += " ".join(o)+"\n"
  with open(sys.argv[2],"w") as outfile:
    outfile.write(outtext)
