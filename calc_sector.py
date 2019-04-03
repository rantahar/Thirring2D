
import os
import sys
import numpy

nruns = 160

sector = 1
if len(sys.argv) > 1 :
  sector = sys.argv[1]

def calc_sector( nruns=nruns ):
  wl_f = []
  for run in range(nruns):
    text_file = open(f"output_{run}", "r")
    wl_f_n = []
    for line in text_file:
      if f"SECTOR {sector} " in line:
        values = line.split()
        wl_f_n.append(values[3])
    wl_f.append(wl_f_n)
    text_file.close()
  
  nmeas = len(wl_f[-1])
  wl_f = [ wl_f_n[0:nmeas-1] for wl_f_n in wl_f ]

  wl_f = numpy.array(wl_f).astype(numpy.float)
  
  mean = wl_f.mean(0)
  std = wl_f.std(0)/numpy.sqrt(nruns)

  return([mean,std])
  
if __name__ == "__main__":
  mean, std = calc_sector()
  for m, d in zip(mean,std):
    print(f"{m} {d}")