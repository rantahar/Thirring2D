
import os
import sys
import numpy

nruns = 160

def calc_sectors( nruns=nruns ):
  wl_f = []
  for run in range(nruns):
    text_file = open(f"WL_F_{run}", "r")
    wl_f_n = []
    for line in text_file:
      values = line.split()
      wl_f_n.append(values[0])
    wl_f.append(wl_f_n[:-1])
    text_file.close()
  
  nmeas = len(wl_f[-1])
  wl_f = [ wl_f_n[0:nmeas-1] for wl_f_n in wl_f ]

  wl_f = numpy.array(wl_f).astype(numpy.float)
  
  mean = wl_f.mean(0)
  std = wl_f.std(0)/numpy.sqrt(nruns-1)

  mean -= numpy.log(numpy.sum(numpy.exp(mean)))

  return([mean,std])
  
if __name__ == "__main__":
  if len(sys.argv) > 1 :
    nruns = int(sys.argv[1])

  mean, std = calc_sectors()
  for m, d in zip(mean,std):
    print(f"{m} {d}")