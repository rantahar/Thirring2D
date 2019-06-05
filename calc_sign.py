
import os
import sys
import numpy
import time

nruns = 160

def calc_sign( nruns=nruns, weights=False ):
  wl_f = []
  for run in range(nruns):
    file_open = False
    while not file_open:
      try:
        text_file = open(f"WL_F_{run}", "r")
        file_open = True
      except:
        time.sleep(1)
        continue
    wl_f_n = []
    for line in text_file:
      values = line.split()
      wl_f_n.append(values[0])
    wl_f.append(wl_f_n[:-1])
    text_file.close()

  wl_f = numpy.array(wl_f).astype(numpy.float)
  wl_w = numpy.exp(wl_f)
  wl_w = ( wl_w.transpose()/wl_w.sum(1) ).transpose()

  sign = numpy.array( [ -(i%2-0.5)*2 for i in range(wl_w.shape[1]) ] )
  wl_w_s = wl_w*sign

  sign_all = wl_w_s.sum(1)
  sign_mean = sign_all.mean()
  sign_std = sign_all.std()/numpy.sqrt(nruns-1)

  if(weights):
    wl_w_mean = wl_w.mean(0)
    wl_w_std = wl_w.std(0)/numpy.sqrt(nruns-1)
    return([sign_mean, sign_std, numpy.array([wl_w_mean, wl_w_std]).T])
  else:
    return([sign_mean, sign_std])
  
if __name__ == "__main__":
  if len(sys.argv) > 1 :
    nruns = int(sys.argv[1])
  
  mean, std = calc_sign( nruns )
  print(f"{mean} {std}")
