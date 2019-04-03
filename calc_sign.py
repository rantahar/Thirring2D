
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
  wl_w = numpy.multiply( 1/numpy.sum(wl_w,1,keepdims=True), wl_w )

  sign = numpy.array( [ -(i%2-0.5)*2 for i in range(101) ] )
  wl_w_s = wl_w*sign

  sign_all = wl_w_s.sum(1)
  sign_mean = sign_all.mean()
  sign_std = sign_all.std()/numpy.sqrt(nruns)

  if(weights):
    wl_w_mean = wl_w.mean(0)
    wl_w_std = wl_w.std(0)/numpy.sqrt(nruns)
    return([sign_mean, sign_std, numpy.array([wl_w_mean, wl_w_std]).T])
  else:
    return([sign_mean, sign_std])
  
if __name__ == "__main__":
  mean, std = calc_sign()
  print(f"{mean} {std}")
