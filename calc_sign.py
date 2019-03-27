
import os
import sys
import numpy

nruns = 160

def calc_sign( nruns=nruns ):
  wl_f = []
  for run in range(nruns):
    text_file = open(f"WL_F_{run}", "r")
    wl_f_n = []
    for line in text_file:
      values = line.split()
      wl_f_n.append(values[0])
    wl_f.append(wl_f_n)
    text_file.close()

  wl_f = numpy.array(wl_f).astype(numpy.float)
  wl_w = numpy.exp(wl_f)
  wl_w = numpy.multiply( 1/numpy.sum(wl_w,1,keepdims=True), wl_w )

  print(wl_w.mean(0))

  sign = numpy.array( [ -(i%2-0.5)*2 for i in range(100) ] )
  wl_w *= sign

  sign_all = wl_w.sum(1)
  sign_mean = sign_all.mean()
  sign_std = sign_all.std()

  return([sign_mean,sign_std])
  
if __name__ == "__main__":
  mean, std = calc_sign()
  print(f"{mean} {std}")