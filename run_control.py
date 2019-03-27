
import os
import sys
import subprocess
import random
import numpy
from time import sleep

from calc_sign import calc_sign

nruns = 160
nnodes = 40

nupdates = 1000000
U = 0.3
m = 0.1
mu = 0.6
X = 32
T = 32

tolerance = 1e-2

if len(sys.argv) == 2 and sys.argv[1] in ["local", "-l"]:
  use_srun = False
else :
  use_srun = True

error = 1
while error > tolerance:
  for run in range(nruns):
    seed = random.randint(0,99999999)
    parameter = f'''{nupdates}
1
{nupdates}
{m}
{U}
{mu}
{seed}
0.01
1000
WL_F_{run}'''
  
    text_file = open(f"parameter_{run}", "w")
    text_file.write(parameter)
    text_file.close()

  n_running = 0
  processes = []
  for run in range(nruns):
    with open(f"parameter_{run}", "r") as inputfile:
      with open(f"output_{run}", "a+") as outputfile:
        if(use_srun):
          p = subprocess.Popen(["srun","-n1","-N1","./worldline_WL"], stdin=inputfile, stdout=outputfile)
        else:
          p = subprocess.Popen(["./worldline_WL"], stdin=inputfile, stdout=outputfile)
        processes.append(p)
        n_running += 1
    while n_running == nnodes:
      sleep(0.05)
      for p in processes:
        if p.poll() is not None:
          n_running -= 1

  
  for p in processes:
    p.wait()
  
  
  sign_mean, sign_std = calc_sign(nruns)
  
  print(sign_mean,sign_std)
  sys.stdout.flush()
  
  error = sign_std/sign_mean
  
