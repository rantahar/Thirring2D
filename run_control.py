
import os
import sys
import subprocess
import random
import numpy
from time import sleep

from calc_sign import calc_sign

nruns = 20
nnodes = 3

nupdates = 1000
nmeasure = 1
U = 0
m = 0.1
mu= 0
step = 10
init_nstep = 1000

tolerance = 0.01
sector_tolerance = 0.001*tolerance
min_accuracy = sector_tolerance
max_iterations = 1000

seed = random.randint(0,99999999)


if len(sys.argv) == 2 and sys.argv[1] in ["local", "-l"]:
  use_srun = False
else :
  use_srun = True


def new_parameter_file( run ):
  global seed
  seed = seed+1
  parameter = f'''{nupdates}
{nmeasure}
{int(nupdates/nmeasure)}
{m}
{U}
{mu}
{seed}
.config_{run}
{step}
{init_nstep}
WL_F_{run}'''
  
  text_file = open(f".parameter_{run}", "w")
  text_file.write(parameter)
  text_file.close()


error = 1
for i in range(max_iterations):
  for run in range(nruns):
    new_parameter_file( run )

  processes = []
  for run in range(nruns):
    with open(f".parameter_{run}", "r") as inputfile:
      with open(f"output_{run}", "a+") as outputfile:
        if(use_srun):
          p = subprocess.Popen(["srun","-n1","-N1","./worldline_WL"], stdin=inputfile, stdout=outputfile)
        else:
          p = subprocess.Popen(["./worldline_WL"], stdin=inputfile, stdout=outputfile)
        processes.append(p)
    while len(processes) >= nnodes:
      sleep(0.5)
      for p in processes:
        if p.poll() is not None:
          processes.remove(p)
  
  for p in processes:
    p.wait()
  
  sign_mean, sign_std, weights = calc_sign(nruns, weights=True)
  sector_done = (weights[:,0]+weights[:,1])/abs(sign_mean) < sector_tolerance
  sector_done = numpy.logical_and( sector_done, (weights[:,0]+weights[:,1]) < min_accuracy )
  sector_done = numpy.append(sector_done, [True,True])
  
  for run in range(nruns):
    sector = 0
    text = ''
    weight_file = open(f"WL_F_{run}", "r")
    for line in weight_file:
      values = line.split()
      if len(values) == 3:
        if sector_done[sector]:
          text += f"{values[0]} {values[1]} 0\n"
        else:
          text += f"{values[0]} {values[1]} 1\n"
        sector += 1
    weight_file.close()
    weight_file = open(f"WL_F_{run}", "w")
    weight_file.write(text)
    weight_file.close()
      
  print(sign_mean,sign_std)
  sys.stdout.flush()
  
  error = sign_std/abs(sign_mean)
  if( error < tolerance ):
    break
  
