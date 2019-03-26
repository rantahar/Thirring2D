
import os
import subprocess
import random

U = 0.3
m = 0.1
mu = 0.6
X = 32
T = 32

# Write start values for f (just 0)
text_file = open("WL_F_1", "w")
for i in range(101):
  text_file.write("0 0\n")
text_file.close()

for run in range(10):
  seed = random.randint(0,99999999)
  parameter = f'''100000
1
100000
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

  processes = []
  with open(f"parameter_{run}", "r") as inputfile:
    with open(f"output_{run}", "a+") as outputfile:
      p = subprocess.Popen(["./worldline_WL"], stdin=inputfile, stdout=outputfile)
      processes.append(p)

for p in processes:
  p.wait()



