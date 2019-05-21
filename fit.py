
import os
import sys
import numpy as np
import matplotlib.pyplot as plot
from scipy.optimize import curve_fit
import time


def fit_function( x, p1, p2, p3, p4 ):
  return p1 + x*p2 + x*x*p3 + x*x*x*p4


def plot_window(wl_f, center, width):
  mean = np.mean(wl_f, axis=0)
  sigma = np.std(wl_f, axis=0)/np.sqrt((wl_f.shape[0]-1))
  x = np.linspace(0, sigma.shape[0]-1, sigma.shape[0])

  window = np.logical_and(x > center-width-1, x < center+width+1)
  x = np.linspace( -width, width, 2*width+1)
  mean = mean[ window ]
  sigma = sigma[ window ]
  plot.errorbar( x, mean, sigma, fmt='o' , capsize=4 )


def fit_window( wl_f, center, width ):
  sigma = np.std(wl_f, axis=0)/np.sqrt((wl_f.shape[0]-1))
  x = np.linspace(0, sigma.shape[0]-1, sigma.shape[0])

  window = np.logical_and(x > center-width-1, x < center+width+1)
  window = np.logical_and(window, x >= 0)
  sigma = np.array(sigma[ window ])
  
  wl_f = np.array(wl_f[ :, window ])
  x = x[window]

  par = []
  for m in range(wl_f.shape[0]):
    parameters, conv = curve_fit(fit_function, x, wl_f[m], [0,0,0,0], sigma=sigma)
    par.append(parameters)
    diff  = (fit_function(x,*par[0]) - wl_f[m])/sigma
    chisq = (diff*diff).sum()
  return par


def evaluate_fit(x, par):
  fit_evaluated = 0*x
  for m in range(len(par)):
    fit_evaluated += fit_function(x, *par[m])
  fit_evaluated /= len(par)
  return fit_evaluated

def plot_fit(par, center, width):
  x = np.linspace( -width, width, 101)
  plot.plot( x, evaluate_fit(x, par) )


def read_data( nruns ):
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

  wl_f = np.array(wl_f).astype(np.float)
  weight_sums = np.sum(np.exp(wl_f), axis=1)
  energy_correction = np.log(weight_sums)
  wl_f = ( wl_f.transpose() - energy_correction ).transpose()

  return wl_f


def window_smooth( x, wl_f, width ):
  par = fit_window(wl_f, x, width)
  return evaluate_fit(x, par)


def plot_window_fit( nruns, center, width ):
  wl_f = read_data( nruns )

  plot_window(wl_f, center, width)

  par = fit_window(wl_f, center, width)
  plot_fit(par, center, width)
  
  plot.show()


def plot_smoothing( nruns, width, max ):
  wl_f = read_data( nruns )
  wl_f = wl_f[:,:max]

  mean = np.mean(wl_f, axis=0)
  sigma = np.std(wl_f, axis=0)/np.sqrt((wl_f.shape[0]-1))
  x = np.linspace(0, sigma.shape[0]-1, sigma.shape[0])
  plot.errorbar( x, mean, sigma, fmt='o' , capsize=4 )

  wl_f_fit = []
  x = np.linspace(0, sigma.shape[0]-1, sigma.shape[0])
  for point in x:
    print(point)
    wl_f_fit.append(window_smooth( point, wl_f, width ))

  plot.plot( x,  np.array(wl_f_fit))

  plot.show()


if __name__ == "__main__":
  plot_smoothing(20, 4, 40)
  
