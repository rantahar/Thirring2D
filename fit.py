
import os
import sys
import numpy as np
import matplotlib.pyplot as plot
from scipy.optimize import curve_fit
import time


def fit_function( x, p1, p2, p3 ):
  return p1 + x*p2 + x*x*p3


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

  x = x[window]-center
  wl_f = wl_f[ :, window ]
  sigma = sigma[ window ] * ( 1+np.abs(x)*2/width )

  par = []
  for m in range(wl_f.shape[0]):
    parameters, conv = curve_fit(fit_function, x, wl_f[m], [0,0,0], sigma=sigma)
    par.append(parameters)
    diff  = (fit_function(x,*par[m]) - wl_f[m])/sigma
    chisq = (diff*diff).sum()
  return par


def evaluate_fit(x, par):
  fit_evaluated = []
  for m in range(len(par)):
    fit_evaluated.append(fit_function(x, *par[m]))
  return np.array(fit_evaluated)


def average_fit(x, par):
  fit = evaluate_fit(x, par)
  mean = np.mean( fit, axis=0)
  return mean


def plot_fit(par, center, width):
  x = np.linspace( -width, width, 101)
  y = average_fit(x, par)
  plot.plot( x, y )


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


def window_smooth( x, wl_f, width, eval = 0 ):
  par = fit_window(wl_f, x, width)
  return evaluate_fit(0, par)


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
    value = window_smooth(point, wl_f, width).mean(axis=0)
    wl_f_fit.append(value)

  plot.plot( x,  np.array(wl_f_fit))

  plot.show()


def average_sign( nruns, width, max ):
  wl_f = read_data( nruns )
  wl_f = wl_f[:,:max]

  weights = []
  for x in range(max):
    free_energy = window_smooth( x, wl_f, width )
    weight = np.exp(free_energy) * ( 1 - x%2*2 )
    weights.append(weight)
  weights = np.array(weights).transpose()

  sign = np.sum(weights, axis=1)
  mean = np.mean(sign)
  sigma = np.std(sign)/(np.sqrt(sign.shape[0]-1))
  return [mean, sigma]
  


if __name__ == "__main__":
  if len(sys.argv) > 3 :
    nruns = int(sys.argv[1])
    width = int(sys.argv[2])
    max_sector = int(sys.argv[3])
    if len(sys.argv) > 4:
      do_plot = (sys.argv[4] == 'plot')
    else :
      do_plot = False
  else:
    print("usage: fit.py Nruns window_width max_sector ")

  if do_plot:
    plot_smoothing(nruns, width, max_sector)
  else:
    mean, sigma = average_sign(nruns, width, max_sector)
    print(mean, sigma)

  
