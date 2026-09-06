import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import interpolate
from collections import defaultdict

### constants
c = 2.998E8   # speed of light in m/s
spectoGHz = c/0.01 * 1e-9   # convert spectrocsopic units of 1/cm to GHz

### functions
def read_transmission_csv(files, finds):
  # read filter transmission data from csv files and return dictionary of transmission data
  transmission = {}
  for ff, filter_file in enumerate(files):
    with open(filter_file, mode ='r', encoding='utf-8-sig')as file:
      csv_reader = csv.DictReader(file)
      first_line = True
      for row in csv_reader:
        if first_line:
          ### meta data
          serial = row['serial']
          filters = np.array(list(row.keys()))[finds]   # filter names are in these columns
          transmission[serial] = {}
          transmission[serial]['description'] = row['description']   # description of measurement
          # transmission[serial]['filters']     = filters   # some metadata that i don't use

          ### initialize data arrays
          for ii in np.arange(7):   # data is split into seven columns, not for any particular reason
            transmission[serial]['freq'+str(ii+1)+' [GHz]'] = np.array([float(row['f'+str(ii+1)])] if row['f'+str(ii+1)] != '' else np.nan)   # freq in GHz
            transmission[serial][filters[ii]]               = np.array([float(row[filters[ii]])] if row[filters[ii]] != '' else np.nan)   # transmission
          first_line = False
        else:   # append arrays
          for ii in np.arange(7):   # data is split into seven columns, not for any particular reason
            transmission[serial]['freq'+str(ii+1)+' [GHz]'] = np.append(transmission[serial]['freq'+str(ii+1)+' [GHz]'], float(row['f'+str(ii+1)]) if row['f'+str(ii+1)] != '' else np.nan)   # some files have fewer than 7 frequencies, so add NaN if f7 is empty
            transmission[serial][filters[ii]]               = np.append(transmission[serial][filters[ii]], float(row[filters[ii]]) if row[filters[ii]] != '' else np.nan)   # some files have fewer than 7 frequencies, so add NaN if transmission is empty
    freqs = np.concatenate((transmission[serial]['freq1 [GHz]'], transmission[serial]['freq2 [GHz]'], transmission[serial]['freq3 [GHz]'], transmission[serial]['freq4 [GHz]'], transmission[serial]['freq5 [GHz]'], transmission[serial]['freq6 [GHz]'], transmission[serial]['freq7 [GHz]']))
    tmiss = np.concatenate((transmission[serial][filters[0]],    transmission[serial][filters[1]],    transmission[serial][filters[2]],    transmission[serial][filters[3]],    transmission[serial][filters[4]],    transmission[serial][filters[5]],    transmission[serial][filters[6]]))
    if ff in [8,9,10,11,12]:   # for DSIR filters, convert frequencies from 1/cm to GHz
        freqs *= spectoGHz
    transmission[serial]['freqs_raw'] = freqs   # combined frequency array
    transmission[serial]['tmiss_raw'] = tmiss   # combined transmission array

  return transmission, serial, filters

def remove_yval_nans(xvals, yvals):   # remove x values and y values where y is nan
  xvals_no_nan = xvals[~np.isnan(yvals)]
  yvals_no_nan = yvals[~np.isnan(yvals)]
  return xvals_no_nan, yvals_no_nan

def find_duplicates(xvals):   # find duplicate x values and return dictionary of x value to list of indices where it occurs
  tally = defaultdict(list)
  for i,item in enumerate(xvals):
    tally[item].append(i)
  return tally

def remove_duplicates(xvals, yvals):   # remove second instance of duplicate x values and corresponding y values
  sorted_xvals, sortinds = np.unique(xvals, return_index=True)   # get unique x values and indices to sort y values
  sorted_yvals = yvals[sortinds]   # y values corresponding to unique x values
  return sorted_xvals, sorted_yvals

def average_duplicates(xvals, yvals):   # find duplicate x values and replace corresponding y values with their average, return unique sorted x and y arrays
  ave_yvals = yvals.copy()
  duplicate_xvals = []
  tally = find_duplicates(xvals)
  for key, locs in tally.items():
    if len(locs)>1:
      ave_yvals[locs] = np.mean(yvals[locs])
      duplicate_xvals.append(key)
  sorted_xvals, sortinds = np.unique(xvals, return_index=True)   # get unique x values and indices to sort y values
  sorted_yvals = ave_yvals[sortinds]   # y values corresponding to unique x values
  return sorted_xvals, sorted_yvals

def clean_and_interp(xvals, yvals, xnew, average=False, check_interp=False):   # average y values of duplicate x values, then interpolate to new x values
  xvals_nonan, yvals_nonan   = remove_yval_nans(xvals, yvals)   # remove x and y where y=nan
  if average:   # average duplicate values
    sorted_xvals, sorted_yvals = average_duplicates(xvals_nonan, yvals_nonan)   # average y values of duplicate x values
  else:   # remove duplicate values
    sorted_xvals, sorted_yvals = remove_duplicates(xvals_nonan, yvals_nonan)   # remove second instance of duplicate x values and corresponding y values

  b = interpolate.make_interp_spline(sorted_xvals, sorted_yvals, k=3)   # k=1 is linear, k=3 is cubic spline, k=4 is quartic spline
  yvals_interp = b(xnew, extrapolate=False)

  yvals_interp[yvals_interp>1] = 1   # max transission is 1

  if check_interp:   # plot for checking interpolation
    plt.figure()
    plt.plot(xvals,        yvals,        '.', markersize=5, alpha=0.5, label='Raw Data')
    plt.plot(sorted_xvals, sorted_yvals, '.', markersize=5, alpha=0.5, label='Cleaned Data')
    plt.plot(xnew,         yvals_interp, linewidth=2, color='k', label='Interpolated')
    plt.legend(markerscale=2)

    # plt.figure()
    # hist, bins, patches = plt.hist(xvals_nonan, bins=50, alpha=0.5)
    # plt.hist(sorted_xvals, bins=bins, alpha=0.5, label='Cleaned Data')
    # plt.title('Removed {} out of {} duplicate x values ({}\\%)'.format(len(xvals_nonan)-len(sorted_xvals), len(xvals_nonan), round((len(xvals_nonan)-len(sorted_xvals))/len(xvals_nonan)*100, 2)))
    plt.pause(0.5); plt.show()

  return sorted_xvals, sorted_yvals, xnew, yvals_interp
