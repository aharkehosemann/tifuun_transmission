# Angi Harke-Hosemann 2026/06
#
# calculates total transmission of TIFUUN filters from Cardiff's filter measurements
#
# meta data
# 1  = FP 3668 ARC = LPE F1? - 50 K thick IR blocker 90-360 GHz
# 2  = FP 3667     = LPE F3? - 1K filter ~360 GHz
# 3  = FP 3662     = LPE F4 360 GHz? - 300 mK filter 180-360 GHz
# 4  = FP 3572     = LPE F4 180 GHz? - 300 mK filter 990-180 GHz
# 5  = FP 3450     = LPE F2? - 4K LPE filter 500 GHz
# 6  = FP 3429 ARC = UHMWPE ARC window at 300 K?
# 7  = FP 3083     = LPE F4 360 GHz as well? - 300 mK filter 180-360 GHz
# 8  = FP 3083     = LPE F4 360 GHz as well? - 300 mK filter 180-360 GHz
# 9  = DSIR5 (phd10) 50 K
# 10 = DSIR5 (phd10) 50 K
# 11 = DSIR3&4 (phd8) 140&50 K
# 12 = DSIR1&2 (phd4) 300&140 K
# 13 = DSIR1&2 (phd4) 300&140 K
#
# duplicate values removed
# 1: 216/7007 data points removed = 3% duplicate x values in 20-330 GHz range, most centered around 150 GHz
# 2: 216/7007 data points removed = 3% duplicate x values in 20-330 GHz range, most centered around 150 GHz
# 3: 194/4004 data points removed = <5% duplicate x values in 20-330 GHz range, most centered around 150 GHz
# 4: 156/7007 data points removed = 2% duplicate x values, some around 30 GHz and some around 150 GHz
# 5: 0/1001 data points removed = 0% duplicate x values
# 6: 0/1001 data points removed = 0% duplicate x values
# 7: 0/1001 data points removed = 0% duplicate x values
# 8: 0/1001 = 0% duplicate x values
# 9: ?/1001 = ?% duplicate x values
# 10: ?/1001 = ?% duplicate x values
# 11: ?/1001 = ?% duplicate x values
# 12: ?/1001 = ?% duplicate x values
# 13: ?/1001 = ?% duplicate x values
#
# notes
# TIFUUN bands are 130–178 GHz and 195–319 GHz, or 90 - 360 GHz according to Akira's SPIE proceedings
# 1/cm = 30 GHz, 130 GHz ~ 4.3 1/cm, 178 GHz ~ 5.9 1/cm, 195 GHz ~ 6.5 1/cm, 319 GHz ~ 10.6 1/cm
# phd10 lowest frequency measurement starts at 600 GHz, 3900
# phd8 lowest frequency measurement starts at 20 1/cm = 600 GHz
# phd4 lowest frequency measurement starts at 4 1/cm = 120 GHz
# F1, 3, and some 4 measurements only go to 330 GHz
# F2, AR window, and some F4 measurements go to 1000-1200 GHz
# 8 (F4) to 1130, 7 (F4) to 1200, 6 (AR window) to 1080, 5 (F2) to 1050, 1-4 (F1,3,4) to 330
# three measurements of F4 180-360GHz, '3' measurement is highest fidelity in band of interest (65-330 GHz). Others could be used to extend that range, but we are still limited to 330 GHz by F1 and F3 measurements
# sept 2026: something is different about 50K transmission around 120 GHz

import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import interpolate
from collections import defaultdict

### user configuration and analysis options
root_dir                = '/Users/angi/tifuun/transmission/tmiss_measurements/'
ind_msmts_to_plot       = [2, 5, 1, 6]   # indices of raw filter measurements to plot, set to [] for no raw measurement plots
ave_dups                = False   # average duplicate frequency measurements when interpolating to new frequencies, otherwise second instance of duplicates will be removed
check_dsir_clean        = True   # plot for checking cleaning of DSIR filter measurements
check_interp            = False   # plot for checking interpolation of cleaned data to new frequencies
freq_min = 0; freq_max = 330   # GHz, frequency range to interpolate total transmission

### plot settings
font = {'family' : 'serif', 'weight' : 'normal', 'size'   : 18}
ticks = {'major.size'   : '5', 'labelsize'     : '16', 'minor.visible' : False}
plt.rc('text', usetex=True); plt.rc('font', **font)
plt.rc('xtick', **ticks);    plt.rc('ytick', **ticks)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}'
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['figure.dpi'] = 100; plt.rcParams['savefig.dpi'] = 300   # higher resolution plots in interactive notebook and when saving pngs
figsize = (10, 6)   # width, height

### constants
c = 2.998E8   # speed of light in m/s
spectoGHz = c/0.01 * 1e-9   # convert spectrocsopic units of 1/cm to GHz

### filter measurement files
filter_file1  = root_dir+'3668data.csv'      # FP 3668 ARC, LPE F1? - 50 K thick IR blocker 90-360 GHz
filter_file2  = root_dir+'3667data.csv'      # FP 3667,     LPE F3? - 1K filter ~360 GHz
filter_file3  = root_dir+'3662data.csv'      # FP 3662,     LPE F4 360 GHz? - 300 mK filter 180-360 GHz
filter_file4  = root_dir+'3572data.csv'      # FP 3572,     LPE F4 180 GHz? - 300 mK filter 990-180 GHz
filter_file5  = root_dir+'S3431R5.csv'       # FP 3450,     LPE F2? - 4K LPE filter 500 GHz
filter_file6  = root_dir+'S3424R7.csv'       # FP 3429 ARC, UHMWPE ARC window at 300 K?
filter_file7  = root_dir+'S3400R11.csv'      # FP 3083, also LPE F4 360 GHz? - 300 mK filter 180-360 GHz
filter_file8  = root_dir+'S3400R13.csv'      # FP 3083, also LPE F4 360 GHz? - 300 mK filter 180-360 GHz
filter_file9  = root_dir+'C0257_3.csv'       # DSIR5 (phd10) 50 K, one measurement, use for 225+ 1/cm
filter_file10 = root_dir+'T1889R7.csv'       # DSIR5 (phd10) 50 K, another measurement, use for 0-225 1/cm
filter_file11 = root_dir+'T1889R5.csv'       # DSIR3&4 (phd8) 140&50 K
filter_file12 = root_dir+'phd4_combined.csv' # DSIR1&2 (phd4) 300&140 K, one measurement, Carole combined low and mid-frequency measurements
filter_file13 = root_dir+'M24779.csv'        # DSIR1&2 (phd4) 300&140 K, another measurement, Carole suggests switching to this around 330 1/cm
filter_files  = [filter_file1, filter_file2, filter_file3, filter_file4, filter_file5, filter_file6, filter_file7, filter_file8, filter_file9, filter_file10, filter_file11, filter_file12, filter_file13]

### read data files
serial_to_filter = {'FP3429ARC ':    'AR Coated Window (300K)',     'FP3668ARC':     'LPE F1 - 90-360GHz (50K)',      # files 6  and 1
                    'FP3450':        'LPE F2 - 500GHz (4K)',        'FP3667':        'LPE F3 - 360GHz (1K)',          # files 5  and 2
                    'FP3662':        'LPE F4 - 180-360GHz (300mK)', 'FP3572':        'LPE F4 - 990-180GHz (300mK)',   # files 3  and 4
                    'FP3083_10p5cm': 'LPE F4 - 180-360GHz (300mK)', 'FP3083_11p5cm': 'LPE F4 - 180-360GHz (300mK)',   # files 7  and 8
                    'C0257_3':       'DSIR5 (50K)',                 'TI889R7':       'DSIR5 (50K)',                   # files 9  and 10
                    'T1889R5':       'DSIR3\\&4 (140\\&50K)',       'M24779':        'DSIR1\\&2 (300\\&140K)',        # files 11 and 13
                    'PHD4_combined': 'DSIR1\\&2 (300\\&140K)',      'PHD4':          'DSIR1\\&2 (300\\&140K)',        # file 12
                    'PHD8':          'DSIR3\\&4 (140\\&50K)',       'PHD10':         'DSIR5 (50K)'}                   # combined later
filter_cols   = [6, 10, 14, 18, 22, 26, 30]   # filter names are in these columns

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

tmiss_all, serial_all, filter_all = read_transmission_csv(filter_files, filter_cols)

serials    = list(tmiss_all)
freq9_raw  = tmiss_all[serials[8]]['freqs_raw']; freq10_raw  = tmiss_all[serials[9]]['freqs_raw']; freq11_raw  = tmiss_all[serials[10]]['freqs_raw']; freq12_raw  = tmiss_all[serials[11]]['freqs_raw']; freq13_raw  = tmiss_all[serials[12]]['freqs_raw']
tmiss9_raw = tmiss_all[serials[8]]['tmiss_raw']; tmiss10_raw = tmiss_all[serials[9]]['tmiss_raw']; tmiss11_raw = tmiss_all[serials[10]]['tmiss_raw']; tmiss12_raw = tmiss_all[serials[11]]['tmiss_raw']; tmiss13_raw = tmiss_all[serials[12]]['tmiss_raw']

### plot raw filter measurements
if len(ind_msmts_to_plot)>0:
  fig, ax = plt.subplots(figsize=figsize, layout='tight')
  for mm in ind_msmts_to_plot:
    plt.plot(tmiss_all[serials[mm-1]]['freqs_raw'], tmiss_all[serials[mm-1]]['tmiss_raw'], '.', markersize=4, label=serial_to_filter[serials[mm-1]], alpha=0.7)
  plt.xlim(0,1000); plt.xlabel('Frequency [GHz]')   # GHz
  plt.legend(markerscale=2, loc='upper right')

  secax = ax.secondary_xaxis('top', functions=(lambda x: x/spectoGHz, lambda x: x*spectoGHz))
  secax.set_xlabel('Wave Number [1/cm]')
  plt.grid(linestyle = '--', which='both', linewidth=0.5)
  plt.ylabel('Filter Transmission')
  plt.ylim(0,1)

### process data
# combine low and high freq DSIR 10 um measurements
phd10_hfinds = np.where(freq9_raw>=225*spectoGHz)[0]   # high freq = large wavenumber
phd10_lfinds = np.where(freq10_raw<225*spectoGHz)[0]   # low freq = small wavenumber
phd10c_freq  = np.concatenate((freq9_raw[phd10_hfinds],  freq10_raw[phd10_lfinds]))
phd10c_tmiss = np.concatenate((tmiss9_raw[phd10_hfinds], tmiss10_raw[phd10_lfinds]))

# combine low and high freq DSIR 4 um measurements
phd4_hfinds = np.where(freq13_raw>=310*spectoGHz)[0]   # high freq = large wavenumber
phd4_lfinds = np.where(freq12_raw<310*spectoGHz)[0]    # low freq = small wavenumber
phd4_freq   = np.concatenate((freq12_raw[phd4_lfinds],  freq13_raw[phd4_hfinds]))
phd4_tmiss  = np.concatenate((tmiss12_raw[phd4_lfinds], tmiss13_raw[phd4_hfinds]))

# tack on DSIR 4 um measurements below 600 GHz to 8 um and 10 um measurements
phd4_maskinds = np.where(phd4_tmiss<min(np.concatenate((freq11_raw, phd10c_freq))))[0]
phd10_freq    = np.concatenate((phd4_freq[phd4_maskinds],  phd10c_freq))
phd10_tmiss   = np.concatenate((phd4_tmiss[phd4_maskinds], phd10c_tmiss))
phd8_freq     = np.concatenate((phd4_freq[phd4_maskinds],  freq11_raw))
phd8_tmiss    = np.concatenate((phd4_tmiss[phd4_maskinds], tmiss11_raw))

# set transmission of DSIR filters to unity below 120 GHz
unitymask_freqs = np.linspace(0, min(phd4_freq))
unitymask_tmiss = np.ones_like(unitymask_freqs)
phd4_freq   = np.concatenate((unitymask_freqs, phd4_freq));  phd4_tmiss  = np.concatenate((unitymask_tmiss, phd4_tmiss))
phd10_freq  = np.concatenate((unitymask_freqs, phd10_freq)); phd10_tmiss = np.concatenate((unitymask_tmiss, phd10_tmiss))
phd8_freq   = np.concatenate((unitymask_freqs, phd8_freq));  phd8_tmiss  = np.concatenate((unitymask_tmiss, phd8_tmiss))

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

  tck = interpolate.splrep(sorted_xvals, sorted_yvals, s=0, k=3)
  yvals_interp = interpolate.BSpline(*tck, extrapolate=False)(xnew)
  yvals_interp[yvals_interp>1] = 1   # max transission is 1

  if check_interp:   # plot for checking interpolation
    plt.figure()
    plt.plot(xvals,        yvals,        '.', markersize=5, alpha=0.5, label='Raw Data')
    plt.plot(sorted_xvals, sorted_yvals, '.', markersize=5, alpha=0.5, label='Cleaned Data')
    plt.plot(xnew,         yvals_interp, '.', markersize=3, color='k', label='Interpolated')
    plt.grid(linestyle = '--', which='both', linewidth=0.5)
    plt.legend(markerscale=2)

    plt.figure()
    hist, bins, patches = plt.hist(xvals_nonan, bins=50, alpha=0.5)
    plt.hist(sorted_xvals, bins=bins, alpha=0.5, label='Cleaned Data')
    plt.title('Removed {} out of {} duplicate x values ({}\\%)'.format(len(xvals_nonan)-len(sorted_xvals), len(xvals_nonan), round((len(xvals_nonan)-len(sorted_xvals))/len(xvals_nonan)*100, 2)))
    plt.show()

  return sorted_xvals, sorted_yvals, xnew, yvals_interp

# average or remove duplicate frequency measurements and interpolate to new frequencies
freqs_interp = np.linspace(freq_min, freq_max, num=1000)   # GHz, frequencies to interpolate total transmission

for ff in np.arange(8):   # for each filter measurement, clean and interpolate to new frequencies
  freq_sorted,  tmiss_sorted,  freq_interp, tmiss_interp = clean_and_interp(tmiss_all[serials[ff]]['freqs_raw'],  tmiss_all[serials[ff]]['tmiss_raw'],  freqs_interp, average=ave_dups, check_interp=check_interp)   # removed 216/7007 = 3% duplicate x values in 20-330 GHz range, most centered around 150 GHz
  tmiss_all[serials[ff]]['freqs_clean'],  tmiss_all[serials[ff]]['tmiss_clean'],  tmiss_all[serials[ff]]['freqs_interp'], tmiss_all[serials[ff]]['tmiss_interp'] = freq_sorted,  tmiss_sorted, freq_interp, tmiss_interp

# freq1 removed 216/7007 = 3% duplicate x values in 20-330 GHz range, most centered around 150 GHz
# freq2 removed 216/7007 = 3% duplicate x values in 20-330 GHz range, most centered around 150 GHz
# freq3 removed 194/4004 = <5% duplicate x values in 20-330 GHz range, most centered around 150 GHz
# freq4 removed 156/7007 = 2% duplicate x values, some around 30 GHz and some around 150 GHz
# freq5 removed 0/1001 = 0% duplicate x values
# freq6 removed 0/1001 = 0% duplicate x values
# freq7 removed 0/1001 = 0% duplicate x values
# freq8 removed 0/1001 = 0% duplicate x values

## use masked arrays for low frequency measurements of DSIR fitlers
phd4_freq_sorted,  phd4_tmiss_sorted,  phd4_freq_interp,  phd4_tmiss_interp  = clean_and_interp(phd4_freq,   phd4_tmiss,   freqs_interp, check_interp=check_interp)   # removed 1/12490 = <<1% duplicate x values, potentially in range but hard to tell because measurements go to 150 THz
phd10_freq_sorted, phd10_tmiss_sorted, phd10_freq_interp, phd10_tmiss_interp = clean_and_interp(phd10_freq,  phd10_tmiss,  freqs_interp, check_interp=check_interp)   # removed 1/26090 = <<1% duplicate x values, probably inherited from phd4
phd8_freq_sorted,  phd8_tmiss_sorted,  phd8_freq_interp,  phd8_tmiss_interp  = clean_and_interp(phd8_freq,   phd8_tmiss,   freqs_interp, check_interp=check_interp)   # removed 1/15093 = <<1% duplicate x values, probably inherited from phd4

if check_dsir_clean:
  ### DSIR Filters - measurements in spectroscopic units of 1/cm
  # # PHD 4, 8, and 10 individually in GHz
  # fig, ax = plt.subplots(figsize=figsize, layout='tight')
  # plt.plot(freq9_raw,      tmiss9_raw,      'o', markersize=4, alpha=0.3, label=serials[8])   # DSIR5 (phd10) 50 K, one measurement
  # plt.plot(freq10_raw,     tmiss10_raw,     'o', markersize=4, alpha=0.3, label=serials[9])   # DSIR5 (phd10) 50 K, another measurement
  # plt.plot(phd10_freq, phd10_tmiss, '.', markersize=2, alpha=0.8, label='PHD10 Combined', color='k')   # DSIR5 (phd10) 50 K, combined measurement
  # plt.plot(freq12_raw,     tmiss12_raw,     'o', markersize=4, alpha=0.3, label='PHD4 Low \\& Mid Freq')   # DSIR1&2 (phd4) 300&140 K, one measurement, Carole combined low and mid-frequency measurements
  # plt.plot(freq13_raw,     tmiss13_raw,     'o', markersize=4, alpha=0.3, label=serials[12])    # DSIR1&2 (phd4) 30OTH&14０ K, another measurement, Carole suggests switching to this around 33０ １/cm = ９９００ THz
  # plt.plot(phd4_freq,  phd4_tmiss,  '.', markersize=2, alpha=0.8, label='PHD4 Combined', color='k')   # DSIR1&2 (phd4) 300&140 K, combined measurement
  # plt.plot(freq11_raw,     tmiss11_raw,     'o', markersize=4, alpha=0.3, label='PHD8')
  # plt.legend(markerscale=2, loc='upper right')
  # plt.xlim(0,1000); plt.xlabel('Frequency [GHz]')   # GHz
  # secax = ax.secondary_xaxis('top', functions=(lambda x: x/spectoGHz, lambda x: x*spectoGHz))

  # all DSIR in GHz
  fig, ax = plt.subplots(figsize=figsize, layout='tight')
  plt.plot(phd4_freq,  phd4_tmiss,  '.', markersize=3, alpha=0.7, label='PHD4')    # DSIR1&2 (phd4) 300&140 K, combined measurements
  plt.plot(freq11_raw,     tmiss11_raw,     '.', markersize=3, alpha=0.7, label='PHD8')    # DSIR3&4 (phd8) 140&50 K
  plt.plot(phd10_freq, phd10_tmiss, '.', markersize=3, alpha=0.7, label='PHD10')   # DSIR5 (phd10) 50 K, combined measurements
  plt.xlim(0,1000); plt.xlabel('Frequency [GHz]')   # GHz
  plt.legend(markerscale=5, loc='lower left')
  secax = ax.secondary_xaxis('top', functions=(lambda x: x/spectoGHz, lambda x: x*spectoGHz))

  plt.ylabel('Filter Transmission')
  secax.set_xlabel('Wave Number [1/cm]')
  plt.grid(linestyle = '--', which='both', linewidth=0.5)
  plt.ylim(0,1)
  plt.tight_layout()

### total transmission to 50 K, 4 K, 1 K, and 300 mK
AR_tmiss        = tmiss_all[serials[5]]['tmiss_interp']
F1_tmiss        = tmiss_all[serials[0]]['tmiss_interp']
F2_tmiss        = tmiss_all[serials[4]]['tmiss_interp']
F3_tmiss        = tmiss_all[serials[1]]['tmiss_interp']
F4_180GHz_tmiss = tmiss_all[serials[3]]['tmiss_interp']
F4_360GHz_tmiss = tmiss_all[serials[2]]['tmiss_interp'] #  ('2' has smallest freq range but highest fidelity)
DSIR1_tmiss     = phd10_tmiss_interp
DSIR2_tmiss     = phd8_tmiss_interp
DSIR4_tmiss     = phd4_tmiss_interp

# # to 50 K: AR Window ('6'), DSIR Filters 1-5 (1x10um '10' ('9' is high freq), 2x8um '11', 2x4um '12' ('13' is high freq)), F1 ('1')
# tmiss_to50K = tmiss6_interp*phd10_tmiss_interp*phd8_tmiss_interp*phd8_tmiss_interp*phd4_tmiss_interp*phd4_tmiss_interp*tmiss1_interp
# # to 4 K: add F2 ('5')
# tmiss_to4K = tmiss_to50K*tmiss5_interp
# # to 1 K: add F3 ('2')
# tmiss_to1K = tmiss_to4K*tmiss2_interp
# # to 300 mK: add F4, one for 90-180 GHz ('4') and another for 180-360 GHz ('3' has smallest freq range but highest fidelity)
# tmiss_to300mK_180GHz = tmiss_to1K*tmiss4_interp
# tmiss_to300mK_360GHz = tmiss_to1K*tmiss3_interp

# to 50 K: AR Window ('5'), DSIR Filters 1-5 (1x10um '9' ('8' is high freq), 2x8um '10', 2x4um '11' ('12' is high freq)), F1 ('0')
tmiss_to50K = AR_tmiss*DSIR1_tmiss*DSIR2_tmiss*DSIR2_tmiss*DSIR4_tmiss*DSIR4_tmiss*F1_tmiss
# to 4 K: add F2
tmiss_to4K = tmiss_to50K*F2_tmiss
# to 1 K: add F3 ('1')
tmiss_to1K = tmiss_to4K*F3_tmiss
# to 300 mK: add F4, one for 90-180 GHz ('3') and another for 180-360 GHz
tmiss_to300mK_180GHz = tmiss_to1K*F4_180GHz_tmiss
tmiss_to300mK_360GHz = tmiss_to1K*F4_360GHz_tmiss

# plt.figure(figsize=figsize)
fig, ax = plt.subplots(figsize=figsize, layout='tight')
plt.plot(freqs_interp, tmiss_to50K,          label='50 K',                 alpha=1, linewidth=2.5)
plt.plot(freqs_interp, tmiss_to4K,           label='4 K',                  alpha=1, linewidth=2.5)
plt.plot(freqs_interp, tmiss_to1K,           label='1 K',                  alpha=1, linewidth=2.5)
plt.plot(freqs_interp, tmiss_to300mK_180GHz, label='300 mK (90-180 GHz)',  alpha=1, linewidth=2.5)
plt.plot(freqs_interp, tmiss_to300mK_360GHz, label='300 mK (180-360 GHz)', alpha=1, linewidth=2.5)
plt.xlabel('Frequency [GHz]')
plt.legend(markerscale=5)

secax = ax.secondary_xaxis('top', functions=(lambda x: x/spectoGHz, lambda x: x*spectoGHz))
secax.set_xlabel('Wave Number [1/cm]')
plt.ylim(0,1); plt.xlim(0, 350)   # GHz
plt.ylabel('Total Transmission')
plt.grid(linestyle = '--', which='both', linewidth=0.5)
plt.tight_layout()
plt.show()