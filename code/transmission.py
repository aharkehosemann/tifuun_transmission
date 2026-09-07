# Angi Harke-Hosemann 2026/06
#
# calculates total transmission of TIFUUN filters from Cardiff's filter measurements
#
# dependencies: numpy, matplotlib, scipy, csv, collections
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
# 9  = DSIR5 (phd10) 50 K, high freq measurement used for freqs ~> 225 1/cm = 6750 GHz (irrelevant to band of interest)
# 10 = DSIR5 (phd10) 50 K, low freq measurement used for freqs ~< 225 1/cm = 6750 GHz (only measurement relevant to band of interest)
# 11 = DSIR3&4 (phd8) 140&50 K
# 12 = DSIR1&2 (phd4) 300&140 K, low freq measurement used for freqs ~< 310 1/cm = 9300 GHz (still irrelevant to band of interest)
# 13 = DSIR1&2 (phd4) 300&140 K, high freq measurement used for freqs ~> 310 1/cm = 9300 GHz (irrelevant to band of interest)
#
# notes
# TIFUUN bands are 130–178 GHz and 195–319 GHz, or 90 - 360 GHz according to Akira's SPIE proceedings
# measurement bands: 1-4 = 17.0-330 GHz, 5 = 24-1049 GHz, 6 = 24-1079 GHz, 7 = 25-1199 GHz, 8 = 25-1129 GHz, 9 = 3900-180000GHz, 10-11 = 600-10500 GHz, 12 = 120-10500 GHz, 13 = 150-10500 GHz
# duplicate values removed: 216/7007 = 3% for 1, 216/7007 = 3% for 2, 194/4004 = 4.9% for 3, 156/7007 = 2.2% for 4, 0/1001 for 5, 0/1001 for 6-8, 1/12490 = 0.01% for 12&13, 1/26090 <<1% for 9/10, 1/15093 <<1% for 11

from transmission_routines import *

### user configuration and analysis options
root_dir                = '/Users/angi/tifuun/transmission/'
# msmts_to_plot           = [12, 13, 11, 10]   # set to [] for no raw measurement plots; 1=F1, 2=F3, 3=F4 360 GHz, 4=F4 180 GHz, 5=F2, 6=AR window, 7=F4 360 GHz, 8=F4 360 GHz, 9=DSIR5 HF (phd10), 10=DSIR5 LF (phd10), 11=DSIR3&4 (phd8), 12=DSIR1&2 LF (phd4), 13=DSIR1&2 HF (phd4)
msmts_to_plot           = []   # set to [] for no raw measurement plots; 1=F1, 2=F3, 3=F4 360 GHz, 4=F4 180 GHz, 5=F2, 6=AR window, 7=F4 360 GHz, 8=F4 360 GHz, 9=DSIR5 HF (phd10), 10=DSIR5 LF (phd10), 11=DSIR3&4 (phd8), 12=DSIR1&2 LF (phd4), 13=DSIR1&2 HF (phd4)
ave_duplicates          = True   # average duplicate frequency measurements when interpolating to new frequencies, otherwise second instance of duplicates will be removed
check_interp            = False   # plot for checking interpolation of cleaned data to new frequencies
analyze_DSIR345         = True   # analyze DSIR 3&4 (phd8) and 5 (phd10) measurements in addition to DSIR 1&2 (phd4) - phd4 is the only measurement relevant to band of interest
check_phd4              = True   # plot for checking cleaning of DSIR 1&2 (phd4) cleaning and interpolation
check_phd8              = True   # plot for checking cleaning of DSIR 3&4 (phd8) cleaning and interpolation
check_phd10             = True   # plot for checking cleaning of DSIR 5 (phd10) cleaning and interpolation
plot_total_transmission = False   # plot total transmission to 50 K, 4 K, 1 K, and 300 mK stages
freq_min = 0; freq_max = 330; num_freqs = 1000   # GHz, frequency range to interpolate total transmission

### plot settings
font = {'family' : 'serif', 'weight' : 'normal', 'size'   : 18}
ticks = {'major.size' : '5', 'labelsize' : '14', 'minor.visible' : False}
grids = {'linestyle' : '--', 'linewidth' : 0.5}; axes = { 'grid' : True }
plt.rc('text', usetex=True); plt.rc('font', **font)
plt.rc('xtick', **ticks);    plt.rc('ytick', **ticks)
plt.rc('grid', **grids);     plt.rc('axes', **axes)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}'
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['figure.dpi'] = 100; plt.rcParams['savefig.dpi'] = 300   # higher resolution plots in interactive notebook and when saving pngs
figsize = (10, 6)   # width, height

### filter measurement files
filter_file1  = root_dir+'tmiss_measurements/3668data.csv'      # FP 3668 ARC, LPE F1? - 50 K thick IR blocker 90-360 GHz
filter_file2  = root_dir+'tmiss_measurements/3667data.csv'      # FP 3667,     LPE F3? - 1K filter ~360 GHz
filter_file3  = root_dir+'tmiss_measurements/3662data.csv'      # FP 3662,     LPE F4 360 GHz? - 300 mK filter 180-360 GHz
filter_file4  = root_dir+'tmiss_measurements/3572data.csv'      # FP 3572,     LPE F4 180 GHz? - 300 mK filter 990-180 GHz
filter_file5  = root_dir+'tmiss_measurements/S3431R5.csv'       # FP 3450,     LPE F2? - 4K LPE filter 500 GHz
filter_file6  = root_dir+'tmiss_measurements/S3424R7.csv'       # FP 3429 ARC, UHMWPE ARC window at 300 K?
filter_file7  = root_dir+'tmiss_measurements/S3400R11.csv'      # FP 3083, also LPE F4 360 GHz? - 300 mK filter 180-360 GHz
filter_file8  = root_dir+'tmiss_measurements/S3400R13.csv'      # FP 3083, also LPE F4 360 GHz? - 300 mK filter 180-360 GHz
filter_file9  = root_dir+'tmiss_measurements/C0257_3.csv'       # DSIR5 (phd10) 50 K, one measurement, use for 225+ 1/cm
filter_file10 = root_dir+'tmiss_measurements/T1889R7.csv'       # DSIR5 (phd10) 50 K, another measurement, use for 0-225 1/cm
filter_file11 = root_dir+'tmiss_measurements/T1889R5.csv'       # DSIR3&4 (phd8) 140&50 K
filter_file12 = root_dir+'tmiss_measurements/phd4_combined.csv' # DSIR1&2 (phd4) 300&140 K, one measurement, Carole combined low and mid-frequency measurements
filter_file13 = root_dir+'tmiss_measurements/M24779.csv'        # DSIR1&2 (phd4) 300&140 K, another measurement, Carole suggests switching to this around 330 1/cm
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

tmiss_all, serial_all, filter_all = read_transmission_csv(filter_files, filter_cols)
serials    = list(tmiss_all)

# DSIR filters need extra processing
freq9_raw  = tmiss_all[serials[8]]['freqs_raw']; freq10_raw  = tmiss_all[serials[9]]['freqs_raw']; freq11_raw  = tmiss_all[serials[10]]['freqs_raw']; freq12_raw  = tmiss_all[serials[11]]['freqs_raw']; freq13_raw  = tmiss_all[serials[12]]['freqs_raw']
tmiss9_raw = tmiss_all[serials[8]]['tmiss_raw']; tmiss10_raw = tmiss_all[serials[9]]['tmiss_raw']; tmiss11_raw = tmiss_all[serials[10]]['tmiss_raw']; tmiss12_raw = tmiss_all[serials[11]]['tmiss_raw']; tmiss13_raw = tmiss_all[serials[12]]['tmiss_raw']

### plot raw filter measurements
if len(msmts_to_plot)>0:
  fig, ax = plt.subplots(figsize=figsize, layout='tight')
  for mm in msmts_to_plot:
    plt.plot(tmiss_all[serials[mm-1]]['freqs_raw'], tmiss_all[serials[mm-1]]['tmiss_raw'], '.', markersize=5, label=serial_to_filter[serials[mm-1]], alpha=0.7)
  plt.xlabel('Frequency [GHz]')   # GHz
  plt.ylabel('Filter Transmission')
  plt.legend(markerscale=2, loc='lower left'); plt.ylim(0,1); plt.xlim(0,700)
  # plt.legend(markerscale=2, loc='lower left'); plt.ylim(0.8,1); plt.xlim(0,1000)
  secax = ax.secondary_xaxis('top', functions=(lambda x: x/spectoGHz, lambda x: x*spectoGHz)); secax.set_xlabel('Wave Number [1/cm]')

### process data
# average or remove duplicate frequency measurements and interpolate to new frequencies
freqs_interp = np.linspace(freq_min, freq_max, num=num_freqs)   # GHz, frequencies to interpolate total transmission

for ff in np.arange(8):   # for each filter measurement, clean and interpolate at new frequencies
  freq_sorted, tmiss_sorted, freq_interp, tmiss_interp = clean_and_interp(tmiss_all[serials[ff]]['freqs_raw'], tmiss_all[serials[ff]]['tmiss_raw'], freqs_interp, average=ave_duplicates, check_interp=check_interp)
  tmiss_all[serials[ff]]['freqs_clean'], tmiss_all[serials[ff]]['tmiss_clean'], tmiss_all[serials[ff]]['freqs_interp'], tmiss_all[serials[ff]]['tmiss_interp'] = freq_sorted, tmiss_sorted, freq_interp, tmiss_interp

## DSIR measurements
# first low and high freq msmts are combined if present, then phd4 msmts are tacked onto phd8 and phd10 to extend their freq range from 600 to 120 GHz, and finally transmission is set to unity for all DSIR filters below 120 GHz
# for frequency range of interest, all filters are characterized using phd4 measurements, so handling of phd8 and phd10 data is unnecessary
# combine low and high freq DSIR 4 um measurements
phd4_hfinds = np.where(freq13_raw>=310*spectoGHz)[0]   # high freq = large wavenumber
phd4_lfinds = np.where(freq12_raw<310*spectoGHz)[0]    # low freq = small wavenumber
phd4_freq   = np.concatenate((freq12_raw[phd4_lfinds],  freq13_raw[phd4_hfinds]))
phd4_tmiss  = np.concatenate((tmiss12_raw[phd4_lfinds], tmiss13_raw[phd4_hfinds]))

# set transmission of DSIR filters to unity below 120 GHz
unitymask_freqs = np.linspace(0, min(phd4_freq)); unitymask_freqs = unitymask_freqs[:-1]   # don't overwrite first data point
unitymask_tmiss = np.ones_like(unitymask_freqs)
phd4_freq_raw   = np.concatenate((unitymask_freqs, phd4_freq)); phd4_tmiss_raw = np.concatenate((unitymask_tmiss, phd4_tmiss))

# handle duplicate frequency measurments and interpolate
phd4_freq_cleaned, phd4_tmiss_cleaned, phd4_freq_interp, phd4_tmiss_interp = clean_and_interp(phd4_freq_raw, phd4_tmiss_raw, freqs_interp, check_interp=check_interp)   # removed 1/12490 = <<1% duplicate x values, potentially in range but hard to tell because measurements go to 150 THz

if analyze_DSIR345:
  # combine low and high freq DSIR 10 um measurements
  phd10_hfinds = np.where(freq9_raw>=225*spectoGHz)[0]   # high freq = large wavenumber
  phd10_lfinds = np.where(freq10_raw<225*spectoGHz)[0]   # low freq = small wavenumber
  phd10c_freq  = np.concatenate((freq9_raw[phd10_hfinds],  freq10_raw[phd10_lfinds]))
  phd10c_tmiss = np.concatenate((tmiss9_raw[phd10_hfinds], tmiss10_raw[phd10_lfinds]))

  # tack on DSIR 4 um measurements below 600 GHz to 8 um and 10 um measurements
  phd4_maskinds = np.where(phd4_tmiss<min(np.concatenate((freq11_raw, phd10c_freq))))[0]
  phd10_freq    = np.concatenate((phd4_freq[phd4_maskinds],  phd10c_freq))
  phd10_tmiss   = np.concatenate((phd4_tmiss[phd4_maskinds], phd10c_tmiss))
  phd8_freq     = np.concatenate((phd4_freq[phd4_maskinds],  freq11_raw))
  phd8_tmiss    = np.concatenate((phd4_tmiss[phd4_maskinds], tmiss11_raw))

  # set transmission of DSIR filters to unity below 120 GHz
  phd10_freq_raw  = np.concatenate((unitymask_freqs, phd10_freq)); phd10_tmiss_raw = np.concatenate((unitymask_tmiss, phd10_tmiss))
  phd8_freq_raw   = np.concatenate((unitymask_freqs, phd8_freq));  phd8_tmiss_raw  = np.concatenate((unitymask_tmiss, phd8_tmiss))

  # handle duplicate frequency measurments and interpolate
  phd10_freq_cleaned, phd10_tmiss_cleaned, phd10_freq_interp, phd10_tmiss_interp = clean_and_interp(phd10_freq_raw,  phd10_tmiss_raw, freqs_interp, check_interp=check_interp)   # removed 1/26090 = <<1% duplicate x values, probably inherited from phd4
  phd8_freq_cleaned,  phd8_tmiss_cleaned,  phd8_freq_interp,  phd8_tmiss_interp  = clean_and_interp(phd8_freq_raw,   phd8_tmiss_raw,  freqs_interp, check_interp=check_interp)   # removed 1/15093 = <<1% duplicate x values, probably inherited from phd4
else: # treat phd8 and phd10 filters as phd4
  phd10_freq_raw, phd10_tmiss_raw, phd10_freq_cleaned, phd10_tmiss_cleaned, phd10_freq_interp, phd10_tmiss_interp = phd4_freq_raw, phd4_tmiss_raw, phd4_freq_cleaned, phd4_tmiss_cleaned, phd4_freq_interp, phd4_tmiss_interp
  phd8_freq_raw, phd8_tmiss_raw, phd8_freq_cleaned, phd8_tmiss_cleaned, phd8_freq_interp, phd8_tmiss_interp = phd4_freq_raw, phd4_tmiss_raw, phd4_freq_cleaned, phd4_tmiss_cleaned, phd4_freq_interp, phd4_tmiss_interp

if check_phd4:   # check cleaning and interpolation of DSIR 1&2 (phd4)
  fig, ax = plt.subplots(figsize=(10,7), layout='tight')
  plt.plot(freq12_raw,        tmiss12_raw,        'o', markersize=5, alpha=1, color='C1', label='raw, LF')    # DSIR1&2 (phd4) 300&140 K, one measurement
  plt.plot(freq13_raw,        tmiss13_raw,        'o', markersize=5, alpha=1, color='C3', label='raw, HF')    # DSIR1&2 (phd4) 300&140 K, another measurement
  # plt.plot(phd4_freq_raw,     phd4_tmiss_raw,     'o', markersize=3, alpha=0.6, color='C2', label='pre-cleaned')    # DSIR1&2 (phd4) 300&140 K, pre-cleaned
  plt.plot(phd4_freq_cleaned, phd4_tmiss_cleaned, 'o', markersize=3, alpha=0.6, color='C0', label='cleaned')    # DSIR1&2 (phd4) 300&140 K, cleaned
  plt.plot(phd4_freq_interp,  phd4_tmiss_interp,  linewidth=2, color='k', label='interpolated')    # DSIR1&2 (phd4) 300&140 K, interpolated
  plt.legend(markerscale=2, loc='lower left')
  plt.xlim(0,700); plt.xlabel('Frequency [GHz]')   # GHz
  plt.ylim(0.88, 1.01); plt.ylabel('Filter Transmission')
  plt.title('DSIR 1\\&2 (PHD4)', pad=10)
  ax.fill_between(phd4_freq_cleaned, 0, 2, where=phd4_freq_cleaned>freq_max, color='gray', alpha=0.2, linewidth=0)   # shade region above max frequency of interest
  secax = ax.secondary_xaxis('top', functions=(lambda x: x/spectoGHz, lambda x: x*spectoGHz)); secax.set_xlabel('Wave Number [1/cm]', labelpad=10)

if check_phd8:   # check cleaning and interpolation of DSIR 3&4 (phd8)
  fig, ax = plt.subplots(figsize=(10,7), layout='tight')
  plt.plot(freq11_raw,        tmiss11_raw,        'o', markersize=5, alpha=1, color='C1', label='raw')    # DSIR3&4 (phd8) 140&50 K, single measurement
  # plt.plot(phd8_freq_raw,     phd8_tmiss_raw,     'o', markersize=3, alpha=0.6, color='C2', label='pre-cleaned')    # DSIR3&4 (phd8) 140&50 K, pre-cleaned
  plt.plot(phd8_freq_cleaned, phd8_tmiss_cleaned, 'o', markersize=3, alpha=0.6, color='C0', label='cleaned')    # DSIR3&4 (phd8) 140&50 K, cleaned
  plt.plot(phd8_freq_interp,  phd8_tmiss_interp,  linewidth=2, color='k', label='interpolated')    # DSIR3&4 (phd8) 140&50 K, interpolated
  plt.legend(markerscale=2, loc='lower left')
  plt.xlim(0,700); plt.xlabel('Frequency [GHz]')   # GHz
  plt.ylim(0.88, 1.01); plt.ylabel('Filter Transmission')
  plt.title('DSIR 3\\&4 (PHD8)', pad=10)
  ax.fill_between(phd8_freq_cleaned, 0, 2, where=phd8_freq_cleaned>freq_max, color='gray', alpha=0.2, linewidth=0)   # shade region above max frequency of interest
  secax = ax.secondary_xaxis('top', functions=(lambda x: x/spectoGHz, lambda x: x*spectoGHz)); secax.set_xlabel('Wave Number [1/cm]', labelpad=10)

if check_phd10:   # check cleaning and interpolation of DSIR 5 (phd10)
  fig, ax = plt.subplots(figsize=(10,7), layout='tight')
  plt.plot(freq10_raw,         tmiss10_raw,         'o', markersize=5, alpha=1, color='C1', label='raw, LF')    # DSIR5 (phd10) 50 K, one measurement
  plt.plot(freq9_raw,          tmiss9_raw,          'o', markersize=5, alpha=1, color='C3', label='raw, HF')    # DSIR5 (phd10) 50 K, one measurement
  # plt.plot(phd10_freq_raw,     phd10_tmiss_raw,     'o', markersize=3, alpha=0.6, color='C2', label='pre-cleaned')    # DSIR5 (phd10) 50 K, pre-cleaned
  plt.plot(phd10_freq_cleaned, phd10_tmiss_cleaned, 'o', markersize=3, alpha=0.6, color='C0', label='cleaned')    # DSIR5 (phd10) 50 K, cleaned
  plt.plot(phd10_freq_interp,  phd10_tmiss_interp,  linewidth=2, color='k', label='interpolated')    # DSIR5 (phd10) 50 K, interpolated
  plt.legend(markerscale=2, loc='lower left')
  plt.xlim(0,700); plt.xlabel('Frequency [GHz]')   # GHz
  plt.ylim(0.88, 1.01); plt.ylabel('Filter Transmission')
  plt.title('DSIR 5 (PHD10)', pad=10)
  ax.fill_between(phd10_freq_cleaned, 0, 2, where=phd10_freq_cleaned>freq_max, color='gray', alpha=0.2, linewidth=0)   # shade region above max frequency of interest
  secax = ax.secondary_xaxis('top', functions=(lambda x: x/spectoGHz, lambda x: x*spectoGHz)); secax.set_xlabel('Wave Number [1/cm]', labelpad=10)

### total transmission to 50 K, 4 K, 1 K, and 300 mK
if plot_total_transmission:
  AR_tmiss        = tmiss_all[serials[5]]['tmiss_interp']   # only measurement
  F1_tmiss        = tmiss_all[serials[0]]['tmiss_interp']   # only measurement
  F2_tmiss        = tmiss_all[serials[4]]['tmiss_interp']   # only measurement
  F3_tmiss        = tmiss_all[serials[1]]['tmiss_interp']   # only measurement
  F4_180GHz_tmiss = tmiss_all[serials[3]]['tmiss_interp']   # only measurement
  F4_360GHz_tmiss = tmiss_all[serials[2]]['tmiss_interp']   # 3 msmts, file '3' has smallest freq range but highest fidelity
  DSIR12_tmiss    = phd4_tmiss_interp
  DSIR34_tmiss    = phd8_tmiss_interp
  DSIR5_tmiss     = phd10_tmiss_interp

  # to 50 K: AR Window ('5'), DSIR Filters 1-5 (1x10um '10' ('9' is high freq), 2x8um '11', 2x4um '12' ('13' is high freq)), F1 ('1')
  tmiss_to50K = AR_tmiss*DSIR12_tmiss*DSIR12_tmiss*DSIR34_tmiss*DSIR34_tmiss*DSIR5_tmiss*F1_tmiss
  # to 4 K: add F2
  tmiss_to4K = tmiss_to50K*F2_tmiss
  # to 1 K: add F3 ('1')
  tmiss_to1K = tmiss_to4K*F3_tmiss
  # to 300 mK: add F4, one for 90-180 GHz ('3') and another for 180-360 GHz
  tmiss_to300mK_180GHz = tmiss_to1K*F4_180GHz_tmiss
  tmiss_to300mK_360GHz = tmiss_to1K*F4_360GHz_tmiss

  fig, ax = plt.subplots(figsize=figsize, layout='tight')
  plt.plot(freqs_interp, tmiss_to50K,          label='50 K',                 alpha=1, linewidth=2.5)
  plt.plot(freqs_interp, tmiss_to4K,           label='4 K',                  alpha=1, linewidth=2.5)
  plt.plot(freqs_interp, tmiss_to1K,           label='1 K',                  alpha=1, linewidth=2.5)
  plt.plot(freqs_interp, tmiss_to300mK_180GHz, label='300 mK (90-180 GHz)',  alpha=1, linewidth=2.5)
  plt.plot(freqs_interp, tmiss_to300mK_360GHz, label='300 mK (180-360 GHz)', alpha=1, linewidth=2.5)
  plt.legend(markerscale=5)
  plt.ylim(0,1); plt.xlim(0, 350)   # GHz
  plt.xlabel('Frequency [GHz]')
  plt.ylabel('Total Transmission')
  secax = ax.secondary_xaxis('top', functions=(lambda x: x/spectoGHz, lambda x: x*spectoGHz)); secax.set_xlabel('Wave Number [1/cm]')

plt.show()