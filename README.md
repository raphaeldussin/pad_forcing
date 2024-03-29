# pad_forcing: pre-processing of forcing files for MOM6 (FMS)

MOM6 (FMS) requires forcing files to be "padded" (i.e. have the last record of the previous file and first record of next file)
in order to properly interpolate in time. This command line tool allows to pad forcing files located in multiple directories, exclude
filenames containing specific patterns (useful for pre-versions or incomplete years), has flags for special treatment of the first or 
last file of the timeserie (i.e. double first or last record).

```
usage: padding.py [-h] [-d DIRS [DIRS ...]] [-f] [-l] -v VAR [-y YEAR] [-x EXCLUDE [EXCLUDE ...]] [--time TIME] INFILE

positional arguments:
  INFILE                forcing file to pad

optional arguments:
  -h, --help            show this help message and exit
  -d DIRS [DIRS ...], --dirs DIRS [DIRS ...]
                        Path(s) to input nc file(s)
  -f, --first           this file is the first of the timeserie
  -l, --last            this file is the last of the timeserie
  -v VAR, --var VAR     variable to process
  -y YEAR, --year YEAR  override for year to process
  -x EXCLUDE [EXCLUDE ...], --exclude EXCLUDE [EXCLUDE ...]
                        exclude some file patterns
  --time TIME           override for time variable name in file
```

* Example:

```
./padding.py -d JRA55-do/v1.5.0/original JRA55do/v1.5.0.1/originals --var vas vas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_201901010000-201912312100.nc --exclude 202001010000-202007152100
>>> current file spans the years 2019 to 2019
>>> looking for a file for vas between 2018 and 2018
>>> found JRA55-do/v1.5.0/original/vas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_201801010000-201812312100.nc
>>> looking for a file for vas between 2020 and 2020
>>> found JRA55do/v1.5.0.1/originals/vas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0-1_gr_202001010000-202012312100.nc
```
