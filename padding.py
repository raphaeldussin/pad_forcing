#!/usr/bin/env python

import argparse
import warnings
from glob import glob
import re
import xarray as xr
import os
import cftime
import numpy as np


def find_file(possible_directories, varname, yearstart, yearend, exclude_pattern=None):
    """find the file for a give varname spanning yearstart to yearend in possible_directories

    Args:
        possible_directories (list of str): directories containing forcing files
        varname (str): variable to look for in files
        yearstart (int): first year of the file
        yearend (int): last year of the file
        exclude_pattern (list of str, optional): patterns in filenames to exclude. Defaults to None.

    Returns:
        str: filename
    """

    all_files = []
    possible_files = []
    # get all files for the time segment
    for d in possible_directories:
        all_files += glob(f"{d}/*{yearstart}*{yearend}*")
        # if filename only have one year mentioned
        if len(all_files) == 0:
            all_files += glob(f"{d}/*{yearstart}*")
    # we should have at least one
    if len(all_files) == 0:
        exit(">>> ERROR: no files found to pad this file")
    # now look for variable name
    for f in all_files:
        out = re.findall(f"{varname}[(_|.)]", f)
        if len(out) > 0:
            possible_files.append(f)
    # we should have found one and only one
    if len(possible_files) > 1:
        if exclude_pattern is not None:
            for tag in exclude_pattern:
                for f in possible_files:
                    if f.find(tag) != -1:
                        possible_files.remove(f)
            if len(possible_files) != 1:
                exit(
                    f">>> ERROR: cannot converge on a single file, current left with {len(possible_files)}. Consider changing exclude pattern"
                )
        else:
            print(">>> found multiple files for padding", possible_files)
            print(">>> consider excluding the undesired pattern(s)")
            exit(">>> could not find single valid file")

    return possible_files[0]


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d",
        "--dirs",
        type=str,
        nargs="+",
        default="./",
        help="Path(s) to input nc file(s)",
    )

    parser.add_argument(
        "infile", metavar="INFILE", type=str, help="forcing file to pad"
    )

    parser.add_argument(
        "-f",
        "--first",
        action="store_true",
        help="this file is the first of the timeserie",
    )

    parser.add_argument(
        "-j",
        "--jra",
        action="store_true",
        help="this file is in the JRA forcing",
    )

    parser.add_argument(
        "-l",
        "--last",
        action="store_true",
        help="this file is the last of the timeserie",
    )

    parser.add_argument(
        "-v", "--var", required=True, type=str, help="variable to process"
    )

    parser.add_argument("-y", "--year", type=str, help="override for year to process")

    parser.add_argument(
        "-x", "--exclude", type=str, nargs="+", help="exclude some file patterns"
    )

    parser.add_argument(
        "--time",
        type=str,
        default="time",
        help="override for time variable name in file",
    )

    args = vars(parser.parse_args())

    if args["first"] and args["last"]:
        warnings.warn(
            "first and last are allowed simultaneouly but make sure this is only in the case of having a single file"
        )

    if args["year"] is None:
        # find the year from the filename
        # use regular expressions to get dates (have more than 4 digits)
        dates_from_fname = re.findall("[0-9][0-9][0-9][0-9]+", args["infile"])
        if len(dates_from_fname) == 2:
            start_date, end_date = dates_from_fname
        elif len(dates_from_fname) == 1:
            start_date = dates_from_fname[0]
            end_date = dates_from_fname[0]
        else:
            raise ValueError(
                ">>> ERROR: year retrieval yield not enough or too many values, use --year to override"
            )

        start_year = int(start_date[:4])
        end_year = int(end_date[:4])
    else:
        # this could be exapnded to handle different start and end year, if needed.
        start_year = args["year"]
        end_year = args["year"]

    print(f">>> current file spans the years {start_year} to {end_year}")

    years_per_file = end_year - start_year + 1
    previous_start_year = start_year - years_per_file
    previous_end_year = start_year - 1

    next_start_year = end_year + 1
    next_end_year = end_year + years_per_file

    varname = args["var"]
    timevar = args["time"]

    # find the directory containing the input file
    cdir = None
    for d in args["dirs"]:
        if os.path.exists(f"{d}/{args['infile']}"):
            cdir = d

    if cdir is None:
        exit(">>> ERROR: forcing file not found in directories")

    # open the file for current year and estimate frequency of data
    current = xr.open_dataset(f"{cdir}/{args['infile']}", decode_times=False)
    dt = current[timevar].diff(dim=timevar).isel({timevar: 0})

    start_at_midnight = float(current[timevar].isel({timevar:0}).data).is_integer()

    # inspect variables
    all_vars = list(set(list(current.variables) + list(current.coords)))
    possible_time_bounds = [timevar + "_bnds", timevar + "_bounds"]
    time_bounds = None
    for var in possible_time_bounds:
        if var in all_vars:
            time_bounds = var

    print(f"found time bounds = {time_bounds}")

    # pad at beginning of file with last record of previous file or doubling first record
    if args["first"]:
        # repeat the first record
        previous = current.isel({timevar: 0})

        if args["jra"]:
            if start_at_midnight:
                previous = None  # no need for a previous record
            else:
                # round down to midnight
                previous[timevar] = np.floor(previous[timevar])
                # not touching time bounds
        else:
            previous[timevar] = previous[timevar] - dt
            previous[timevar].attrs.update(current[timevar].attrs)
            if time_bounds is not None:
                previous[time_bounds] = previous[time_bounds] - dt
                previous[time_bounds].attrs.update(current[time_bounds].attrs)
    else:
        print(
            f">>> looking for a file for {varname} between {previous_start_year} and {previous_end_year}"
        )
        previous_file = find_file(
            args["dirs"],
            varname,
            previous_start_year,
            previous_end_year,
            exclude_pattern=args["exclude"],
        )
        print(f">>> found {previous_file}")
        previous = xr.open_dataset(previous_file, decode_times=False).isel(
            {timevar: -1}
        )

    # pad at end of file with first record of last file or doubling last record
    if args["last"]:
        # repeat the last record
        nexttime = current.isel({timevar: -1})
        nexttime[timevar] = nexttime[timevar] + dt
        nexttime[timevar].attrs.update(current[timevar].attrs)
        if time_bounds is not None:
            nexttime[time_bounds] = nexttime[time_bounds] + dt
            nexttime[time_bounds].attrs.update(current[time_bounds].attrs)
    else:
        print(
            f">>> looking for a file for {varname} between {next_start_year} and {next_end_year}"
        )
        next_file = find_file(
            args["dirs"],
            varname,
            next_start_year,
            next_end_year,
            exclude_pattern=args["exclude"],
        )
        print(f">>> found {next_file}")
        next_data = xr.open_dataset(next_file, decode_times=False).isel({timevar: 0})

    # concatenate records and force time to be first dimension
    if previous is not None:
        out = xr.concat([previous, current, next_data], dim=timevar)
    else:
        out = xr.concat([current, next_data], dim=timevar)

    out = out.transpose(*(timevar, ...))

    # override lon/lat bnds
    out["lon_bnds"] = current["lon_bnds"].copy(deep=True)
    out["lon_bnds"].encoding.update({"_FillValue": None, "coordinates": None})
    out["lat_bnds"] = current["lat_bnds"].copy(deep=True)
    out["lat_bnds"].encoding.update({"_FillValue": None, "coordinates": None})

    out["lon"] = current["lon"].copy(deep=True)
    out["lon"].encoding.update({"_FillValue": None, "coordinates": None})
    out["lat"] = current["lat"].copy(deep=True)
    out["lat"].encoding.update({"_FillValue": None, "coordinates": None})

    if "height" in list(out.variables):
        out["height"].encoding.update({"_FillValue": None})

    # clean up some dimensions/attributes
    for var in all_vars:
        if "comment" in out[var].attrs:
            if len(out[var].attrs["comment"]) > 500:
                out[var].attrs.pop("comment")

    # remove fill value on time
    out[timevar].encoding.update({"_FillValue": None, "chunksizes": (1,)})

    # remove fill value on time bounds, if appropriate
    if time_bounds is not None:
        out[time_bounds].encoding.update({"_FillValue": None, "coordinates": None})

    # write output file
    fileout = args["infile"].replace(".nc", ".padded.nc")
    out.to_netcdf(fileout, format="NETCDF4_CLASSIC", unlimited_dims=args["time"])
