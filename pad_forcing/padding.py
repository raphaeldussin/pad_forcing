import argparse
import warnings
from glob import glob
import re


def find_file(possible_directories, varname, yearstart, yearend):
    """ find the file for variable and time segment

    Args:
        possible_directories ([type]): [description]
        varname ([type]): [description]
        yearstart ([type]): [description]
        yearend ([type]): [description]
    """

    all_files = []
    possible_files = []
    # get all files for the time segment
    for d in possible_directories:
        all_files += glob(f"{d}/*{yearstart}*{yearend}*")
    # we should have at least one
    if len(all_files) == 0:
        raise ValueError("no files found to pad this file")
    # now look for variable name
    print(all_files)
    for f in all_files:
        out = re.findall(f"{varname}[(_|.)]", f)
        if len(out) > 0:
            possible_files.append(f)
    # we should have found one and only one
    if len(possible_files) > 1:
        print(possible_files)
        raise ValueError("could not find single valid file")

    return possible_files[0]


#def find_similar_file(possible_directories, fname, current_start, current_end, new_start, new_end):
#    expected_name = fname.replace(current_start, new_start)




parser = argparse.ArgumentParser()

parser.add_argument(
    "-d",
    "--dirs",
    type=str,
    nargs="+",
    default="./",
    help="Path(s) to input nc file(s)",
)

parser.add_argument("infile", metavar="INFILE", type=str, help="forcing file to pad")

parser.add_argument(
    "--first", action="store_true", help="this file is the first of the timeserie"
)
parser.add_argument(
    "--last", action="store_true", help="this file is the last of the timeserie"
)
parser.add_argument("--var", type=str, help="variable to process")
parser.add_argument("--year", type=str, help="override for year to process")

args = vars(parser.parse_args())

print(args)

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
            "year retrieval yield not enough or too many values, use --year to override"
        )

    start_year = int(start_date[:4])
    end_year = int(end_date[:4])
else:
    # this could be exapnded to handle different start and end year, if needed.
    start_year = args["year"]
    end_year = args["year"]

print(start_year, end_year)

years_per_file = end_year - start_year + 1
previous_start_year = start_year - years_per_file
previous_end_year = start_year - 1

next_start_year = end_year + 1
next_end_year = end_year + years_per_file

varname = args["var"]

if not args["first"]:
    print(f">>> looking for a file for {varname} between {previous_start_year} and {previous_end_year}")
    previous_file = find_file(args["dirs"], varname, previous_start_year, previous_end_year)
    print(f">>> found {previous_file}")

if not args["last"]:
    print(f">>> looking for a file for {varname} between {next_start_year} and {next_end_year}")
    next_file = find_file(args["dirs"], varname, next_start_year, next_end_year)
    print(f">>> found {next_file}")

