#!/usr/bin/env python3
# coding: utf-8
# /************************************************************************/
# /* Package:  ARTI                                                       */
# /* Module:   extractor.py                                               */
# /************************************************************************/
# /* Authors:  Hernán Asorey                                              */
# /* e-mail:   hernan.asoreyh@iteda.cnea.gov.ar                           */
# /************************************************************************/
# /* Comments: Get GDAS instantaneous and monthly averages atmosheric     */
# /*           profiles for a site. Site definitions are stored in a      */
# /*           public read-only google spreadsheet. If you need a new     */
# /*           site please contact us.                                    */
# /************************************************************************/
# /*
# LICENSE BSD-3-Clause
# Copyright (c) 2015
# The LAGO Collaboration
# https://lagoproject.net
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#    1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#    2. Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
# NO EVENT SHALL LAB DPR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation are
# those of the authors and should not be interpreted as representing
# official policies, either expressed or implied, of the LAGO Collaboration.
# 
# */
# /************************************************************************/


# for now, using instantaneous profiles at local noon and mignight for each site
# ## Preliminaries
# required:
# sudo pip install pytz gspread numpy scipy pandas timezonefinder oauth2client 

from datetime import datetime
import sys
import getopt
import subprocess
import pytz
from timezonefinder import TimezoneFinder
from oauth2client.service_account import ServiceAccountCredentials
import gspread
from calendar import monthrange
from collections import namedtuple
import os.path
import numpy as np
import pandas as pd


def days_in_year(y, date):
    for month in range(1, 13):
        for day in range(1, monthrange(y, month)[1] + 1):
            yield date(y, month, day)


def gdas_get_time(site, gdas_date, gdas_hour, row):
    gdas_time = datetime(gdas_date.year, gdas_date.month, gdas_date.day, gdas_hour.hour, gdas_hour.min, gdas_hour.sec)
    tz = pytz.timezone(TimezoneFinder().timezone_at(lat=site["LAT"], lng=site["LONG"]))
    utc = str(int(datetime.timestamp(tz.localize(gdas_time).astimezone(pytz.UTC))))
    file = "atmg{}{}.atm".format('{:04d}'.format(row["SiteId"]), utc)
    return utc, file


def get_layer(an_h, an_atm):
    layer = 4
    for i in range(1, 5):
        if an_h < an_atm[0][i]:
            layer = i - 1
            break
    return layer


def depth(an_h, an_atm):
    layer = get_layer(an_h, an_atm)
    if layer == 4:
        t = an_atm[1][layer] - an_atm[2][layer] * an_h / an_atm[3][layer]
    else:
        t = an_atm[1][layer] + an_atm[2][layer] * np.exp(- an_h / an_atm[3][layer])
    return t


def density(an_h, an_atm):
    layer = get_layer(an_h, an_atm)
    if layer == 4:
        d = (an_atm[2][layer] / an_atm[3][layer]) * (an_atm[0][layer] / an_h)
    else:
        d = an_atm[2][layer] * np.exp(- an_h / an_atm[3][layer]) / an_atm[3][layer]
    return d


def n_index():
    # for now, we avoid this calculation as it is not needed
    return 3.e-3


def main(argv):
    # Some auxiliary variables
    gdas_atm_path = '../atm/'
    gdas_local_hour = [00, 12]
    gdas_min_height = str(0.)
    gdas_max_height = str(0.)

    # arrays to hold atmospheric variables, height, density, thick, refraction index
    # for height, we use the same altitudes as in CORSIKA external variables
    h = np.concatenate((np.arange(0, 25, 1), np.arange(25, 50, 2.5), np.arange(50, 125, 5)))
    rho = np.zeros_like(h)
    thick = np.zeros_like(h)
    n = np.zeros_like(h)
    # number of files that are being averaged
    n_files = 0
    # Default configuration
    verbose = False
    extract = True
    average = True
    gdas_site = False
    start_year = 2018
    end_year = start_year
    default_last = True
    # Command line arguments
    try:
        opts, args = getopt.getopt(argv, "hveas:y:d:",
                                   ["verbose", "extract", "average", "site=", "year=", "end="])
    except getopt.GetoptError:
        print('Error parsing arguments. Try extractor.py -h for help')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('extractor.py -h(elp) -v(erbose) -e(xtract) -a(verage) -s <site> -y <year> -d <last year>')
            print("Defaults:")
            print(f"    verbosity :  {verbose}")
            print(f"    extract   :  {extract}")
            print(f"    average   :  {average}")
            print(f"    site      :  {gdas_site}")  
            print(f"    year      :  {start_year}")
            print(f"    last year :  {end_year}")
            sys.exit()
        elif opt in ("-v", "--verbose"):
            verbose = not verbose
            print(f"    verbosity :  {verbose}")
        elif opt in ("-e", "--extract"):
            extract = not extract
            print(f"    extract   :  {extract}")
        elif opt in ("-a", "--average"):
            average = not average
            print(f"    average   :  {average}")
        elif opt in ("-s", "--site"):
            gdas_site = arg
            print(f"    site      :  {gdas_site}")
        elif opt in ("-y", "--year"):
            start_year = int(arg)
            print(f"    year      :  {start_year}")
        elif opt in ("-d", "--end"):
            end_year = arg
            default_last = False
            print(f"    last year :  {end_year}")

    if default_last: 
        end_year = start_year
        print(f"    last year :  {end_year}")

    # time and date
    Date = namedtuple("Date", ["year", "month", "day"])
    Time = namedtuple("Time", ["hour", "min", "sec"])

    # ## Sites definitions
    # Interacting with Sites Definition at
    # [Google Sheet]
    # (https://docs.google.com/spreadsheets/d/1lBCheS1KuygPQU_j_VDySrm23gTU8QuF9snXPtHm2mI/edit?usp=sharing)
    # trough Google API. 
    scope = [
        'https://spreadsheets.google.com/feeds',
        'https://www.googleapis.com/auth/drive'
    ]
    if verbose:
        print("Extracting data from sites google sheet...", end='')
    credentials = ServiceAccountCredentials.from_json_keyfile_name('access-data-64a03ecb59d0.json', scope)
    client = gspread.authorize(credentials)
    sheet = client.open("DatosRC").sheet1
    data = sheet.get_all_records()
    if verbose:
        print("done.")

    # ## GDAS atmospheric integration
    # We are using the `gdastool` included into CORSIKA regular installation since v7.46.
    # The first time it runs for a specific month it could take a couple of minutes as it has
    # to download the corresponding gdas file (250MB-500MB per GDAS file, only the first time)
    # Warning: GDAS files are updated on days 1, 7, 14, 21 and 28.
    for row in data:
        if not gdas_site or int(row['SiteId']) == int(gdas_site):
            for year in range(start_year, end_year + 1):
                for date in days_in_year(year, Date):
                    for hour in gdas_local_hour:
                        utctime, gdas_file = gdas_get_time(row, date, Time(hour, 0, 0), row)
                        gdas_file = gdas_atm_path + gdas_file
                        if verbose:
                            print(f"{row['SiteId']} {date} {hour} ", end="")
                        if extract:
                            if not os.path.isfile(gdas_file):
                                print("Extracting GDAS atmospheric profiles")
                                print("(it could take a while if the gdas file must be downloaded, ~500MB each)")
                                gdas_cmd = [
                                    'python3',
                                    'gdastool',
                                    '-t', utctime,
                                    '-o', gdas_file,
                                    '-m', gdas_min_height,
                                    '-M', gdas_max_height,
                                    '-c', str(row["LAT"]), str(row["LONG"])
                                ]
                                gdas = subprocess.run(gdas_cmd)
                                print(gdas)
                            elif verbose:
                                print("... file exists... ", end="")
                        if average:
                            if os.path.isfile(gdas_file):
                                atm = np.loadtxt(fname=gdas_file, max_rows=4, comments="#")
                                rho += np.array([density(hi * 1.0e5, atm) for hi in h])
                                thick += np.array([depth(hi * 1.0e5, atm) for hi in h])
                                n += np.array([n_index() for hi in h])
                                n_files += 1
                                if date.day == monthrange(date.year, date.month)[1] and hour == gdas_local_hour[-1]:
                                    if n_files > 0:
                                        # averaging
                                        rho /= n_files
                                        thick /= n_files
                                        n /= n_files
                                        # data frame
                                        avg_atm = pd.DataFrame(
                                            {"a": h, "r": rho, "t": thick, "n": n}
                                        )
                                        # taking care of negative values and last values/
                                        avg_atm[avg_atm < 0] = 1e-5
                                        avg_atm.at[49, 't'] = 0
                                        # output
                                        avg_id = "{}{}{}".format(
                                            '{:d}'.format(row["SiteId"]),
                                            '{:02d}'.format(date.year-int(date.year/100)*100),
                                            '{:02d}'.format(date.month)
                                        )
                                        avg_file = f"atmprof{avg_id}.dat"
                                        f = open(gdas_atm_path + avg_file, "w")
                                        f.write(f"# Atmospheric Model {avg_id}\n")
                                        f.write("# Col. #1          #2           #3            #4\n")
                                        f.write("# Alt [km]    rho [g/cm^3] thick [g/cm^2]    n-1\n")
                                        f.close()
                                        avg_atm.to_csv(gdas_atm_path + avg_file,
                                                       float_format='%.5E', header=None, index=None, sep=' ', mode='a')
                                        print(f"new {avg_file} created")
                                        rho = np.zeros_like(h)
                                        thick = np.zeros_like(h)
                                        n = np.zeros_like(h)
                                        n_files = 0
                                    else:
                                        print(f"Warning: no files for {date.month}/{date.year}")
                            else:
                                print(f"No file {gdas_file}")
                        if verbose:
                            print("done")


if __name__ == "__main__":
    main(sys.argv[1:])
