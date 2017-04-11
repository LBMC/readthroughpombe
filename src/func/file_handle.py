#!/usr/bin/python3
# -*-coding:Utf-8 -*

import sys
import os.path
import argparse
import datetime
import re
import glob
if sys.version_info[0] == 2:
    print("file_handle.py is only compatible with python3.\
    Please run file_handle.py as an executable or with the command\
    'python3 file_handle.py'")
    exit(1)


class Dated_file:
    '''
    Dated_file class to manage file prefixed with a date in the format \
    yyyy_mm_dd_filemames.
    '''
    def __init__(self, file_name, date=None):
        file_name = os.path.abspath(file_name)
        self.file_name = os.path.basename(file_name)
        self.file_path = os.path.abspath(os.path.dirname(file_name))
        self.date = datetime.date(1, 1, 1)
        self.set_date(date)
        self.__truncate_file_name()

    def __test_date(self, date):
        '''
        test if the file_name contain a date tag
        '''
        format_test = re.match(
            r"\d{4}\_\d{2}\_\d{2}.*",
            date)
        return(format_test is not None)

    def __extract_date(self, date, return_date=False):
        '''
        extract date from a str beginning with yyyy_mm_dd
        '''
        format_search = re.search(
                r"(\d{4})\_(\d{2})\_(\d{2}).*",
                date)
        date = datetime.date(
            int(format_search.group(1)),
            int(format_search.group(2)),
            int(format_search.group(3))
        )
        if return_date:
            return date
        self.date = date

    def __truncate_file_name(self):
        '''
        remove date tag from file name
        '''
        format_search = re.search(
                r"\d{4}\_\d{2}\_\d{2}\_(.*)",
                self.file_name)
        if format_search is not None:
            self.file_name = format_search.group(1)

    def get_date(self):
        '''
        output date of the current file in yyyy_mm_dd format
        '''
        return(str(self.date.strftime("%Y_%m_%d")))

    def set_date(self, date):
        '''
        test if the date is in the yyyy_mm_dd format and is a real date
        if true, set date to current date, else exctact it
        '''
        if date is None:
            date = self.file_name
        if not self.__test_date(date):
            self.date = datetime.date.today()
        else:
            self.__extract_date(date)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog="file_handle.py",
        description="script to handle date in file name in the format \
        'yyyy_mm_dd_filemames'. By default return the last file corresponding \
        to the filename if existing or create it with the current date.")
    parser.add_argument(
        "-f",
        help="input filename",
        default=None,
        action="store",
        dest="input_file",
        required=True,
        nargs=1)
