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
        self.date_list = list()
        self.__list_files()
        self.__set_to_last(date)

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

    def __list_files(self):
        '''
        we get the list the date of the different versions of the file
        '''
        date_list = glob.glob(
            self.get_file_path() +
            "/*" +
            self.get_file_name())
        if len(date_list) > 0:
            format_search = re.compile(
                r".*" +
                re.escape(str(self.file_name)))
            for i in range(len(date_list)):
                if format_search.match(date_list[i]):
                    date = os.path.basename(date_list[i])
                    if self.__test_date(date):
                        date = self.__extract_date(date, True)
                        self.date_list.append(date)
            self.date_list.sort(reverse=True)

    def __set_to_last(self, date):
        '''
        if the file exist in different version and we don't have a date
        set file to the last version
        '''
        test = date is None and \
            self.date == datetime.date.today()
        if test and len(self.date_list) > 0:
            self.date = self.date_list[0]

    def get_file_name(self):
        '''
        return the file name without date
        '''
        return(str(self.file_name))

    def get_file_path(self):
        '''
        return absolute path toward the file
        '''
        return(str(self.file_path))

    def __getitem__(self, key):
        '''
        return the full name of the file of indice key in the list of files
        ordering from the newest to the oldest file
        '''
        if key >= 0 and key < len(self.date_list):
            return(
                str(self.date_list[key].strftime("%Y_%m_%d_")) +
                str(self.file_name))
        else:
            return(None)

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
