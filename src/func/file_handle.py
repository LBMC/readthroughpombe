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
        file_name = os.path.abspath(str(file_name))
        self.file_name = os.path.basename(str(file_name))
        self.file_path = os.path.dirname(str(file_name))
        self.date = datetime.date(1, 1, 1)
        self.set_date(date)
        self.__truncate_file_name()
        self.date_list = list()
        self.__list_files()
        self.__set_to_last(date)
        self.__date_file()

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
        if not self.__test_date(str(date)):
            self.date = datetime.date.today()
        else:
            self.__extract_date(str(date))

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

    def __date_file(self):
        '''
        if the file exist but it's not dated, we date it
        '''
        if len(self.date_list) == 0:
            path_in = self.file_path + "/" + self.get_file_name()
            if os.path.isfile(path_in):
                path_out = self.file_path + "/" + self.get_full_file_name()
                os.rename(path_in, path_out)

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

    def get_full_file_name(self):
        '''
        get the full name (date tag and file_name) of sur current file
        '''
        return(
            str(self.get_date()) +
            "_" +
            str(self.get_file_name())
            )

    def __str__(self):
        return(self.get_file_path() + "/" + self.get_full_file_name())


class Dated_file_list:
    def __init__(self, file_name_list, date_list=None):
        self.file_name_list = file_name_list
        self.date_list = date_list
        self.file_date_list = list()
        self.__load_files()

    def __list_len(self):
        if type(self.date_list) is list and len(self.date_list) > 1:
            list_len = min(len(self.file_name_list), len(self.date_list))
        else:
            list_len = len(self.file_name_list)
        return(int(list_len))

    def __load_files(self):
        for i in range(self.__list_len()):
            if type(self.date_list) is list:
                date = self.date_list[i]
            else:
                date = self.date_list
            self.file_date_list.append(
                Dated_file(self.file_name_list[i], date)
            )

    def __getitem__(self, key):
        if key >= 0 and key <= len(self.file_date_list):
            return(self.file_date_list[key])
        else:
            return(None)

    def __str__(self):
        output = ""
        for i in range(self.__list_len()):
            output += str(self[i])
        return(output)


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
        nargs='*')
    parser.add_argument(
        "-d",
        help="date to write if not file with this date exist, otherwise return \
        the file corresponding to this date.",
        default=None,
        action="store",
        dest="date",
        required=False,
        nargs='*')
    parser.add_argument(
        "-v",
        help="version information",
        default=False,
        action="store_true",
        dest="version")

    args = parser.parse_args()

    if args.version:
        print("0.0.1")
        exit(0)

    files_handle = Dated_file_list(args.input_file, args.date)
    print(str(files_handle))
