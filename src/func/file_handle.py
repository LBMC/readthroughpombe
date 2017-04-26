#!/usr/bin/python3
# -*-coding:Utf-8 -*
# Copyright laurent modolo for the LBMC UMR 5239 ©.
# contributor(s) : laurent modolo (2017)
#
# laurent.modolo@ens-lyon.fr
#
# This software is a computer program whose purpose is to manage dated file
# names in complience with the bioinformatic good practices used in the LBMC.
#
# This software is governed by the CeCILL  license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

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
    def __init__(self, file_name, date=None, redate=False, escape=False):
        file_name = os.path.abspath(str(file_name))
        self.redate = redate
        self.escape = escape
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
        if self.redate:
            path_in = self.file_path + "/" + self.get_full_file_name()
            self.date = datetime.date.today()
            path_out = self.file_path + "/" + self.get_full_file_name()
            if os.path.isfile(path_in):
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
        if self.escape:
            return(
                self.get_file_path().replace(" ", "\ ") +
                "/" + self.get_full_file_name().replace(" ", "\ "))
        else:
            return(
                '"' + self.get_file_path() +
                "/" + self.get_full_file_name() + '"')

    def check(self):
        """            print("test")

        check if file exist
        """
        file_name = self.get_file_path() + \
            "/" + self.get_full_file_name()
        if os.path.isfile(file_name) and file_name != "":
            return True
        else:
            print("error: " + file_name + " not found")
            return False


class Dated_file_list:
    def __init__(
            self, file_name_list, date_list=list(),
            check=False, redate=False, escape=False):
        self.check = check
        self.redate = redate
        self.escape = escape
        self.file_name_list = list()
        self.file_name_list.extend(file_name_list)
        self.date_list = list()
        self.date_list.extend(date_list)
        self.__expend_files_name_from_wildcard()
        self.file_date_list = list()
        self.__load_files()

    def __handle_wildcard(self, file_name):
        '''
        list file corresponding to the UNIX wildcard of a file_name
        '''
        file_name = str(file_name)
        file_name = os.path.abspath(file_name)
        file_path = os.path.dirname(str(file_name))
        file_name = os.path.basename(str(file_name))
        file_names = glob.glob(file_path + "/*" + file_name)
        if self.check and len(file_names) == 0:
            print("error: no file found")
            exit(1)
        for i in range(len(file_names)):
            if self.check and not os.path.isfile(file_names[i]):
                print("error: " + str(file_names[i]) + " not found")
                exit(1)
            file_names[i] = os.path.basename(file_names[i])
            # if we don't have specified the date in the file name we remove it
            if re.match(r"\d{4}\_\d{2}\_\d{2}\_", file_name) is None:
                file_names[i] = re.search(
                    r"(\d{4}\_\d{2}\_\d{2}\_){0,1}(.*)",
                    file_names[i]).group(2)
            file_names[i] = file_path + "/" + file_names[i]
        file_names = list(set(file_names))
        file_names.sort()
        return(file_names)

    def __expend_files_name_from_wildcard(self):
        '''
        get all the files corresponding to the UNIX wildcard for a list of
        file_names
        '''
        if not isinstance(self.file_name_list, list):
            print(type(self.file_name_list))
            print("error : file_name_list is not a list")
        if not isinstance(self.date_list, list):
            print(type(self.date_list))
            print("error : date_list is not a list")
        file_name_list = list()
        for i in range(self.__list_len()):
            file_name_list.extend(
                self.__handle_wildcard(self.file_name_list[i]))
        self.file_name_list = file_name_list[:]

    def __list_len(self):
        if type(self.date_list) is list and len(self.date_list) > 1:
            list_len = min(len(self.file_name_list), len(self.date_list))
        else:
            list_len = len(self.file_name_list)
        return(int(list_len))

    def __load_files(self):
        for i in range(self.__list_len()):
            if len(self.date_list) > 1:
                date = self.date_list[i]
            elif len(self.date_list) == 1:
                date = self.date_list[0]
            else:
                date = None
            self.file_date_list.append(
                Dated_file(
                    self.file_name_list[i],
                    date,
                    self.redate,
                    self.escape)
            )
            if self.check:
                if not self.file_date_list[len(self.file_date_list)-1].check():
                    print("error")
                    exit()

    def __getitem__(self, key):
        if key >= 0 and key <= len(self.file_date_list):
            return(self.file_date_list[key])
        else:
            return(None)

    def __str__(self):
        output = ""
        for i in range(self.__list_len()):
            if i > 0:
                if self.escape:
                    output += "\\n"
                else:
                    output += " "
            output += str(self[i])
        return(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog="file_handle.py",
        description="script to handle date in file name in the format \
        'yyyy_mm_dd_filemames'. By default return the last file corresponding \
        to the filename if existing or create it with the current date.")
    parser.add_argument(
        "-f", "--file",
        help="input filename",
        default=None,
        action="store",
        dest="input_file",
        required=True,
        nargs='+')
    parser.add_argument(
        "-d", "--date",
        help="date to write if not file with this date exist, otherwise return \
        the file corresponding to this date.",
        default=list(),
        action="store",
        dest="date",
        required=False,
        nargs='*')
    parser.add_argument(
        "-r", "--redate",
        help="date file with the current date even if dated",
        default=False,
        action="store_true",
        dest="redate",
        required=False)
    parser.add_argument(
        "-c", "--check",
        help="Return an error if the dated file is not found",
        default=False,
        action="store_true",
        dest="check",
        required=False)
    parser.add_argument(
        "-e", "--escape",
        help="escape file name instead of coting them",
        default=False,
        action="store_true",
        dest="escape",
        required=False)
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

    files_handle = Dated_file_list(
        args.input_file,
        args.date,
        args.check,
        args.redate,
        args.escape)
    print(str(files_handle))
