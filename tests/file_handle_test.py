#!/usr/bin/python3
# -*-coding:Utf-8 -*
# Copyright laurent modolo for the LBMC UMR 5239 ©.
# contributor(s) : laurent modolo (2017)
#
# laurent.modolo@ens-lyon.fr
#
# This software is a computer program whose purpose is to test the
# src/func/file_handle.py program of this project
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

import unittest
import datetime
import os
import sys
sys.path.append(os.path.abspath("src/func/"))
from file_handle import Dated_file
from file_handle import Dated_file_list


class Dated_file_TestCase(unittest.TestCase):
    def test_test_no_date_found(self):
        '''
        if there is no date and older file don't exist we should set the \
        current date.
        '''
        datefile = Dated_file("/path/test_file.txt")
        current_date = datetime.date.today()
        self.assertEqual(
            datefile.get_date(),
            current_date.strftime("%Y_%m_%d")
        )
        self.assertEqual(
            datefile.get_full_file_name(),
            current_date.strftime("%Y_%m_%d_") + "test_file.txt"
        )

    def test_test_date_found(self):
        '''
        if there is no date and older file don't exist we should set the \
        current date.
        '''
        datefile = Dated_file("/path/2003_10_02_test_file.txt")
        self.assertEqual(
            datefile.get_date(),
            "2003_10_02"
        )

    def test_date_setting(self):
        datefile = Dated_file("/path/2003_10_02_test_file.txt")
        datefile.set_date("2005_11_04")
        self.assertEqual(
            datefile.get_date(),
            "2005_11_04"
        )

    def test_tuncate_file_name(self):
        datefile = Dated_file("/path/2003_10_02_test_file.txt")
        self.assertEqual(
            datefile.get_file_name(),
            "test_file.txt"
        )
        self.assertEqual(
            datefile.get_full_file_name(),
            "2003_10_02_test_file.txt"
        )

    def test_abs_path(self):
        datefile = Dated_file("/path/2003_10_02_test_file.txt")
        self.assertEqual(
            datefile.get_file_path(),
            os.path.abspath("/path/")
        )

    def test_list_files_empty(self):
        datefile = Dated_file("./data/examples/2004_10_02_test_file.txt")
        self.assertEqual(
            datefile[2],
            "2004_10_02_test_file.txt"
        )
        self.assertEqual(
            datefile[1],
            "2004_12_02_test_file.txt"
        )
        self.assertEqual(
            datefile[0],
            "2006_02_08_test_file.txt"
        )

    def test_set_to_last_existing_file(self):
        datefile = Dated_file("./data/examples/test_file.txt")
        self.assertEqual(
            datefile.get_date(),
            "2006_02_08"
        )
        self.assertEqual(
            datefile.get_file_name(),
            "test_file.txt"
        )
        self.assertEqual(
            datefile.get_full_file_name(),
            "2006_02_08_test_file.txt"
        )

    def test_set_to_dated_existing_file(self):
        datefile = Dated_file("./data/examples/2004_12_02_test_file.txt")
        self.assertEqual(
            datefile.get_date(),
            "2004_12_02"
        )
        self.assertEqual(
            datefile.get_file_name(),
            "test_file.txt"
        )
        self.assertEqual(
            datefile.get_full_file_name(),
            "2004_12_02_test_file.txt"
        )

    def test_date_existed_file(self):
        with open(os.path.abspath("./data/examples/test_file2.txt"), 'w'):
            datefile = Dated_file("./data/examples/test_file2.txt")
            datefile = Dated_file("./data/examples/test_file2.txt")
            current_date = datetime.date.today()
            self.assertTrue(
                os.path.isfile(
                    "./data/examples/" +
                    current_date.strftime("%Y_%m_%d_") +
                    "test_file2.txt"))
            self.assertEqual(
                datefile.get_full_file_name(),
                current_date.strftime("%Y_%m_%d_") + "test_file2.txt"
            )
            os.remove(os.path.abspath(
                "./data/examples/" +
                current_date.strftime("%Y_%m_%d_") +
                "test_file2.txt"))

    def test_set_date_existed_file(self):
        with open(os.path.abspath("./data/examples/test_file2.txt"), 'w'):
            datefile = Dated_file(
                "./data/examples/test_file2.txt",
                "2008_12_02")
            datefile = Dated_file("./data/examples/2008_12_02_test_file2.txt")
            self.assertEqual(
                datefile.get_full_file_name(),
                "2008_12_02_test_file2.txt"
            )
            os.remove(os.path.abspath(
                "./data/examples/2008_12_02_test_file2.txt"))


class Dated_file_list_TestCase(unittest.TestCase):
    def test_read_list(self):
        file_list = [
            "./data/examples/2004_10_02_test_file.txt",
            "./data/examples/2004_12_02_test_file.txt",
            "./data/examples/2006_02_08_test_file.txt"]
        datefile_list = Dated_file_list(file_list)
        for i in range(len(file_list)):
            self.assertEqual(
                datefile_list[i].get_full_file_name(),
                os.path.basename(file_list[i])
            )

    def test_read_list_with_date(self):
        file_list = [
            "./data/examples/test_file.txt",
            "./data/examples/test_file.txt",
            "./data/examples/test_file.txt"]
        date_list = ["2004_10_02", "2004_12_02", "2006_02_08"]
        datefile_list = Dated_file_list(file_list, date_list)
        for i in range(len(file_list)):
            self.assertEqual(
                datefile_list[i].get_full_file_name(),
                date_list[i] + "_" + os.path.basename(file_list[i])
            )

    def test_read_list_with_wildcard(self):
        file_list = ["./data/examples/test_file.*"]
        file_list_check = [
            "./data/examples/2008_04_12_test_file.csv",
            "./data/examples/2006_02_08_test_file.txt"]
        date_list = list()
        datefile_list = Dated_file_list(file_list, date_list)
        for i in range(len(file_list_check)):
            self.assertEqual(
                datefile_list[i].get_full_file_name(),
                os.path.basename(file_list_check[i])
            )

    def test_read_list_with_wildcard_and_date(self):
        file_list = ["./data/examples/test_file.*"]
        file_list_check = [
            "./data/examples/2008_04_12_test_file.csv",
            "./data/examples/2004_12_02_test_file.txt"]
        date_list = ["2008_04_12", "2004_12_02"]
        datefile_list = Dated_file_list(file_list, date_list)
        for i in range(len(file_list_check)):
            self.assertEqual(
                datefile_list[i].get_full_file_name(),
                os.path.basename(file_list_check[i])
            )


if __name__ == '__main__':
    unittest.main()
