#!/usr/bin/python3
# -*-coding:Utf-8 -*

import unittest
import datetime
import os
import sys
sys.path.append(os.path.abspath("src/func/"))
from file_handle import Dated_file


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

if __name__ == '__main__':
    unittest.main()
