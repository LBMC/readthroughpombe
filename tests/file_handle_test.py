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

if __name__ == '__main__':
    unittest.main()
