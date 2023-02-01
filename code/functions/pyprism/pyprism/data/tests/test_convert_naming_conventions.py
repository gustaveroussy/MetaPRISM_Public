# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat

Tests for _convert_naming_conventions.py module.
"""

from pyprism.data import any_to_cu, any_to_cw, any_to_mc, any_to_us

def test_any_to_cu():
    # convert CapitalizedWords to Capitalized_Under_Scores
    assert any_to_cu("PrintHTML") == "Print_Html"
    assert any_to_cu("IOError") == "Io_Error"
    assert any_to_cu("SetXYPosition") == "Set_Xy_Position"

    # convert mixedCase to Capitalized_Under_Scores
    assert any_to_cu("printHTML") == "Print_Html"
    assert any_to_cu("ioError") == "Io_Error"
    assert any_to_cu("setXYPosition") == "Set_Xy_Position"

    # convert other to Capitalized_Under_Scores
    assert any_to_cu("DFI.time") == "Dfi_Time"
    assert any_to_cu("DFI-time") == "Dfi_Time"

    # keep identical
    assert any_to_cu("Print_Html") == "Print_Html"
    assert any_to_cu("Set_Xy_Position") == "Set_Xy_Position"


def test_any_to_cw():
    # convert mixedCase to Capitalized_Under_Scores
    assert any_to_cw("printHTML") == "PrintHtml"
    assert any_to_cw("ioError") == "IoError"
    assert any_to_cw("setXYPosition") == "SetXyPosition"

    # convert Capitalized_Under_Scores to CapitalizedWords
    assert any_to_cw("Print_Html") == "PrintHtml"
    assert any_to_cw("Set_Xy_Position") == "SetXyPosition"

    # convert other to CapitalizedWords
    assert any_to_cw("DFI.time") == "DfiTime"
    assert any_to_cw("DFI-time") == "DfiTime"

    # keep identical
    assert any_to_cw("PrintHTML") == "PrintHtml"
    assert any_to_cw("IoError") == "IoError"
    assert any_to_cw("SetXYPosition") == "SetXyPosition"


def test_any_to_mc():
    # convert Capitalized_Under_Scores to mixedCase
    assert any_to_mc("Print_Html") == "printHtml"
    assert any_to_mc("Set_Xy_Position") == "setXyPosition"

    # convert CapitalizedWords to mixedCase
    assert any_to_mc("PrintHTML") == "printHtml"
    assert any_to_mc("IoError") == "ioError"
    assert any_to_mc("SetXYPosition") == "setXyPosition"

    # convert other to CapitalizedWords
    assert any_to_mc("DFI.time") == "dfiTime"
    assert any_to_mc("DFI-time") == "dfiTime"

    # keep identical
    assert any_to_mc("printHTML") == "printHtml"
    assert any_to_mc("ioError") == "ioError"
    assert any_to_mc("setXyPosition") == "setXyPosition"



def test_any_to_us():
    # convert CapitalizedWords to underscores
    assert any_to_us("PrintHTML") == "print_html"
    assert any_to_us("IOError") == "io_error"
    assert any_to_us("SetXYPosition") == "set_xy_position"

    # convert mixedCase to underscores
    assert any_to_us("printHTML") == "print_html"
    assert any_to_us("ioError") == "io_error"
    assert any_to_us("setXYPosition") == "set_xy_position"

    # convert other to underscores
    assert any_to_us("DFI.time") == "dfi_time"
    assert any_to_us("DFI-time") == "dfi_time"

    # keep identical
    assert any_to_us("print_html") == "print_html"
    assert any_to_us("set_xy_position") == "set_xy_position"
