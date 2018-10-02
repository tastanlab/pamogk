#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from clinical_ov import read_ov

def read_clinical(dtype):
    if dtype == 'ov':
        read_ov()
    elif dtype == 'brca':
        read_brca()
    elif dtype == 'kirc':
        read_kirc()
    elif dtype == 'coad':
        read_coad()

if __name__ == '__main__':
    read_clinical(dtype='ov')
