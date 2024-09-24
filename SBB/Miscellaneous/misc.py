#!/bin/env/python
#! -*- coding: utf-8 -*-

def yes_or_no(question,default_ans='y'):
    while "the answer is invalid":
        reply = str(input(question+' (y/n):').encode('utf-8')).lower().strip() or default_ans
        if reply[:1] == 'y':
            return True
        elif reply[:1] == 'n':
            return False
        else :
            return True if default_ans=='y' else False