#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 13:35:25 2020

@author: hamzamahdi

This script generates a single MATLAB file containing all functions from the 
Modern Robotics Textbook: https://github.com/NxRLab/ModernRobotics/tree/master/packages/MATLAB
Place this script in a folder with all the matlab functions you would like to include.
This should work for any set of MATLAB functions and is not specific to
the repository referenced

To load the generated function in MATLAB, simply place loadFuns.m in your MATLAB directory
and type loadFuns in your script
"""


import glob

function_names = [f[:-2] for f in glob.glob("*.m")]

f= open("loadFuns.m","w+")


'''
create code to import the functions directly into the workspace using the assignin
function. Reference: https://stackoverflow.com/a/16870010
'''
f.write("function message = loadFuns\r\n")
for i, name in enumerate(function_names):
     f.write("  assignin('base','"+name+"',@"+name+");\n")

f.write("  message='Done importing functions to workspace';\n")
f.write("end\n")


# add the actual functions
for i, name in enumerate(function_names):
     m = open(name+'.m','r')
     contents =m.read()
     f.write(contents)
     f.write('\n')

f.close()
