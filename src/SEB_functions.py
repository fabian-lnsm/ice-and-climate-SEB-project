#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 17:41:56 2023

@author: berg
"""

import datetime as dt
import numpy as np
import re
import os

# this file contains
# Methods to read and process data obtained from the UU/IMAU AWS network
# Notably from stations on the K-transect, Greenland, and Antarctic Peninsula.
# These data are published (or will be soon) on the Pangaea data repository

# the data themselves are received form M. van Tiggelen by personal communication.

# the included functions are
# the class SEB_data
# more info? -> help(SEB_data)


version = "2.1"
class SEB_data:
    '''Version 2.1 of a python class that reads and organizes SEB model output,
    derived from observations from the K-transect, Greenland, and Antarctic Peninsula'''
    

    def __init__(self, FileName=""):
        '''For the initialisation of this class, only the filename (including path) is needed.
        During the initialisation all data is read.
        The function works on all four provided data sets.
        
        Output is a SEB_data-class object, containing the variables:
            ok             bool     data properly readed
            version        string   version of SEB_data
            FileName       string   name of the SEB data file, without path
            nval           integer  number of 
        '''
        
        self.ok = False
        self.version = "SEB_hourly_data version "+version
        self.FileName = re.split(r'/+', FileName)[-1] # remove path
        
        if FileName=="":
            print("A FileName is required!")
            return
        if not os.path.isfile(FileName):
            print("SEB data file "+FileName+" does not exist, return")
            return
        
        # first count number of lines and read the header
        with open(FileName) as AWSfile:
            self.nval = -1
            for line in AWSfile:
                li=line.strip()
                if not li.startswith("#"):
                    if self.nval == -1:
                        # remove trailing line break and split line on the spaces
                        splitted_header = (line.strip()).split(",") 
                        splitted_filename = (FileName.strip()).split("_")
                        print("The header has {0:3d} entries.".format(np.size(splitted_header)))    
                        if splitted_header[0]=="time":                            
                            if splitted_filename[2]=="hourly":
                                self.TimeStep = 3600.
                                self.Hourly   = True
                                DTGFormat = '%Y-%m-%d %H:%M:%S'
#                                ivstart = 3
#                                ntime   = 3
                            else:    # assume daily data
                                self.TimeStep = 86400.
                                self.Hourly   = False
                                DTGFormat = '%Y-%m-%d'
#                                ivstart = 2
#                                ntime   = 2
                            ivstart = 2
#                            ntime   = 2
                        else:
                            print("We assumed that every header starts with 'time' (and comma delimited data)")
                            print("However, this file starts with '"+splitted_header[0]+
                                  "' and '"+splitted_header[1]+"'.")
                            print("Fix this! (or ask someone to fix this)")
                            return
    
                        if splitted_header[ivstart]=="Time":
                            ivstart+=1 # neglect this entry
                            
                        self.nvar = np.size(splitted_header) - ivstart
                        self.Variables = splitted_header[ivstart:]
    
                    self.nval += 1
            print("AWS file '{0:s}' has {1:5d} lines of data for {2:2d} variables, start reading it.".format(self.FileName, 
                        self.nval, self.nvar))

        
        self.AllData    = np.zeros( [ self.nvar, self.nval ])                
        self.DateTime   = np.zeros( self.nval, dtype=dt.datetime )
        self.nvarDate   = np.zeros( self.nval, dtype=int)    # the number of variables per time entry. Should be constant, but isn't for AWS5 data :-(
        self.varvalid   = np.ones( self.nval, dtype=bool )
        
        # read data
        print("now start reading the data")
        with open(FileName) as AWSfile:
            ival = -1
            for line in AWSfile:
                li=line.strip()
                if not li.startswith("#"):
                    if ival>=0:
                        splitted_line = (line.strip()).split(",")
                        splitted_line = ["-999." if x == '' else x for x in splitted_line]

                        # the time stap has errors, neglect is as a whole except first entry
                        self.DateTime[ival] = dt.datetime.strptime(splitted_line[0], DTGFormat)
                       
                        if np.size(splitted_line)-ivstart == self.nvar:
                            self.AllData[:, ival] = [float(value) for value in splitted_line[ivstart:]]
                            self.nvarDate[ival]   = self.nvar
 
                        else:
                            print("Errors should not occur in other dataset than the one for S6.")
                    ival += 1
      
        print("Reading completed.")
        
        # turn invalid data into NaN
        self.AllData = np.where(self.AllData==-999., np.NaN, self.AllData)

        return
    
    def List_Variables(self):
        print("  #: Variable")
        for ivar in range(self.nvar):
            if self.varvalid[ivar]:
                print("{0:3d}:  {1:11s}".format(ivar+1, self.Variables[ivar]))
            else:
                print("{0:3d}: ERASED {1:11s}".format(ivar+1, self.Variables[ivar]))
    
    
    def Find_Variable(self, VarName):
        '''This function finds the entry number of a Variable in the AllData array. 
        Please note this function is case sensitive - your request should be exactly matching.'''
        
        lok      = True
        VarIndex = -1
        for v in range(self.nvar):
            if VarName.strip() == self.Variables[v].strip():
                VarIndex = v
        if VarIndex == -1:
            print("Variable '"+VarName+"' not found.")
            lok = False
        return VarIndex, lok  


    def Extract_Variable(self, VarName, NoDataForFail=True):
        '''This function extracts the data of a Variable from the main dataset AllData.
        Unless the optional parameter NoDataForFail is set to False, no data is 
          given if the requested variable name is not found.'''
        
        VarIndex, lok = self.Find_Variable(VarName)
        if not self.varvalid[VarIndex] and lok:
            print("WARNING: The data for variable '"+VarName+"' has been erased due to consistency problems.")
            print("         NO DATA IS PROVIDED")

        if lok and self.varvalid[VarIndex]:
            return self.AllData[VarIndex, :]
        elif NoDataForFail:
            return
        else:
            return np.zeros(self.nval)*np.NaN
#%%                
def get_doy(year, month, day):
    timediff = dt.datetime(year=year, month=month, day=day) - dt.datetime(year=year, month=1, day=1)
    return timediff.days + 1

def convert_T_in_LWout(Tskin, Celcius=True):
    if Celcius:
        Tadd = 273.16
    else:
        Tadd = 0.
        
    return -( (Tskin+Tadd)**4 ) * 5.67E-8

def convert_LWout_in_T(LWout, Celcius=True):
    if Celcius:
        Tadd = 273.16
    else:
        Tadd = 0.
              
    return ( -LWout/5.67E-8 )**(0.25) - Tadd

def get_running_melt_sum(MeltE, TimeStep, ResetAtNan=True):
    '''The meltsum is in m w.e. per year'''
    RmeltSum = np.zeros(np.size(MeltE))
    for i in range(np.size(MeltE)):
        if not np.isnan(MeltE[i]):
            RmeltSum[i] = RmeltSum[i-1] + MeltE[i]*TimeStep/334000.
        elif not ResetAtNan:
            RmeltSum[i] = RmeltSum[i-1]
        elif np.isnan(MeltE[i-1]):
            RmeltSum[i-1] = np.NaN
    
    RmeltSum = RmeltSum/1000.
    
    return RmeltSum
            
def get_daily_average(Var, DateTime, GiveDatesBack=False, PrintInfo=False):
    '''This function derives the daily averages of a variable.
    It start at the first full day.'''
    Istart = (24-DateTime[0].hour)%24
    Ndays  = (DateTime[-1]-DateTime[0]).days
    
    if PrintInfo:
        print("The date and hour of the first entry is {}.".format(DateTime[0].strftime("%B %d, %Y, %H:%M")))
        print("The date and hour of the last entry is {}.".format(DateTime[-1].strftime("%B %d, %Y, %H:%M")))
        print("The dataset contains data of {0:5d} days.".format(Ndays))

# In order to do this in one call, I reshape the array into a 2D array and average over the second axis, thus the data of one day.    
    if Ndays+1 == len(Var):
        VarDay = Var[Istart:Istart+Ndays:1]
        dh = 1
    else:
        VarDay = np.mean(np.reshape(Var[Istart:Istart+Ndays*24], [Ndays, 24]),1)
        dh = 24
    
    if GiveDatesBack:
        return VarDay, DateTime[Istart:Istart+Ndays*dh:dh]
    else:
        return VarDay
    
def get_next_month(DateTime):
    if DateTime.month==12:
        return dt.datetime(year=DateTime.year+1, month=1, day=1)
    else:
        return dt.datetime(year=DateTime.year, month=DateTime.month+1, day=1)
    
    
def get_monthly_average(Var, DateTime, GiveDatesBack=False, PrintInfo=False):
    '''This function derives the monthly averages of a variable.
    It start at the first full day.'''
    
    # find the start of the first full month
    
    dtime = (DateTime[1] - DateTime[0]).seconds/3600 
    dh = 1
    if dtime == 1:
        dh = 24
        
    NextMonth = get_next_month(DateTime[0])
    iStart = int( (NextMonth-DateTime[0]).days*dh + (NextMonth-DateTime[0]).seconds/3600 )

    # find the number of months
    nData  = np.size(DateTime)
    iMS    = iStart
    iME    = iMS+(get_next_month(DateTime[iMS])-DateTime[iMS]).days * dh
    nMonth = 0
       
    while iME<nData:
        iMS     = iME
        iME    += (get_next_month(DateTime[iME])-DateTime[iME]).days * dh
        nMonth += 1

    iFinal = iMS    
    
    if PrintInfo:
        print("The date and hour of the first entry is {}.".format(DateTime[0].strftime("%B %d, %Y, %H:%M")))
        print("First full month starts at {}".format(DateTime[iStart].strftime("%B %d, %Y, %H:%M")))
        print("The date and hour of the last entry is {}.".format(DateTime[-1].strftime("%B %d, %Y, %H:%M")))
        print("Last full month end at {}".format(DateTime[iFinal-1].strftime("%B %d, %Y, %H:%M")))
        print("The dataset contains data of {0:3d} months.".format(nMonth))
    
    VarMonth = np.zeros(nMonth)
    VarDate  = np.zeros(nMonth, dtype=dt.datetime )
    iMS = iStart
    for iMonth in range(nMonth):
        iME = iMS + (get_next_month(DateTime[iMS])-DateTime[iMS]).days * dh
        VarMonth[iMonth] = np.mean(Var[iMS:iME])
        VarDate[iMonth]  = DateTime[iMS]
        iMS = iME
        
    if GiveDatesBack:
        return VarMonth, VarDate
    else:
        return VarMonth
    
    
def get_avg_monthly_value(Var, DateTime, PrintInfo=False, GiveNumberOfMonths=False):
    '''This function derives the mean monthly value of a variable'''
    VarMonth, VarDate = get_monthly_average(Var, DateTime, GiveDatesBack=True, PrintInfo=PrintInfo)
    
    VarMeanMonth = np.zeros(14)
    Ndata        = np.zeros(13)
    for iMonth in range(np.size(VarMonth)):
        iM = VarDate[iMonth].month
        # leave out months with no valid data
        if not np.isnan(VarMonth[iMonth]):
            VarMeanMonth[iM] += VarMonth[iMonth]
            Ndata[iM]        += 1
    
    VarMeanMonth[1:13] = VarMeanMonth[1:13]/Ndata[1:13]
    VarMeanMonth[0]    = VarMeanMonth[12]
    VarMeanMonth[13]   = VarMeanMonth[1]
    
    if GiveNumberOfMonths:
        return VarMeanMonth, Ndata
    else:
        return VarMeanMonth
   
#%%
def get_epsilon_from_LWin_Temp(LWin,Temp,Celcius=True):
    if Celcius:
        Tadd = 273.16
    else:
        Tadd = 0.
    
    return LWin/(5.67E-8 * (Temp+Tadd)**4)

