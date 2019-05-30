# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 20:43:51 2018

@author: ErykD
"""

import numpy as np

def loadnpyfile(filename):
	with open(filename,'rb') as f:
		data = np.load(f,encoding='latin1')
	nd1 = len(data)
	nd2 = len(data[0])
	return (data,nd1,nd2)

import pickle

def loadpklfile(filename,chan):
	data={}
	with open(filename,'rb') as f:
		data = pickle.load(f,encoding='latin1')
	print(data['xTime'].dtype)
	return (data['xTime'], data['yWaves'][0], data['yWaves'][1+ chan%2])
