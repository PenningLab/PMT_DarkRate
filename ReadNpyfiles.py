# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 20:43:51 2018

@author: ErykD
"""

import numpy as np

def loadnpyfile(filename):
	data = np.load(filename)
	nd1 = len(data)
	nd2 = len(data[0])
	return (data,nd1,nd2)
