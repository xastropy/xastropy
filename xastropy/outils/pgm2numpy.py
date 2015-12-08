from __future__ import print_function, absolute_import, division, unicode_literals

"""
Created on 2011-3-3

@author: summit

Taken by JXP from :
https://code.google.com/p/jolly/source/browse/trunk/jolly/src/jolly/imp/pgm2numpy.py?spec=svn43&r=41

Slightly modified for format of K1DM3 data
"""

import re
import numpy
import math
import array

def pgm_nextline(fin):
    while True:
        header  = fin.readline().strip()
        if header.startswith('#'):
            continue
        else:
            break
    return header
            
def pgm2numpy_p5(fin, debug = True):
    try:
        
        rows, cols = 0, 0
        raw_data = re.split('\W+', pgm_nextline(fin))
        rows = int(raw_data[1])
        cols = int(raw_data[0])
        raw_data = re.split('\W+', pgm_nextline(fin))
        bits = (int(raw_data[0])+1)/256
        
        assert (rows, cols) != (0,0)
        
        if debug:
            print('Rows: %d, cols: %d' %(rows, cols))
        
        assert bits > 0
        
        if debug:
            print('color bit is %d' % bits)
        
        if bits<=1:
            binvalues = array.array('c') # unsigned char 
        elif bits<=2 and bits>1:
            binvalues = array.array('H') # unsigned short
        elif bits<=4 and bits>2:
            binvalues = array.array('f') # float
        else:
            print('color bit is %d too big' % bits)
        result = numpy.fromfile(fin, dtype=numpy.uint8)
        result = numpy.reshape(result, (rows, cols))
        
        return result
            
    finally:
        if fin != None:
            fin.close()
        fin = None
    
    return None
    
    
def pgm2numpy_p2(fin, debug = True, K1DM3=True):
    '''Algorithm to parse a PGM P2 filename
    I also modified the minimum image type to be int16
    Parameters:
    ----------
    fin: file 
    K1DM3: bool, optional [True]
      Removes first and last characters which are '' for 
      K1DM3 data (Alignment Telesope)

    Returns:
    ---------
    numpy array of the image
    '''

    try:
        rows, cols = 0, 0
        
        raw_data = re.split('\W+', pgm_nextline(fin))
        
        rows = int(raw_data[1])
        cols = int(raw_data[0])
        raw_data = re.split('\W+', pgm_nextline(fin))
        bits = (int(raw_data[0])+1)/256
        
        assert (rows, cols) != (0,0)
        
        if debug:
            print('Rows: %d, cols: %d' %(rows, cols))
        
        assert bits > 0
        
        if debug:
            print('color bit is %d' % bits)
        #
        # Initialise a 2D numpy array
        #
        if bits<=1:
            #result = numpy.zeros((rows, cols), numpy.int8) # Original
            result = numpy.zeros((rows, cols), numpy.int16)
        elif bits<=2 and bits>1:
            result = numpy.zeros((rows, cols), numpy.int16)
        elif bits<=4 and bits>2:
            result = numpy.zeros((rows, cols), numpy.int32)
        else:
            print('color bit is %d too big' % bits)
        pxs = []
        raw_data = re.split('\W+', fin.read())
        # Remove first and last character
        if K1DM3:
          raw_data.pop(0)
          raw_data.pop(-1)
        #
        # Read to EOF
        #
        #xdb.set_trace()
        for v in raw_data:
            if v == '':
                break
            if int(v) <0:
                xdb.set_trace()
            pxs.append(int(v))
        
        if len(pxs) != rows*cols:
            if debug:
                print('Insufficient image data:', len(pxs))
            return None
        
        for r in range(rows):
            for c in range(cols):
                #
                # Index into the numpy array and set the pixel value
                #
                result[r, c] = pxs[r*cols+c]
        return result
    
    finally:
        if fin != None:
            fin.close()
        fin = None
    
    return None

def pgm2numpy(filename):
    '''
    Read a PGM into a numpy array. 
    '''
    fin = None 
    debug = True
    r = None

    print('Reading: {:s}'.format(filename))
    fin = open(filename, 'rb')
    header = pgm_nextline(fin)
 
    if header == 'P1':
        pass
    elif header == 'P2':
        if debug:
            print('PBM is p2 format')
        r =  pgm2numpy_p2(fin, debug)
    elif header == 'P4':
        assert False, 'Raw PBM reading not implemented yet'
    elif header == 'P5':
        if debug:
            print('PBM is p5 format')
        r =  pgm2numpy_p5(fin, debug)
    else:
        #
        # unexpected header.
        #
        if debug:
            print('Bad mode: ', header)
        return None
    return r    
        
