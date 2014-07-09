#!/usr/bin/env python
# -*- coding: utf-8 -*-
#==============================================================================
#title           :aggregation.py
#description     :This function reduces the resolution of a raster aggregating 
#                 pixels by a multiple factor of cells (integer > 1) according 
#                 to one of these statistical parameters: mean, max, min, sum 
#                 and standard deviation.
#author          :Giovanni Caudullo
#date            :20140626
#version         :0.1
#usage           :$python aggregation.py input.tif output.tif cellfactor
#notes           :Output is a float 32 TIFF raster
#python_version  :2.7.7
#license         :Creative Commons Attribution-ShareAlike 3.0 Unported License
#==============================================================================

import numpy
import argparse
from osgeo import gdal


def ras_aggregation(inras, outras, cellfact, aggr_type='AVG', expand=True, prop_nodata=True):
    gdal.AllRegister()
    
    functions = {'MAX':'max', 'MIN':'min', 'AVG':'mean', 'SUM':'sum', 'STD':'std'}
    if prop_nodata==True:
        ign = ''
    else:
        ign = 'nan'
    
    ## GET INPUT INFO
    dataset = gdal.Open(inras, gdal.GA_ReadOnly)
    geot = dataset.GetGeoTransform()
    proj = dataset.GetProjection()
    band = dataset.GetRasterBand(1)
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    ND = band.GetNoDataValue()
    if ND is None:
        ND = -3.4028234663852886e+38
    
    ## SET OUTPUT EXTENT
    if expand == True:
        outRows = (rows + cellfact - rows%cellfact)/cellfact
        outCols = (cols + cellfact - cols%cellfact)/cellfact
    else:
        outRows = (rows - rows%cellfact)/cellfact
        outCols = (cols - cols%cellfact)/cellfact
    
    ## CREATE OUTPUT
    driver = dataset.GetDriver()
    outDataset = driver.Create(outras, outCols, outRows, 1, gdal.GDT_Float32)
    outDataset.SetGeoTransform((geot[0], geot[1]*cellfact, geot[2], geot[3], geot[4], geot[5]*cellfact))
    outDataset.SetProjection(proj)
    outBand = outDataset.GetRasterBand(1)
    outBand.SetNoDataValue(ND)
    
    ## BLOCK ITERATION
    xBlockSize = 200 * cellfact
    yBlockSize = 200 * cellfact
    for i in range(0, rows, yBlockSize):
        if i + yBlockSize < rows:
            numRows = yBlockSize
        else:
            numRows = rows - i
        for j in range(0, cols, xBlockSize):
            if j + xBlockSize < cols:
                numCols = xBlockSize
            else:
                numCols = cols - j
    
            ## INPUT BLOCK ARRAY
            array = band.ReadAsArray(j , i, numCols, numRows)
            data = numpy.where(array==ND, numpy.nan, array)
    
            ## EXPAND OR REMOVE ROWS AND COLS
            if numRows%cellfact != 0:
                if expand == True:
                    diffRows = cellfact - numRows%cellfact
                    addRows = numpy.empty((diffRows, numCols))
                    addRows[:] = numpy.nan
                    data = numpy.append(data, addRows, 0)
                else:
                    data = data[0:-(numRows%cellfact),:]
            tmpRows, tmpCols = data.shape
    
            if numCols%cellfact != 0:
                if expand == True:
                    diffCols = cellfact - numCols%cellfact
                    addCols = numpy.empty((tmpRows, diffCols))
                    addCols[:] = numpy.nan
                    data = numpy.append(data, addCols, 1)
                else:
                    data = data[:,0:-(numCols%cellfact)]
    
            #olderr = numpy.seterr(all='ignore')
    
            ## AGGREGATION
            nrows, ncols = data.shape
            aggr_array = numpy.array([[]])
            for x in range(0, nrows, cellfact):
                for y in range(0, ncols, cellfact):
                    cmd = "func = numpy.%s%s(data[x:cellfact+x,y:cellfact+y])"
                    exec cmd % (ign, functions[aggr_type])
                    aggr_array = numpy.append(aggr_array, func)
            aggr_array = aggr_array.reshape(nrows/cellfact, ncols/cellfact)
    
            ## BLOCK WRITING
            outDataset.GetRasterBand(1).WriteArray(aggr_array, j/cellfact, i/cellfact)
    
    outDataset = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='aggregation.py', description= 
    '''This function reduces the resolution of a raster aggregating pixels by a
       multiple factor of cells (integer > 1) according to one of these statistical 
       parameters: mean, max, min, sum and standard deviation.''')
    
    parser.add_argument('inraster', help='input raster file')
    parser.add_argument('outraster', help='output TIFF raster (as float 32)')
    parser.add_argument('cellfactor', type=int, nargs=1, \
                         help='factor of cell multiplication (i.e. output cell ' + \
                              'is "cellfactor" times larger than input cell)')
    parser.add_argument('-at', dest='aggreg', choices=['MAX','MIN','AVG','SUM','STD'], default='AVG', \
                        help='aggregation technique: AVG for mean (default), SUM for sum, ' +\
                             'MAX for maximum, MIN for minimum, STD for standard deviation')
    parser.add_argument('-e', dest='exp', action='store_true', default=False, \
                        help='expands of the output to fit a multiple of the cell factor, '+\
                             'otherwise output is reduced (default)')
    parser.add_argument('-nd', dest='nodata', action='store_true', default=False, \
                        help='propagates NoData in cell computation, otherwise NoData are ognored (default)')
    
    args = parser.parse_args()    

    ras_aggregation(args.inraster, args.outraster, args.cellfactor[0], \
                    aggr_type=args.aggreg, expand=args.exp, prop_nodata=args.nodata)
