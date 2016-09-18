#!/usr/bin/python
# RasterClipper.py - clip a geospatial image using a shapefile

import operator
from osgeo import gdal, gdalnumeric, ogr, osr
import Image, ImageDraw
from numpy import *
import sys


# This function will convert the rasterized clipper shapefile 
# to a mask for use within GDAL.    
def imageToArray(i):
    """
    Converts a Python Imaging Library array to a 
    gdalnumeric image.
    """
    a=gdalnumeric.fromstring(i.tostring(),'b')
    a.shape=i.im.size[1], i.im.size[0]
    return a

def arrayToImage(a):
    """
    Converts a gdalnumeric array to a 
    Python Imaging Library Image.
    """
    i=Image.fromstring('L',(a.shape[1],a.shape[0]),
            (a.astype('b')).tostring())
    return i
     
def world2Pixel(geoMatrix, x, y):
  """
  Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
  the pixel location of a geospatial coordinate 
  """
  ulX = geoMatrix[0]
  ulY = geoMatrix[3]
  xDist = geoMatrix[1]
  yDist = geoMatrix[5]
  rtnX = geoMatrix[2]
  rtnY = geoMatrix[4]
  pixel = int((x - ulX) / xDist)
  line = int((ulY - y) / xDist)
  return (pixel, line) 

def histogram(a, bins=range(0,256)):
  """
  Histogram function for multi-dimensional array.
  a = array
  bins = range of numbers to match 
  """
  fa = a.flat
  n = gdalnumeric.searchsorted(gdalnumeric.sort(fa), bins)
  n = gdalnumeric.concatenate([n, [len(fa)]])
  hist = n[1:]-n[:-1] 
  return hist

def stretch(a):
  """
  Performs a histogram stretch on a gdalnumeric array image.
  """
  hist = histogram(a)
  im = arrayToImage(a)   
  lut = []
  for b in range(0, len(hist), 256):
    # step size
    step = reduce(operator.add, hist[b:b+256]) / 255
    # create equalization lookup table
    n = 0
    for i in range(256):
      lut.append(n / step)
      n = n + hist[i+b]
  im = im.point(lut)
  return imageToArray(im)

def Usage():
    print('Usage: raster_clipper.py raster_file clipping_shapefile clipped_raster_file[GeoTiff]')
    print('')

argv = sys.argv

#print argv
if len(argv) != 4:
    Usage()
    sys.exit(1)
    
# Raster image to clip
raster = argv[1]

# Polygon shapefile used to clip
shp = argv[2]

# Name of clip raster file(s)
output = argv[3]
     

# Raster image to clip
#raster = "E:/jjc/WebSites/WS_RasterClipper/App_Data/altitude.tif"

# Polygon shapefile used to clip
#shp = "E:/jjc/WebSites/WS_RasterClipper/App_Data/water.shp"

# Name of clip raster file(s)
#output = "E:/jjc/WebSites/WS_RasterClipper/App_Data/clip.tif"

# Load the source data as a gdalnumeric array
srcArray = gdalnumeric.LoadFile(raster)

# Also load as a gdal image to get geotransform 
# (world file) info
srcImage = gdal.Open(raster)
nrow = srcImage.RasterYSize
ncol = srcImage.RasterXSize
band = srcImage.GetRasterBand(1)
nodata = band.GetNoDataValue()
#print nrow, ncol, nodata
geoTrans = srcImage.GetGeoTransform()
#print geoTrans

clip = srcArray[0:nrow, 0:ncol]

# Create an OGR layer from a boundary shapefile
shapef = ogr.Open(shp)
lyr = shapef.GetLayer(0)
poly = lyr.GetNextFeature()

# Convert the layer extent to image pixel coordinates
minX, maxX, minY, maxY = lyr.GetExtent()
minx = geoTrans[0]
maxx = geoTrans[0] + ncol * geoTrans[1]
miny = geoTrans[3] + nrow * geoTrans[5]
maxy = geoTrans[3]
if not(maxY < miny or maxX < minx or minY > maxy or minX > maxx):
    #    sys.exit(1)
    #print minX, maxX, minY, maxY
    #print minx, maxx, miny, maxy
    ulX, ulY = world2Pixel(geoTrans, minX, maxY)
    lrX, lrY = world2Pixel(geoTrans, maxX, minY)
    #print ulX, lrX, ulY, lrY

    # Calculate the pixel size of the new image
    #pxWidth = int(lrX - ulX)
    #pxHeight = int(lrY - ulY)
    pxWidth = ncol
    pxHeight = nrow
    #print pxWidth, pxHeight

    # Create a new geomatrix for the image
    geoTrans = list(geoTrans)
    #print geoTrans

    # Map points to pixels for drawing the 
    # boundary on a blank black and white, mask image.
    mask = zeros((nrow, ncol))
    while poly:
        points = []
        pixels = []
        
        geom = poly.GetGeometryRef()
        pts = geom.GetGeometryRef(0)
        numpts = pts.GetPointCount()
        for p in range(numpts):
          points.append((pts.GetX(p), pts.GetY(p)))
        
        for p in points:
            pixels.append(world2Pixel(geoTrans, p[0], p[1]))

        rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
        rasterize = ImageDraw.Draw(rasterPoly)
        rasterize.polygon(pixels, 0)
        mask_tmp = imageToArray(rasterPoly)

        mask = mask + mask_tmp
        poly = lyr.GetNextFeature()
           
    mask = where(mask < 12, 1, 0)          

    # Clip the raster using the mask
    clip = gdalnumeric.choose(mask,(clip, nodata)).astype(gdalnumeric.float32)

driver = gdal.GetDriverByName('GTiff')
out_ds = driver.Create(output, ncol, nrow, 1,
gdal.GDT_Float32)
out_ds.GetRasterBand(1).WriteArray(clip,0,0)
out_ds.GetRasterBand(1).SetNoDataValue(nodata)
out_ds.SetGeoTransform(geoTrans)
out_ds.SetProjection(srcImage.GetProjection())
out_ds.FlushCache()
