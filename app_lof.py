import streamlit as st
from matplotlib.pylab import *
import numpy as np

# Process conditions

widthAv = st.slider('Ave melt-pool width', min_value = 0, max_value = 100, step=10, value = 50)
# widthAv = 60.0   # average melt pool half width
widthdev = 0.001    # relative standard deviation of width
dwratioAv = 0.0   # average (D-R)/W ratio of melt pool
dwdev = 0.001       # relative standard deviation of depth:width ratio
hwratioAv = 0.000   # average R/W ratio of melt pool
hwdev = 0.001       # relative standard deviation of height:width ratio
hatch = 110.0     # hatch spacing
layer = 40.0      # layer thickness
anglestep = radians(45.0)  # angle between rasters in subsequent layers
angleoffset = radians(15.1) # To ensure that cos(angle) is never zero
minplot = 0       # minimum and maximum
maxplot = 3.5       # values for plot range of number of melting cycles

# Definition of calculation domain
dims = 501        # lateral dimension of square calculation domain in pixels
xwidth=500      # lateral dimension of square calculation domain in microns

xpixel = xwidth / (dims-1) # pixel size in microns
dimsplus = dims + int(round((hwratioAv+dwratioAv)*4*widthAv/xpixel,0))
offset = int((dimsplus-dims)/2)

# Expected average melt count:
exptAv = np.pi * (widthAv**2) * (hwratioAv + dwratioAv) / (hatch*layer)

# Definition of position of first pass
originx = 0.0
originz = -hwratioAv*widthAv*2
layers = int((xwidth-originz+(hwratioAv+dwratioAv)*2*widthAv)/layer)+3

a = zeros((dimsplus,dims)) # zero all pixels to begin with
pool=zeros((360,2))  # array which contains coordinates of melt pool boundary
collo_u=zeros(dimsplus,dtype=np.int)    # array with indices of first column of melt pool in a given row
collo_l=zeros(dimsplus,dtype=np.int)    # array with indices of first column of melt pool in a given row
colhi_u=zeros(dimsplus,dtype=np.int)    # as collo, but last column
colhi_l=zeros(dimsplus,dtype=np.int)    # as collo, but last column

# matrix for drawing color scale
b = zeros((dims,dims//10))

fig = plt.figure(figsize=(8, 6))


figure(1,facecolor="white")   # set up window for plotting melt pool shapes
axes()

for layercount in range(layers):  # iterate over layers
   angle = anglestep * (layercount-3) + angleoffset
   xstep = fabs (hatch / cos(angle))
   originx = -1.0 * np.random.random() * hatch/fabs(cos(angle))  # randomize position of leftmost melt pool
   # originx = (0.0 * fabs (hatch / cos(angle)) + xwidth) * iseven
   row = dimsplus-int(round(originz/xpixel,0)) # position of top of current layer in pixels
   # print row
   hatchcount = int((xwidth+originx)/hatch*fabs(cos(angle))+ 10) # number of scans across section in current layer
   for k in range(hatchcount): # iterate over scans
      column = int(round(originx/xpixel,0))
      hwidth = np.random.normal(widthAv,widthdev*widthAv) # randomize pool half width
      dwratio = np.random.normal(dwratioAv,dwdev*dwratioAv) # randomize pool depth:width ratio
      hwratio = np.random.normal(hwratioAv,hwdev*hwratioAv) # randomize pool height:width ratio
      xradius = fabs (hwidth / cos(angle))   # horizontal half-width of melt pool on section, in microns
      depth = dwratio * hwidth * 2
      height = hwratio * hwidth * 2

      # Find the average melt count in the region to be covered by the pool was before the current melting event
      pixelcount = 0.0
      poolaverage = 0.0

      # First define the lower part of melt pool
      indexr = int(hwidth*2.0*dwratio/xpixel) # melt pool depth in pixels
      rowmin_l=max(0,row)
      rowmax_l=min(row+indexr,dimsplus)  # Deleted +1
      for i in range(rowmin_l,rowmax_l):
         zradius = xradius * sqrt(1 - ((i-row)*xpixel/(hwidth*2.0*dwratio))**2)
         locindexr = int(zradius/xpixel)
         collo_l[i] = max(-locindexr+column,0)
         colhi_l[i]  = min(locindexr+1+column,dims)
         for j in range(collo_l[i],colhi_l[i]):
            pixelcount += 1
            poolaverage += a[i,j]

      # Now define the upper part of melt pool ("reinforcement)
      indexr = int(hwidth*2.0*hwratio/xpixel) # melt pool height in pixels
      rowmin_u=min(row,dimsplus)
      rowmax_u=max(row-indexr,0)
      # Find the average melt count of the reinforcement (should be zero)
      for i in range(rowmax_u,rowmin_u):
         zradius = xradius * sqrt(1 - ((i-row)*xpixel/(hwidth*2.0*hwratio))**2)
         locindexr = int(zradius/xpixel)
         collo_u[i] = max(-locindexr+column,0)
         colhi_u[i]  = min(locindexr+1+column,dims)
         for j in range(collo_u[i],colhi_u[i]):
            pixelcount += 1
            poolaverage += a[i,j]

      pixelcount = max(1,pixelcount)
      poolaverage = poolaverage / pixelcount

      # Now increase melt count by 1
      for i in range(rowmin_l,rowmax_l):
         for j in range(collo_l[i],colhi_l[i]):
            a[i,j]=poolaverage + 1

      # Now increase melt count by 1
      for i in range(rowmax_u,rowmin_u):
         for j in range(collo_u[i],colhi_u[i]):
            a[i,j]=poolaverage + 1

      # Calculate melt pool boundary
      for i in range(180):
            calcAngle = radians(1.0*i)
            radius1 = xradius*height/sqrt((height*cos(calcAngle))**2 + (xradius*sin(calcAngle))**2)
            radius2 = xradius*depth/sqrt((depth*cos(calcAngle))**2 + (xradius*sin(calcAngle))**2)
            pool[i,0]=originx+radius1*cos(calcAngle)
            pool[i+180,0]=originx-radius2*cos(calcAngle)
            pool[i,1]=originz+radius1*sin(calcAngle)-(offset-1)*xpixel
            pool[i+180,1]=originz-radius2*sin(calcAngle)-(offset-1)*xpixel
      poly1=Polygon(pool,facecolor='white',edgecolor='black', linewidth=0.5)
      gca().add_patch(poly1)  # add melt pool shape to plot

      # originx += fabs (hatch / cos(angle))
      originx += xstep
   originz += layer

#remove transition layers at top and bottom
a=a[offset:offset+dims,0:dims]
allAv = np.average(a)
allMin = np.min(a)
if allMin == 0:
   allMin = allAv
   for i in range(dims):
      for j in range(dims):
         if (a[i,j]>0) and (a[i,j]<allMin):
            allMin = a[i,j]

# set up plot of melt pool outlines
axis('scaled')             # ensures correct aspect ratio
axis([0,xwidth,0,xwidth])  # ranges for melt pool shape plot
axis('off')                # do not show axes

#print counts/(dims*dims)
"average count, expected av, max count, min excluding zero, min"             # average melt count
round(allAv,2), round(exptAv,2), round(np.max(a),2), round(allMin,2), round(np.min(a),2)
countsdist = np.histogram(a,bins=20,range=(0,5))
totalcounts = np.sum(countsdist[0])
# "Distribution of number of melting cycles:", countsdist[0]*1.0/totalcounts
# "Bins used to determine distribution:", countsdist[1]

# Plot map of number of melting cycles
# matshow(a,cmap=cm.afmhot, vmin=minplot, vmax=maxplot) # other possible colormaps: afmhot, coolwarm, bone, gray, RdGy, hot, afmhot

# Calculate and draw color scale
# for i in range(dims-2):
#    for j in range(dims//10-2):
#       b[i,j]=i*1.0/(dims-1)*(maxplot - minplot) + minplot
#    for j in (0,1,dims//10-2,dims//10-1):
#       b[i,j]=minplot
# for i in (dims-3,dims-2,dims-1):
#    for j in range(dims//10):
#       b[i,j]=minplot

# matshow(b,cmap=cm.afmhot, vmin=minplot, vmax=maxplot, origin="lower") # other possible colormaps: afmhot, coolwarm, bone, gray, RdGy, hot, afmhot
# axis('off')
# scalerange = maxplot-minplot
# for i in range(int(scalerange+1)):
#    text(-10,int(round(i/scalerange*(dims-1),0)),str(round(minplot+i,1)),horizontalalignment='right',verticalalignment='center',fontsize=40)
# text(-10,(dims-1),str(round(maxplot,1)),horizontalalignment='right',verticalalignment='center',fontsize=40)

show()                     # show all figures
st.pyplot(fig)
