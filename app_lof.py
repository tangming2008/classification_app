import streamlit as st
from matplotlib.pylab import *
import numpy as np
import pandas as pd


st.title('Predicting Lack-of-Fusion Porosity in Additive Manufacturing')
st.write('Authors: Ming Tang and Chris Pistorius, Carnegie Mellon University')
# st.markdown("Contact: mingt@alumni.cmu.edu, pcp@andrew.cmu.edu")

# Create two columns: 1) input, processing conditions, 2) output, lof fraction and meso-structure
col1, col2 = st.columns([1, 2])

with col1:
    # Process conditions - sliders in the first column
    widthAv = st.slider('Melt-Pool Width', min_value=0, max_value=400, step=10, value=200)/2 # half width
    depthAv = st.slider('Melt-Pool Depth', min_value=0, max_value=200, step=10, value=100)/2 # half depth
    hatch = st.slider('Hatch Spacing', min_value=0, max_value=200, step=10, value=100)
    layer = st.slider('Layer Thickness', min_value=0, max_value=100, step=10, value=50)
    angleradians = st.slider('Hatch Rotation Angle', min_value=0, max_value=90, step=5, value=45)

# Rest of your processing conditions
widthdev = 0.001
dwratioAv = depthAv/widthAv
dwdev = 0.001
hwratioAv = 0.001
hwdev = 0.001
anglestep = radians(angleradians)
angleoffset = radians(15.1)
minplot = 0
maxplot = 3.5

dims = 1001
xwidth = 1000.0
xpixel = xwidth / (dims-1)
dimsplus = dims + int(round((hwratioAv+dwratioAv)*4*widthAv/xpixel,0))
offset = int((dimsplus-dims)/2)

exptAv = np.pi * (widthAv**2) * (hwratioAv + dwratioAv) / (hatch*layer)

originx = 0.0
originz = -hwratioAv*widthAv*2
layers = int((xwidth-originz+(hwratioAv+dwratioAv)*2*widthAv)/layer)+3

a = zeros((dimsplus,dims))
pool = zeros((360,2))
collo_u = zeros(dimsplus, dtype=int)
collo_l = zeros(dimsplus, dtype=int)
colhi_u = zeros(dimsplus, dtype=int)
colhi_l = zeros(dimsplus, dtype=int)

b = zeros((dims,dims//10))
df_all = pd.DataFrame()

for layercount in range(layers):
    angle = anglestep * (layercount-3) + angleoffset
    xstep = fabs(hatch / cos(angle))
    originx = -1.0 * np.random.random() * hatch/fabs(cos(angle))
    row = dimsplus-int(round(originz/xpixel,0))
    hatchcount = int((xwidth+originx)/hatch*fabs(cos(angle))+ 10)

    for k in range(hatchcount):
        column = int(round(originx/xpixel,0))
        hwidth = np.random.normal(widthAv,widthdev*widthAv)
        dwratio = np.random.normal(dwratioAv,dwdev*dwratioAv)
        hwratio = np.random.normal(hwratioAv,hwdev*hwratioAv)
        xradius = fabs(hwidth / cos(angle))
        depth = dwratio * hwidth * 2
        height = hwratio * hwidth * 2

        pixelcount = 0.0
        poolaverage = 0.0

        indexr = int(hwidth*2.0*dwratio/xpixel)
        rowmin_l = max(0,row)
        rowmax_l = min(row+indexr,dimsplus)
        for i in range(rowmin_l,rowmax_l):
            zradius = xradius * sqrt(1 - ((i-row)*xpixel/(hwidth*2.0*dwratio))**2)
            locindexr = int(zradius/xpixel)
            collo_l[i] = max(-locindexr+column,0)
            colhi_l[i] = min(locindexr+1+column,dims)
            for j in range(collo_l[i],colhi_l[i]):
                pixelcount += 1
                poolaverage += a[i,j]

        indexr = int(hwidth*2.0*hwratio/xpixel)
        rowmin_u = min(row,dimsplus)
        rowmax_u = max(row-indexr,0)
        for i in range(rowmax_u,rowmin_u):
            zradius = xradius * sqrt(1 - ((i-row)*xpixel/(hwidth*2.0*hwratio))**2)
            locindexr = int(zradius/xpixel)
            collo_u[i] = max(-locindexr+column,0)
            colhi_u[i] = min(locindexr+1+column,dims)
            for j in range(collo_u[i],colhi_u[i]):
                pixelcount += 1
                poolaverage += a[i,j]

        pixelcount = max(1,pixelcount)
        poolaverage = poolaverage / pixelcount

        for i in range(rowmin_l,rowmax_l):
            for j in range(collo_l[i],colhi_l[i]):
                a[i,j] = poolaverage + 1

        for i in range(rowmax_u,rowmin_u):
            for j in range(collo_u[i],colhi_u[i]):
                a[i,j] = poolaverage + 1

        for i in range(180):
            calcAngle = radians(1.0*i)
            radius1 = xradius*height/sqrt((height*cos(calcAngle))**2 + (xradius*sin(calcAngle))**2)
            radius2 = xradius*depth/sqrt((depth*cos(calcAngle))**2 + (xradius*sin(calcAngle))**2)
            pool[i,0] = originx+radius1*cos(calcAngle)
            pool[i+180,0] = originx-radius2*cos(calcAngle)
            pool[i,1] = originz+radius1*sin(calcAngle)-(offset-1)*xpixel
            pool[i+180,1] = originz-radius2*sin(calcAngle)-(offset-1)*xpixel

        df = pd.DataFrame(pool)
        df.columns = ['x', 'y']
        df['layer'] = layercount
        df['hatch'] = k
        df_all = pd.concat([df_all, df])
        originx += xstep
    originz += layer

a = a[offset:offset+dims,0:dims]
allAv = np.average(a)
allMin = np.min(a)
if allMin == 0:
    allMin = allAv
    for i in range(dims):
        for j in range(dims):
            if (a[i,j]>0) and (a[i,j]<allMin):
                allMin = a[i,j]

df = df_all
fig, ax = plt.subplots(figsize=(3,3))
col_groupby = ['layer', 'hatch']
for index, (name, group) in enumerate(df.groupby(col_groupby)):
    x, y = group['x'], group['y']
    ax.fill(x, y, facecolor='grey', edgecolor='black', linewidth=0.5)
    ax.set_xlim(left=0, right=1000)
    ax.set_ylim(bottom=0, top=1000)
    ax.set_aspect('equal', adjustable='box')

# Display the figure in the second column
with col2:
    countsdist = np.histogram(a, bins=20, range=(0,5))
    totalcounts = np.sum(countsdist[0])
    st.write("""Total fraction of lack-of-fusion porosity:""", countsdist[0][0]*1.0/totalcounts)
    st.pyplot(fig)
