#Image analysis notes 20190109 Karin
#
#NOTE: in Jython interpreter you can get the ROI via proi = IJ.getImage().getRoi()
#
#Working protocol:
#step 1: Take slice of zstack you want to analyse with ECM labelled and imaged
#
#step 2: draw ROIs rectangles of fixed width of the ECM you want to analyse. Aim for 15um thick rectangles with 150um long. Save the rois.
#
#step 3: isolate ECM channel, select ROI you want to start with, Edit-> selection -> Straighten. Now it will be flat on one line. 
# Double check that the organoid should be located on the left via flipping horizontally if needed.
# Save the straightened ROI
# 
#step 4: analyze the straightened ROI with orientationJ, (parameters depend on the image), gaussian window
# save the orientation file and the coherency file (make sure these images are shown)
#
#step5: Run the following code in steps in Jython interpreter:

from __future__ import with_statement, division
import math
from ij import IJ, ImagePlus
from ij.process import FloatProcessor

#remove the old date file if it is there
import os
try:
  os.remove('ECM_results.csv')
except:
  pass


#function get_order_param to calculate order parameter, using a threshold for the coherency
#to reduce noise. This function also gives the orientation distribution and average orientation

def get_order_param_v2(orientation,coherency,threshold): #Or=orientation image, Coh= coherency image
	import math
	from ij import IJ, ImagePlus
	from ij.process import FloatProcessor

	orient_array=orientation.getProcessor().getPixels()
	coher_array=coherency.getProcessor().getPixels()

	mask=map(lambda x: x>=set_coh, coher_array)

	sin_or=[]
	cos_or=[]
	weight_sin_or=[]
	weight_cos_or=[]
	coh_list=[]
	or_list=[]

	for i in range(len(mask)):
		if mask[i]==True:
			sin_tmp= math.sin(float(orient_array[i])*2) # here I do sin(2*angle) where angle is in radian
			cos_tmp=math.cos(float(orient_array[i]) *2)
			weight_sin_or_tmp=float(coher_array[i])*math.sin(float(orient_array[i])*2)
			weight_cos_or_tmp=float(coher_array[i])*math.cos(float(orient_array[i])*2)
			cos_or.append(cos_tmp)
			sin_or.append(sin_tmp)	
			weight_sin_or.append(weight_sin_or_tmp)
			weight_cos_or.append(weight_cos_or_tmp)
			coh_list.append(coher_array[i])
			or_list.append(orient_array[i]*360/(2*math.pi))

	avg_sin_or=float(sum(sin_or)) / float(len(sin_or))
	avg_cos_or=float(sum(cos_or))/float(len(cos_or))
	avg_weight_sin_or=float( sum(weight_sin_or))/ float(sum(coh_list))
	avg_weight_cos_or=float( sum(weight_cos_or))/ float(sum(coh_list))

	order_param=math.sqrt(avg_sin_or**2+avg_cos_or**2)
	weighted_order_param=math.sqrt(avg_weight_sin_or**2+avg_weight_cos_or**2)
	#print (order_param,weighted_order_param)	

	return(order_param,weighted_order_param,coh_list,or_list)


#############START WITH DETERMINING ORIENTATION,COHERENCY AND INTENSITY#############
####################################################################################
#######################################################################################
##########################################################################################
#SELECT COHERENCY IMAGE from OrientationJ in IMAGEJ####################################

Coh= IJ.getImage()

############################################################################################
#SELECT Orientation IMAGE from OrientationJ#########################

Or= IJ.getImage()


#############################################################################################
### Determine order parameter by 3 different regions: close to the organoid, middle and far away
#########################################################################################
#########################################################################################
###SELECT the ORIGINAL ROI (ECM signal), that is callibrated#########################
#######################################################################################

imp = IJ.getImage()
X=imp.getWidth()
Y=imp.getHeight()

Border1=50 #first border in microns
Border2=100 #second border in microns
Border3=150
##############################
###SELECT YOUR X0#############
###SELECT CUTOFF COHERENCY####
##############################

##for 20181031 H7 from 12-10 T37 +FGF2 plane8 ROI2.tif 
#x0=7.44 #in microns
##for 20181031 H7 from 12-10 T37 +FGF2 plane8 ROI1.tif 
#x0=10.24 #in microns
##for 20181031 H7 from 12-10 T37 +FGF2 plane6 ROI2.tif
#x0=4.58
##for 20181031 H7 from 12-10 T37 +FGF2 plane5 ROI3.tif
#x0=12.8
##for 20181031 H7 from 12-10 T37 +FGF2 plane5 ROI2.tif
#x0=11.67
##for 20181031 H7 from 12-10 T37 +FGF2 plane5 ROI1.tif
#x0= 4.76
##for 20181031 H7 from 12-10 T37 +FGF2 plane1 ROI2
#x0= 5.75
##for 20181031 H7 from 12-10 T37 -FGF2 plane1 ROI1
x0= 13.76

set_coh=0.4 #Value between 0 and 1

####################
#####################

cal = imp.getCalibration()
#calibration = [cal.pixelWidth, cal.pixelHeight, cal.pixelDepth]

Xrange=range(0, X)
Xrange= [x * cal.pixelWidth for x in Xrange] #now we have a range of x-values which are calibrated

indexB0=[ n for n,i in enumerate(Xrange) if i>x0 ][0] #Select the last index 
Xrange= [x-x0 for x in Xrange]#adjusted for x0 offset
try:   #since I do not know how big the original ROI is, first try to draw the roi
    indexB1=[ n for n,i in enumerate(Xrange) if i>Border1 ][0]-1 #Select the last index 
except:
    indexB1=len(Xrange)-1

try:
    indexB2=[ n for n,i in enumerate(Xrange) if i>Border2 ][0]-1 #Select the last index 
except:
    indexB2=len(Xrange)-1

try:
	indexB3=[ n for n,i in enumerate(Xrange) if i>Border3 ][0]-1 #Select the last index
except:
    indexB3=len(Xrange)-1


##select region of interests (roi)
rm = RoiManager.getInstance()
if not rm:
  rm = RoiManager()
rm.reset()
roi = Roi(indexB0,0,indexB1-indexB0,Y) #Roi(x0,y0,x_width,y_width)
rm.addRoi(roi)
roi = Roi(indexB1+1,0,indexB2-indexB1,Y)
rm.addRoi(roi)
roi = Roi(indexB2+1,0,indexB3-indexB2,Y)
rm.addRoi(roi)


#Cool, making automated rois work, now let's determine 
#1) average intensity
#2) the average orientation and the histogram of orientation per region.
#3) order parameter per region 
#for point 2 and 3 I want to collect all the pixel values that have coherency > set_coh, to limit noise

### point 1 (average intensity)#####
means=[]

for i in RoiManager.getInstance().getRoisAsArray():
  imp.setRoi(i)
  stats = imp.getStatistics(Measurements.MEAN)
  means.append(stats.mean)



##### point 2 and 3 combined
param=[]
param_weighted=[]

coh_roi_array=[] #these lists will contain arrays of the coherency and orientation of the roi's, where the coherency is above the set threshold
or_roi_array=[]
means_angle=[]

for k in RoiManager.getInstance().getRoisAsArray():

	Coh.setRoi(k)
	Or.setRoi(k)

	image_or=Or.crop()
	image_coh=Coh.crop()


	[order_param_tmp,order_param_weighted_tmp,coh_tmp,or_tmp]=get_order_param_v2(image_or,image_coh, set_coh)  
	param.append(order_param_tmp)
	param_weighted.append(order_param_weighted_tmp)
	coh_roi_array.append(coh_tmp) #coherency values are above the set threshold 
	or_roi_array.append(or_tmp)  #orientation values are angle values in degrees
	image_coh.close()
	image_or.close()
	avg_orientation_tmp=float(sum(or_tmp)) / float(len(or_tmp))
	means_angle.append(avg_orientation_tmp)
	

## now exporting the data to csv files


try: #if the file already exists, remove it to create a new one
	os.remove('order_param_results.csv')
except:
  	pass

g = open('order_param_results.csv', 'w+')
g.write('up to border nr#'+','+'param'+ ','+ 'param_weighted'+ ','+ 'mean intensity' + ','+ 'average angle' '\n')

index=0
for i in param:
	
	
	g.write(str(index)+','+str(i)+ ','+ str(param_weighted[index])+ ','+ str(means[index]) +',' + str(means_angle[index])+'\n')
	index=index+1 
g.close()	


try:  #if the file already exists, remove it to create a new one
	os.remove('Angle_distribution_3borders.csv')
except:
  	pass

g = open('Angle_distribution_3borders.csv', 'w+')
g.write('angle_border0' + ','+ 'coherency_border0' +','+'\n')



index=0
for i in coh_roi_array[0]:
	
	g.write(str(or_roi_array[0][index])+','+str(i)+'\n')

	index=index+1 

g.write(''+','+''+','+'\n')

g.write('angle_border1' + ','+ 'coherency_border1'+'\n')

index=0
for i in coh_roi_array[1]:
		
	g.write(str(or_roi_array[1][index])+','+str(i)+'\n')

	index=index+1 

g.write(''+','+''+','+'\n')

g.write('angle_border2' + ','+ 'coherency_border2'+'\n')

index=0
for i in coh_roi_array[2]:
		
	g.write(str(or_roi_array[2][index])+','+str(i)+'\n')

	index=index+1 


g.close()	












#############################################################################################
##Determine order parameter with a sliding window############################################
#
#rm_order = RoiManager.getInstance()
#if not rm_order:
 # rm_order = RoiManager()
#rm_order.reset()
#
#create Rois for coherency and orientation image
#for i in range(0, X-window):
#	roi_order = Roi(i,0,window,Y)
#	rm_order.addRoi(roi_order)
#
#rois=rm_order


# CHANGE WINDOW FOR DETERMINING THE ORDER PARAMETER
window=200 #this window is the width of the sliding window in pixels. Can be adjusted. Would recommend that it is bigger than twice the observed pore-size.


rm = RoiManager.getInstance()
if not rm:
  rm = RoiManager()
rm.reset()



for i in range(0, X-window):
	roi = Roi(i,0,window,Y)
	rm.addRoi(roi)
 

#function get_order_param to calculate order parameter
def get_order_param(orientation,coherency): #Or=orientation image, Coh= coherency image
	import math
	from ij import IJ, ImagePlus
	from ij.process import FloatProcessor

	orient_array=orientation.getProcessor().getPixels()
	coher_array=coherency.getProcessor().getPixels()
	
	sin_or=[]
	cos_or=[]
	weight_sin_or=[]
	weight_cos_or=[]
	coh_list=[]

	for i in range(len(orient_array)):
		sin_tmp= math.sin(float(orient_array[i])*2) # sin(2*angle) where angle is in radian
		cos_tmp=math.cos(float(orient_array[i]) *2)
		weight_sin_or_tmp=float(coher_array[i])*math.sin(float(orient_array[i])*2)
		weight_cos_or_tmp=float(coher_array[i])*math.cos(float(orient_array[i])*2)
		cos_or.append(cos_tmp)
		sin_or.append(sin_tmp)	
		weight_sin_or.append(weight_sin_or_tmp)
		weight_cos_or.append(weight_cos_or_tmp)
		coh_list.append(coher_array[i])

	avg_sin_or=float(sum(sin_or)) / float(len(sin_or))
	avg_cos_or=float(sum(cos_or))/float(len(cos_or))
	avg_weight_sin_or=float( sum(weight_sin_or))/ float(sum(coh_list))
	avg_weight_cos_or=float( sum(weight_cos_or))/ float(sum(coh_list))

	order_param=math.sqrt(avg_sin_or**2+avg_cos_or**2)
	weighted_order_param=math.sqrt(avg_weight_sin_or**2+avg_weight_cos_or**2)

	return(order_param,weighted_order_param)






param=[]
param_weighted=[]

for k in RoiManager.getInstance().getRoisAsArray():

	Coh.setRoi(k)
	Or.setRoi(k)
	#image_or=Or.getProcessor().duplicate  #here somehow the file type changes or something, so get_order_param does not work at getProcessor().getPixels
	#image_coh=Coh.getProcessor().duplicate
	image_or=Or.crop()
	image_coh=Coh.crop()
	[order_param_tmp,order_param_weighted_tmp]=get_order_param(image_or,image_coh)  
	param.append(order_param_tmp)
	param_weighted.append(order_param_weighted_tmp)
	image_coh.close()
	image_or.close()
	

os.remove('order_param_results.csv')
g = open('order_param_results.csv', 'w+')
g.write('pixel'+','+'um'+','+'param'+ ','+ 'param_weighted'+ '\n')

index=0
for i in param:
	
	
	g.write(str(index)+','+ str(index*cal.pixelWidth+window/2*cal.pixelWidth)+ ','+str(i)+ ','+ str(param_weighted[index])+'\n')
	index=index+1 
g.close()	


################################################################################################################################################################
################################################################################################################################################################
#Determine allignment on the whole image:

################################################################################################################################################################
################################################################################################################################################################
from __future__ import with_statement, division
import math
from ij import IJ, ImagePlus
from ij.process import FloatProcessor

or_array=Or.getProcessor().getPixels()
coh_array=Coh.getProcessor().getPixels()

sin_or=[]
cos_or=[]
weight_sin_or=[]
weight_cos_or=[]
coh_list=[]

for i in range(len(or_array)):
	sin_tmp= math.sin(float(or_array[i])*2) # sin(2*angle) where angle is in radian
	cos_tmp=math.cos(float(or_array[i]) *2)
	weight_sin_or_tmp=float(coh_array[i])*math.sin(float(or_array[i])*2)
	weight_cos_or_tmp=float(coh_array[i])*math.cos(float(or_array[i])*2)
	cos_or.append(cos_tmp)
	sin_or.append(sin_tmp)	
	weight_sin_or.append(weight_sin_or_tmp)
	weight_cos_or.append(weight_cos_or_tmp)
	coh_list.append(coh_array[i])

avg_sin_or=float(sum(sin_or)) / float(len(sin_or))
avg_cos_or=float(sum(cos_or))/float(len(cos_or))
avg_weight_sin_or=float( sum(weight_sin_or))/ float(sum(coh_list))
avg_weight_cos_or=float( sum(weight_cos_or))/ float(sum(coh_list))

order_param=math.sqrt(avg_sin_or**2+avg_cos_or**2)
weighted_order_param=math.sqrt(avg_weight_sin_or**2+avg_weight_cos_or**2)
print "order parameter= ", order_param, "weighted order parameter= ", weighted_order_param



--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------below is old code, redoing it now!!!!

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#now, determine the bit-depth of your image:

type=imp.getType
string_value = type.__str__()
str1="bit"
bit=string_value.find(str1)
#now select the bit depth
bit_depth=string_value[bit-3:bit-1]
bit_depth=int(bit_depth)




def split_list(alist, wanted_parts=1):
    """Split a list to the given number of parts."""
    length = len(alist)
    # alist[a:b:step] is used to get only a subsection of the list 'alist'.
    # alist[a:b] is the same as [a:b:1].
    # '//' is an integer division.
    # Without 'from __future__ import division' '/' would be an integer division.
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
	for i in range(wanted_parts) ]

matrix_coh = split_list(Coh.getProcessor().getPixels())
matrix_or = split_list(Or.getProcessor().getPixels())
sin_or=[]
matrix_or= [list(x) for x in zip(*matrix_or)]
for y in range(Y): #Y is the size in y
	for x in range(X): #X is the size in x
		
		# When assigning, we multiply the value by the amplitude.
		sin_tmp= math.sin(matrix_or[x][y])
		sin_or.append(sin_tmp)

in matlab:
coh=double(Image_coherency)./(2^16); 
or=double(Image_orientation)./(2^16);


order_parameter=( mean(  mean(sin(2.*(double(Image_orientation)./(2.^16)).*pi))).^2  +mean(mean(cos(2.*(double(Image_orientation)./(2.^16)).*pi))).^2 ).^(1/2)
%This one seems to work, now why doesn't the weighted one work?

order_parameter_weighted=sqrt((sum(sum(coh.*cos(2*or*pi)))/sum(sum(coh)))^2 + (sum(sum(coh.*sin(2*or*pi)))/sum(sum(coh)))^2)

