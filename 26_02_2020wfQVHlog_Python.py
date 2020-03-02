# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:13:26 2020

@author: burghelea
"""

#!/usr/bin/env python
# coding: utf-8

# In[17]:


import os.path
import os
import numpy as np
import pydicom
import matplotlib.pyplot as plt
import csv
import math

from tkinter import filedialog

#DVH calculation
import matplotlib.path
from six import iteritems


# In[18]:


#select the directoy and import the MAA and Y90 RTdose files

dirname = filedialog.askdirectory(initialdir=os.getcwd(),title='Please select directory of the pacient')
print (dirname)
RefDsMAA = filedialog.askopenfilename(initialdir = "dirname",title = "Select MAA DICOM File")
RefDsMAA = pydicom.read_file(RefDsMAA)
imgMAA=RefDsMAA.pixel_array
imgMAA.setflags(write=1)
ScaleMAA=RefDsMAA.DoseGridScaling

DsY90 = filedialog.askopenfilename(initialdir = "dirname",title = "Select Y90 DICOM File")
DsY90 = pydicom.read_file(DsY90)
imgY90=DsY90.pixel_array
imgY90.setflags(write=1)
ScaleY90=DsY90.DoseGridScaling

#plot the same slice for the MAA and Y90
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.imshow(imgMAA[45,:,:],cmap='binary') #affiche l'image correspondant à la coupe 26
ax1.set_title('MAA Dose')
ax2.imshow(imgY90[45,:,:],cmap='binary') #affiche l'image correspondant à la coupe 26
ax2.set_title('Y90 Dose')

#Verify if the two images have the same dimensions 
#print('MAA:',imgMAA.shape, ScaleMAA)
#print('Y90:', imgY90.shape, ScaleY90)
#print('MAA max dose Gy', imgMAA.max()*ScaleMAA)
#print('Y90 max dose Gy', imgY90.max()*ScaleY90)


# In[19]:


#open the directory of the masks and print the number of files to verify if the nr of lesions is correct
dirMask = filedialog.askdirectory(title='Please select dir of the masks')
nr_masks = sum(len(files) for _, _, files in os.walk(dirMask))
print("nr files:", nr_masks)

# list of the mask_names from dir without the extension .dcm
l=os.listdir(dirMask)
mask_names=[mask.split('.')[0] for mask in os.listdir(dirMask)] 
print(mask_names)

#create a list with all the dicom masks, this will be used to for loop through each mask for the QVH calc
masks = []
for i in range(nr_masks):
    masks.append(pydicom.read_file(os.path.normpath(os.path.join(dirMask,l[i]))))
print('shape:',masks[0].pixel_array.shape)


# In[20]:


class MaskResults:
   
    def __init__(self, imgMAA, ScaleMAA, imgY90, ScaleY90, mask, mask_name): #creating a class using as arguments all the paramaters requiered to calculate the QVH
        self.mask = mask
        self.mask_name = mask_name
        self.Qarraysorted, self.sumWF = self.qvh_calculation (imgMAA, ScaleMAA, imgY90, ScaleY90, mask) #the first and second outputs of the qvh_calculation function are saved
        self.x = [q[0] for q in self.Qarraysorted]
        self.y = [100]
        for i in range (1,len(self.Qarraysorted)):       # %volume of a voxel that is ponderated by the WF
            self.y.append(self.y[i-1]-self.Qarraysorted[i][1]/self.sumWF*100)
    
#to extand the lenght of each row based on the largest lenth of data # necessary to be able to export all the data in the same csv file
    def extend_xy(self, length):  
        x_len = len(self.x)
        num_empty_needed = length - x_len
        self.x.extend(["" for i in range(num_empty_needed)])#create an emplty string and add it to the list
        self.y.extend(["" for i in range(num_empty_needed)])
    
    
    def qvh_calculation(self,imgMAA, ScaleMAA, imgY90, ScaleY90, mask):
 # create a initial matrix for each factor
        wf_MAA=np.zeros_like(np.copy(imgMAA),dtype=float)
        wf_Y90=np.zeros_like(np.copy(imgMAA),dtype=float)
        wf_QVH=np.zeros_like(np.copy(imgMAA),dtype=float)
        Qarray = list()

        sumWF = 0
        qFWeightedSum = 0

    
# loop through each voxel of the image 
        for i in (range((imgMAA.shape[0]))):
            for j in (range((imgMAA.shape[1]))):
                for k in (range ((imgMAA.shape[2]))):         
                    if mask.pixel_array[i][j][k] == 65535:  #Calculate the quality division factor (qF) for each voxel inside the preselected mask
                        Dt = 0.1 if (imgMAA[i,j,k]*ScaleMAA)==0 else imgMAA[i,j,k]*ScaleMAA
                        Ds = 0.1 if (imgY90[i,j,k]*ScaleY90)==0 else imgY90[i,j,k]*ScaleY90
                        
                        qF = math.log10(Dt/Ds)
                                                         

                    #Calculate the weighting factor for MAA(wf_MAA) based on the predefined dose thresholds
                        if (imgMAA[i,j,k]*ScaleMAA)<=10 or (imgMAA[i,j,k]*ScaleMAA)>=200: 
                            wf_MAA=0
                        elif (imgMAA[i,j,k]*ScaleMAA)>10 and (imgMAA[i,j,k]*ScaleMAA)<40:
                            wf_MAA=((1/30)*(imgMAA[i,j,k]*ScaleMAA))-(1/3)
                        elif (imgMAA[i,j,k]*ScaleMAA)>120 and (imgMAA[i,j,k]*ScaleMAA)<200:
                            wf_MAA=((-2/160)*(imgMAA[i,j,k]*ScaleMAA))+(20/8) 
                        else:
                            wf_MAA=1

                    #Calculate the weighting factor for Y90(wf_Y90) based on the predefined dose thresholds
                        if (imgY90[i,j,k]*ScaleY90)<=10 or (imgY90[i,j,k]*ScaleY90)>=200:
                            wf_Y90=0
                        elif (imgY90[i,j,k]*ScaleY90)>10 and (imgY90[i,j,k]*ScaleY90)<40:
                            wf_Y90=((1/30)*(imgY90[i,j,k]*ScaleY90))-(1/3) 
                        elif (imgY90[i,j,k]*ScaleY90)>120 and (imgY90[i,j,k]*ScaleY90)<200:
                            wf_Y90=((-2/160)*(imgY90[i,j,k]*ScaleY90))+(20/8)  
                        else:
                            wf_Y90=1

                     #calculate the max wf between Y90 and MAA (wf_QVH)
                        if (imgMAA[i,j,k]*ScaleMAA)<80 and (imgY90[i,j,k]*ScaleY90)>80:
                            wf_QVH = 1
                        elif (imgY90[i,j,k]*ScaleY90)<80 and (imgMAA[i,j,k]*ScaleMAA)>80:
                            wf_QVH = 1
                        else:   
                            wf_QVH=max(wf_Y90,wf_MAA) 

                        sumWF += wf_QVH
                        qFWeightedSum += abs(qF)*wf_QVH #at the level of the voxel, the quality factor is multiplied with the coresponding wighting factor

                        Qarray.append((qF,wf_QVH))  # crate a 2D array:first colum result of the dose division, second column the weighting factors

        Qarraysorted = sorted(Qarray) #sort the data based on the dose division (lowest to highest)

#calculate QVH quality factor
        qualityFactor= qFWeightedSum/sumWF
        print(self.mask_name,'qualityFactor:', qualityFactor)

        return Qarraysorted, sumWF
    
    
#create a dictionary with all the masks    
    mask_dict = dict() 
    
#create a list with the outputs of the class MaskResults, for each mask     
res_list = [MaskResults(imgMAA, ScaleMAA, imgY90, ScaleY90, masks[i], mask_names[i]) for i in range(nr_masks)] 


# In[23]:



#plot QVH for each mask and save it in the pacient directory
fig = plt.figure(figsize = (20, 10))
print(nr_masks)
for m_res_index in range(nr_masks):
    plt.subplot(nr_masks, nr_masks, m_res_index+1)
    plt.plot(res_list[m_res_index].x,res_list[m_res_index].y)
    plt.xlabel('qF ')
    plt.ylabel('Volume (%)')
    plt.xlim(0,5)
    plt.ylim(0,100)
    plt.grid()
    plt.title (res_list[m_res_index].mask_name)

FigFile=os.path.join(os.path.realpath(dirname),"Patient.png")
fig.savefig(FigFile)

plt.show()


#create a csv file that includes the dose pondereted and volume for each mask
Xs = [res.x for res in res_list]
Ys = [res.y for res in res_list]
max_length = max([len(x) for x in Xs])

print(max_length)

for res in res_list:
    res.extend_xy(max_length)
    
Xs = [res.x for res in res_list]
Ys = [res.y for res in res_list]

#y=ponderated volume
#x= sorted dose division

# delim = ";"
# colnames = ["%s%s%s" % ("Dose_"+mask_names[i], delim,
#                           "Volume_"+mask_names[i])
#            for i in range(nr_masks)]

# with open(os.path.join(dirname, "Patient"+'.csv'), "w+") as csv_file:
#     csv_file.write(delim.join(colnames)+"\n")
    
#     for i in range(max_length):
#         for mask_idx in range(nr_masks):
#             csv_file.write(str(res_list[mask_idx].x[i]) + delim + str(res_list[mask_idx].y[i]) + delim)
#         csv_file.write("\n")


# In[ ]:





# In[ ]:




