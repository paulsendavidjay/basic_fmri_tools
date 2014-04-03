#! /usr/bin/env python
''' RUN FROM WITHIN SCRIPTS FOLDER
This script will collect the time course intensity from ventricle voxels specified in '../Analysis/ventricle_voxels_L.txt'
and '../Analysis/ventricle_voxels_R.txt' for the subject and run numbers passed in from the command line. 
Output will be single column txt files placed in the single_column_files in the Analysis directory, will be used at first level. 
fMRI data are from filtered_func in 0_prestats.
'''
import sys, os, glob, shutil, string, nibabel, subprocess
from numpy import *
from optparse import OptionParser


# SET UP PARSER FOR PASSING IN SOURCE AND DESTINATION DIRECTORY OPTIONS
parser = OptionParser()
parser.add_option("-s", "--subject", dest="subject")
parser.add_option("-y", "--year", dest="year")
parser.add_option("-r", "--run", dest="run")
(options, args) = parser.parse_args()

if options.run==None:
	raise NameError, "Must specify run: -r #"
if options.subject==None:
	raise NameError, "Must specify subject: -s #"
if options.year==None:
	raise NameError, "Must specify year: -y #"


if (options.subject == "subjectID"):
	sys.exit()

subject = options.subject
subject = "%03d" % (int(subject),)
year = options.year
run = options.run


# how many voxels in each direction from the voxel coordinate
clusterSize=1

# Define ventricle masks
left_ventricle_mask = '/Users/ncanda/Documents/Research/NCANDA/ROI_masks/ventricle_masks/ventricle_mask_L.nii.gz'
right_ventricle_mask = '/Users/ncanda/Documents/Research/NCANDA/ROI_masks/ventricle_masks/ventricle_mask_R.nii.gz'

data_dir = "/Users/ncanda/Documents/Research/NCANDA/data_MR/"
run_dir = data_dir + "A" + subject + "_" + year + "/run" + run

filteredFuncData = run_dir + "/nfswkmtd_" + subject + "_run" + run + ".nii.gz"
meanFilteredFuncData = run_dir + "/mean_nfswkmtd_" + subject + "_run" + run + ".nii.gz"

# Create mean of functional data
os.system(('fslmaths %s -Tmean %s' % (filteredFuncData, meanFilteredFuncData)) )


#####		RIGHT			######
xyzR = subprocess.check_output( ('fslstats %s -k %s -x' % (filteredFuncData, right_ventricle_mask)) , shell=True).strip().split(' ')
xyzR = [int(x) for x in xyzR]
[xR, yR, zR] = xyzR

sys.stdout.write(("Subject:" + subject + " Year:" + str(year) + " Run:" + run))
sys.stdout.write(('   Right:' + str(xR) + "," + str(yR) + "," + str(zR)))


#####		LEFT			######
xyzL = subprocess.check_output( ('fslstats %s -k %s -x' % (filteredFuncData, left_ventricle_mask)) , shell=True).strip().split(' ')
xyzL = [ int(x) for x in xyzL]
[xL, yL, zL] = xyzL

sys.stdout.write(('   Left:' + str(xL) + "," + str(yL) + "," + str(zL)))
sys.stdout.write('\n')

if xL == '': # if x was not set, then corrdinates not found, exit
	sys.stdout.write(("\tNo Left xyz coordinates found for %s run %s\n" % (subject, run)))
	sys.exit()


# WRITE VOXEL OUTPUT TO TEXT FILE
voxel_tc_output_file_name = run_dir + "/ventrical_coords.txt" # output file name
coord_out = open(voxel_tc_output_file_name, 'w')
coord_out.write( ('Right: %s %s %s\nLeft: %s %s %s' % (xR, yR, zR, xL, yL, zL) ))
coord_out.close()


# LOAD IMAGE
image_header = nibabel.load(filteredFuncData)
image_data = image_header.get_data()

# GET TIME COURSE OF VOXEL RIGHT
voxel_tc_a = image_data[xR-clusterSize:xR+clusterSize+1, yR-clusterSize:yR+clusterSize+1,zR-clusterSize:zR+clusterSize+1, :].mean(axis=1)
voxel_tc_b = voxel_tc_a.mean(axis=1)
voxel_tc_R = voxel_tc_b.mean(axis=0)

# GET TIME COURSE OF VOXEL LEFT
#print 'voxel_tc_a'
voxel_tc_a = image_data[xL-clusterSize:xL+clusterSize+1, yL-clusterSize:yL+clusterSize+1, zL-clusterSize:zL+clusterSize+1, :].mean(axis=1)
voxel_tc_b = voxel_tc_a.mean(axis=1)
voxel_tc_L = voxel_tc_b.mean(axis=0)


# WRITE TIME COURSE TO OUTPUT FILE
voxel_tc_output_file_name = run_dir + "/ventrical_tc.txt" # output file name
voxel_tc_output = open(voxel_tc_output_file_name, 'w')
for i in range(0, len(voxel_tc_L)):
	output_line = str(voxel_tc_L[i]) + " " + str(voxel_tc_R[i]) + "\n"
	voxel_tc_output.write(output_line)

voxel_tc_output.close()


