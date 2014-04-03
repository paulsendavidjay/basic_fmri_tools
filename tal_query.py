#! /usr/bin/env python
'''
tal_query.py
08-25-2009 
David Paulsen
Edited to include functional mask making 08-18-2010
Edited for nibabel 03-26-2013


note that output from fsl's cluster command is reverse coded for cluster number

This script is set up to 
	(1) query an FSL image using "cluster" with a specified threshold 
	(2) use the talairach demon to collect labels for cluster peaks.
	(3) Optionally:
		(i) create an image containing 3x3x3 cubes centered on each peak
		(ii) create a single image for each 3x3x3 cube centered on each peak
		(iii) create a functional mask for each cluster from which peaks are derrived
		-note that option iii may not work well for extremly large clusters of abnormal shape.
		 E.g. a cluster with a peak in occipital cortex but extends throughout the striatum and 
		 back along the cingulate. CHECK MASKS!!!

Example: 
	
	python tal_query.py fsl_output_file.nii.gz 3.2 [optional talairach search area parameter (see below)]
	

A path to talairach.jar file much be specified in the path check section. The jar file can be found at http://www.talairach.org/talairach.jar .

Output will be in tab-delimited format, for easy manipulation with Excel, should you direct the output to a file.

TAL can be set to False to turn off talairach querry and just use "cluster"
min_cluster_size can also be set to limit talairach querries to larger cluster sizes.

If you know how to set the search area for talairach, you may, as optional additional parameter ( http://www.talairach.org/manual.html )
default value is "3:7"
'''

import sys, os, string
from numpy import *
from datetime import date
from time import sleep
sys.setrecursionlimit(100000) #
if len(sys.argv) < 3:
	raise NameError, "must specify file and threshold"


input_file = str(sys.argv[1])
threshold = sys.argv[2]


# USR PARAMETERS
TAL = True # set to True to search for Talairach labels
min_cluster_size = 2 # for region reporting ; not statistical
BRAIN = False # set to True to create image with points at peak coordinates
MASKS = 0
# set to 0 to create 1 image with a 3x3 mask colored for every peak point, intensity value = peak voxel number (e.g. 1,2,3,...,n)
# set to 1 to create multiple images, each with 3x3 mask for every peak point, intensity = peak intensity
# CURRENTLY NOT WORKING set to 2 to create mutliple images, each with one mask per cluster point



#################################################
#				PATH CHECK
#################################################

# jar file is needed to get tailorach labels ( http://www.talairach.org/talairach.jar )
talairach_jar = '/opt/ni_tools/talairach.jar' 

# if you know how to set the search area for talairach, you may, as optional additional parameter ( http://www.talairach.org/manual.html )
if TAL is True:
	if len(sys.argv) == 4:
		tal_search_area = str(sys.argv[3])
	else:
		tal_search_area = "3:7"
	if not os.path.exists(input_file): # make sure data exists ( or that subject is entered correctly )
		print "\nfile not found: " + input_file + "\n"
		os.sys.exit(0)


if not os.path.exists(talairach_jar): # make sure data exists ( or that subject is entered correctly )
	print "\nfile not found: " + talairach_jar + "\n"
	print "Jar file can be found at http://www.talairach.org/talairach.jar\n"
	print "Download jar file and specify location in this script (tal_query.py)"
	os.sys.exit(0)


#################################################
#				FUNCTIONS
#################################################

# function for converting FSL_MNI to TAL ( http://brainmap.org/icbm2tal/index.html )
def _mni2tal_(x, y, z):
	incoord = array([x, y, z, 1])
	icbm_fsl = array([[0.9464, 0.0034, -0.0026, -1.0680],
		[-0.0083, 0.9479, -0.0580, -1.0239],
		[0.0053, 0.0617,  0.9010,  3.1883],
		[0.0000, 0.0000,  0.0000,  1.0000]])
	outmat = icbm_fsl*incoord
	tal_coord = []	
	for i in range(0,3):
		tal_coord.append(sum(outmat[i]))
	return tal_coord

def _tal2mni_(x, y, z):
	incoord = array([x, y, z, 1])
	icbm_fsl = array([[0.9464, 0.0034, -0.0026, -1.0680],
		[-0.0083, 0.9479, -0.0580, -1.0239],
		[0.0053, 0.0617,  0.9010,  3.1883],
		[0.0000, 0.0000,  0.0000,  1.0000]])
	outmat = linalg.inv(icbm_fsl)*incoord # use the inverse
	tal_coord = []	
	for i in range(0,3):
		tal_coord.append(sum(outmat[i]))
	return tal_coord





def _tal_names_(jar_file, search_area, x, y, z):
	java_cmd = 'java -cp ' + jar_file + ' org.talairach.PointToTD ' + search_area + ', ' + str(x) + ', ' + str(y) + ', ' + str(z)
	dummy, f = os.popen2(java_cmd) # f will become the output
	sleep(0.5)
	output = f.readlines()
 	if len(output) < 2:
		print "missed connection 1"
		dummy, f = os.popen2(java_cmd) # f will become the output
		sleep(0.5)		
		output = f.readlines()
		if len(output) < 2:
			print "missed connection 2"
			dummy, f = os.popen2(java_cmd) # f will become the output
			sleep(0.5)		
			output = f.readlines()
		elif output[1].find("TalairachUtility") == 0:
			print "missed connection 2"
			dummy, f = os.popen2(java_cmd) # f will become the output
			sleep(0.5)		
			output = f.readlines()
 	elif output[1].find("TalairachUtility") == 0:
 		print "missed connection 1"
		dummy, f = os.popen2(java_cmd) # f will become the output
		sleep(0.5)		
		output = f.readlines()
		if len(output) < 2:
			print "missed connection 2"
			dummy, f = os.popen2(java_cmd) # f will become the output
			sleep(0.5)		
			output = f.readlines()
		elif output[1].find("TalairachUtility") == 0:
			print "missed connection 2"
			dummy, f = os.popen2(java_cmd) # f will become the output
			sleep(0.5)		
			output = f.readlines()
	return output[2:]


''' get_functional_roi_steps & get_functional_roi work together '''

def get_functional_roi_steps(image_data, mask_data, z, y, x, iteration):
	mask_data[z,y, x] = 1
	if mask_data[z, y, x+1] == 0:
		if image_data[z, y, x+1] > 0:
			mask_data[z, y, x+1] = 1
			get_functional_roi_steps(image_data, mask_data, z, y, x+1, iteration)
	if iteration == 1:
		return mask_data
	if mask_data[z, y, x-1] == 0:
		if image_data[z, y, x-1] > 0:
			mask_data[z, y, x-1] = 1
			get_functional_roi_steps(image_data, mask_data, z, y, x-1, iteration)
	if iteration == 2:
		return mask_data
	if mask_data[z, y+1, x] == 0:
		if image_data[z, y+1, x] > 0:
			mask_data[z, y+1, x] = 1
			get_functional_roi_steps(image_data, mask_data, z, y+1, x, iteration)
			get_functional_roi_steps(image_data, mask_data, z, y+1, x, iteration)
	if iteration == 3:
		return mask_data
	if mask_data[z, y-1, x] == 0:
		if image_data[z, (y-1), x] > 0:
			mask_data[z, (y-1), x] = 1
			get_functional_roi_steps(image_data, mask_data, z, y-1, x, iteration)
			get_functional_roi_steps(image_data, mask_data, z, y-1, x, iteration)
	if iteration == 4:
		return mask_data
	if mask_data[z+1, y, x] == 0:
		if image_data[z+1, y, x] > 0:
			mask_data[z+1, y, x] = 1
			get_functional_roi_steps(image_data, mask_data, z+1, y, x, iteration)
			get_functional_roi_steps(image_data, mask_data, z+1, y, x, iteration)
	if iteration == 5:
		return mask_data	
	if mask_data[z-1, y, x] == 0:
		if image_data[z-1, (y), x] > 0:
			mask_data[z-1, (y), x] = 1
			get_functional_roi_steps(image_data, mask_data, z-1, y, x, iteration)
			get_functional_roi_steps(image_data, mask_data, z-1, y, x, iteration)
	return mask_data

def get_functional_roi(image_data, z, y, x):
	mask_data = zeros_like(image_data)
	mask_data[z,y, x] = 1
	mask_data = get_functional_roi_steps(image_data, mask_data, z, y, x, 1)
	mask_data = get_functional_roi_steps(image_data, mask_data, z, y, x, 2)
	mask_data = get_functional_roi_steps(image_data, mask_data, z, y, x, 3)
	mask_data = get_functional_roi_steps(image_data, mask_data, z, y, x, 4)
	mask_data = get_functional_roi_steps(image_data, mask_data, z, y, x, 5)
	mask_data = get_functional_roi_steps(image_data, mask_data, z, y, x, 6)
	return mask_data


#################################################
#				BEGIN CLUSTER SEARCH
#################################################

print '\ninput file =', input_file
print 'threshold =', threshold
print ""
print "\nnumber\tsize\tz_val\tx(vox)\ty(vox)\tz(vox)" # output header
cluster_cmd = 'cluster -i ' + input_file + ' -t ' + str(threshold) + ' --mm'
cluster_cmd_vox = 'cluster -i ' + input_file + ' -t ' + str(threshold)
dummy, f = os.popen2(cluster_cmd_vox) # f will become the output of cluster_cmd
count = 0
for line in f.readlines()[1:]:
	line = line.split()
	count += 1
	if line[1] < min_cluster_size:
		continue
	line[0] = str(count)
	print string.join(line[0:6], '\t')

print "\n\n"


if BRAIN is True:
	import nibabel

	file_loc, name = os.path.split(input_file)
	name = os.path.splitext(os.path.splitext(name)[0])[0]
	thresholded_image_name = name + "_z" + str(threshold) + ".nii.gz"
	cluster_cmd2 = 'cluster -i ' + input_file + ' -t ' + str(threshold) + ' -o ' + os.path.join(file_loc, thresholded_image_name)
	brain_output = os.path.join(file_loc, (name + "_z" + str(threshold) + "_clstrs") )
	
	mask_image_header = nibabel.load(input_file)
	mask_image_data = mask_image_header.get_data()
	mask_image_data = zeros_like(mask_image_data) # set mask data matrix to 0
	# mask_image = NiftiImage(input_file)
	# mask_image.data = zeros_like(mask_image.data) # set mask data matrix to 0
	dummy, f = os.popen2(cluster_cmd2) # f will become the output of cluster_cmd
	sleep(1) # need to wait for the fsl cluster command to finish creating output file
	stat_image_header = mask_image_header
	stat_image_data = stat_image_header.get_data()
	# stat_image = NiftiImage( os.path.join(file_loc, thresholded_image_name) )
	
	print "CLUSTER LOCI IN IMAGE SPACE"
	i = 0
	for line in f.readlines()[1:]:
		
		i += 1
		line = line.split()
		current_cluster_size = line[1]
		
		if int(current_cluster_size) < min_cluster_size: # continue loop if cluster size is too small, will eventually terminate
			continue
		
		current_z_max, x, y, z  = line[2:6]
		x, y, z = int(x), int(y), int(z)
		print "Cluster", (i), "x", x, "y", y, "z", z
		
		if MASKS == 1:
			
			brain_output = os.path.join(file_loc, (str(date.today().year) + str(date.today().month) + str(date.today().day) + "_" + name + "_z" + str(threshold) + "_cls" + str("%.2i" % i) + "_vx" + str(x) + str(y) + str(z)) )
			#mask_image_data[x, y, z] = current_z_max # 1x1x1 vox
			mask_image_data[(x-1):(x+2), (y-1):(y+2), (z-1):(z+2)] = current_z_max # 3x3x3 vox
			output_data = nibabel.Nifti1Image(mask_image_data, stat_image_header.get_affine()) # affine transformation info is required
			sleep(0.5)
			nibabel.save(output_data, brain_output)
			mask_image_data = zeros_like(mask_image_data) # reset mask data matrix to 0
			
		elif MASKS == 2:
			
			brain_output = os.path.join(file_loc, (str(date.today().year) + str(date.today().month) + str(date.today().day) + "_" + name + "_z" + str(threshold) + 
				"_func_mask_cls" + str("%.2i" % i)) )
			mask_image_data = get_functional_roi(stat_image_data, z, y, x)
			mask_image_data = nibabel.Nifti1Image(mask_image_data, stat_image_header.get_affine()) # affine transformation info is required
			sleep(0.5)
			nibabel.save(mask_image_data, brain_output)
			
		else:
			mask_image_data[(x-1):(x+2), (y-1):(y+2), (z-1):(z+2)] = i #current_z_max # 3x3x3 vox

	if MASKS == 0:
		mask_image_data = nibabel.Nifti1Image(mask_image_data, stat_image_header.get_affine()) # affine transformation info is required
		nibabel.save(mask_image_data, brain_output)




#################################################
#				BEGIN TALAIRAC SEARCH
#################################################

#	COLLECTING TALAIRAC LABELS - MUST HAVE ONLINE CONNECTION AND JAR FILE
if TAL is True:
	dummy, f = os.popen2(cluster_cmd) # f will become the output of cluster_cmd
	count = 0
	for line in f.readlines()[1:]:
		count += 1
		line = line.split()
		current_cluster_size = line[1]
		if int(current_cluster_size) < min_cluster_size: # continue loop if cluster size is too small, will eventually terminate
			print "CLUSTER", count, "LESS THAN", min_cluster_size, "VOXELS"
			continue
		current_z_max = line[2]
		x, y, z = line[3:6]
		tal_xyz = _mni2tal_(int(x), int(y), int(z))
		names_list = _tal_names_(talairach_jar, tal_search_area, tal_xyz[0], tal_xyz[1], tal_xyz[2])
		print "CLUSTER NUMBER ", count
		print "size	z_max	MNI_x	MNI_y	MNI_z	TAL_x	TAL_y	TAL_z"
		print "%s	%s	%s	%s	%s	%.1f	%.1f	%.1f" % (current_cluster_size, current_z_max, x, y, z, tal_xyz[0], tal_xyz[1], tal_xyz[2])
		print "%s\t%s" % (names_list[0].split(':')[1].split(',')[0], names_list[0].split(',')[1]) # major regions are usually the same, just print here
		print 'Hits\tLevel 3\tLevel 4\tLevel 5'
		for name in names_list:
			# the output string from tal is separated by a : and commas. Output will be tab-delimited for use with Excel.
			current_name = name.rstrip().split(',')
			if len(current_name) > 0:
				if len(current_name) >3: # some output won't have a Level 5, this makes sure there is one before writing, or else...
					print "%s\t%s\t%s\t%s" % (current_name[0].split(':')[0], current_name[2], current_name[3], current_name[4])
				elif len(current_name) > 2:
					print "%s\t%s\t%s" % (current_name[0].split(':')[0], current_name[2], current_name[3])
		print "\n"
