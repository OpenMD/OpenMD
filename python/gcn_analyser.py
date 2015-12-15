__author__ = 'jmichalk'


import sys
import os.path
import re
import math

#region Vector manipulation and distance manipulation




# A class to hold information about each frame
class Frame:
    """Frame class"""
    def __init__(self, regions):
        self.regions = regions


#endregion


#region file input
#region accepting file
if len(sys.argv) == 1:
    filename = ""
else:
    try_filename = sys.argv[1]
    if "gcn" in try_filename:
        filename = try_filename
        print "Filename: %s" % filename
    else:
        print "(" + try_filename + ") contained errors, please submit a valid gcn file."
        filename = ""

while(filename == ""):
    print "Please enter a filename to process: "
    try_filename = raw_input("Filename: ---> ")
    if "gcn" in try_filename:
        filename = try_filename
    else:
        print "(" + try_filename + ") contained errors, please submit a valid gcn file."

file_exists = os.path.isfile(filename)
if (file_exists == False):
    print "File %s does not exist" % (filename)
    print "Quitting"
    exit(0)

#endregion


my_frames = []
my_regions = []
my_bins = []

#region reading file
file_handle = open(filename, 'r')

for line in file_handle:
    if "Frame" in line:
        currentFrame = int(re.findall(r'\d+', line)[0])
        if currentFrame > 0:
            newFrame = Frame(my_regions)
            my_frames.append(newFrame)
            my_regions = []


        #print "Current Frame in read in: %d" % currentFrame

    if "Region" in line:
        regionList = line.split()
        my_regions = [[] for i in range(0,len(regionList)-1)]
        #print "Number of Regions: %d" % (len(regionList)-1)

    if "#" not in line:
        data = line.split()
        if currentFrame == 0:
            my_bins.append(float(data[0]))

        for i in range(0,len(regionList)-1):
            my_regions[i].append(float(data[i+1]))

newFrame = Frame(my_regions)
my_frames.append(newFrame)
#endregion


#endregion

#region Analyze/Shrink data
#TODO possibility of expanding past just collapsing number of frames (but this makes the most sense at first)
#TODO perhaps collapse number of regions
collapsedFrames = []
collapsedRegion = [[0 for j in range(0, len(my_frames[0].regions[0]))] for k in range(0, len(regionList)-1)]
ratio = float(raw_input("Please enter the ratio of frames you wish to have (10->2=5): "))
newFrameNumber = len(my_frames)/ratio

print "Current frames: %d\t New Frames: %d\t Approximate Frame ratio: %f\n" % (len(my_frames), newFrameNumber, ratio)

count = 0
for i in range(0, len(my_frames)):

    if count >= int(ratio):
        count = 0
        for r in collapsedRegion:
            for index in range(0,len(r)):
                r[index] = float(r[index]/ratio)
        newFrame = Frame(collapsedRegion)
        collapsedFrames.append(newFrame)
        #collapsedRegion = [[] for i in range(0, len(regionList)-1)]
        collapsedRegion = [[0 for j in range(0, len(my_frames[0].regions[0]))] for k in range(0, len(regionList)-1)]

    currentFrame = my_frames[i]
    for regionIndex in range(0, len(currentFrame.regions)):
        for dataIndex in range(0, len(currentFrame.regions[regionIndex])):
            collapsedRegion[regionIndex][dataIndex] += currentFrame.regions[regionIndex][dataIndex]

    count = count + 1

for r in collapsedRegion:
    for index in range(0, len(r)):
        r[index] = float(r[index]/count)

newFrame = Frame(collapsedRegion)
collapsedFrames.append(newFrame)
print "Collapsed Frames: %d" % len(collapsedFrames)
#endregion


print collapsedFrames
print collapsedFrames[0].regions
print collapsedFrames[0].regions[0]



#region Output
#TODO Possibility of multiple output files
outfilename=raw_input("Please enter name of file(s) to create: ")
for findex in range(0, len(collapsedFrames)):
    tempfilename = outfilename + ".gcnf" + str(findex)
    out=open(tempfilename, "w")
    out.write("#Made from file " + filename + "\n")
    out.write("#Frame " + str(findex) + "\n")
    out.write("#Original Frames: " + str(len(my_frames)) + "\tCollapsed Frame: " + str(len(collapsedFrames)) + "\n")
    frame = collapsedFrames[findex]
    region = frame.regions[0]
    for index in range(0,len(region)):
        out.write(str(my_bins[index]) + "\t")
        for r in frame.regions:
            out.write(str(r[index]) + "\t")
        out.write("\n")

    out.close()


#endregion
