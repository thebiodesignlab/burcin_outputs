"""This script was written by Burcin Acar and is designed to process CSV files containing data
from a microplate reader. The script reads in the CSV file specified by the user, removes any
rows after the "Analysis Result" row, and formats the remaining data into a new CSV file with
a specific filename. The resulting CSV files contain data in a format suitable for downstream
analysis and visualization. The script uses the Python csv and numpy libraries to read and
manipulate the data, and requires a properly formatted input file to run.If you have any questions
or comments about this script, please contact Burcin Acar."""

import csv
import sys
import numpy as np
from datetime import datetime,timedelta

# Check if the filename is provided as an argument
if len(sys.argv) < 2:
    print('Usage: python script.py filename')
    sys.exit()

# Get the filename from command line arguments
filename = sys.argv[1]
print("Reading input data...\n")

# Define the time interval between measurements
interval=10
row_str=""
# Read the CSV file
scount=0
fn = [None] * 2
with open(filename, 'rbU') as csvfile:
    # Create a CSV reader object
	reader = csv.reader(csvfile, delimiter=',')

    # Initialize an empty list to hold the text data
	text_array = []
	barcode_count = 0
	results_count = 0
    # Loop through each row in the CSV file
	for row in reader:
		row_str = ",".join(row)
		barcode_count += row_str.count("Barcode:")
		results_count += row_str.count("Results")
		# If the row contains the string "Analysis Result", stop reading the file
		if row and row[0] and "Analysis Result" in row[0]:
			break
		# If the row contains neither "Barcode" nor "Results", append it to the text array
		if row and row[0] and "Results" not in row[0] and "Barcode" not in row[0] and row[0].isalpha():
			text_array.append(row)
		# If the row contains "Barcode" and no data has been read yet, get the timestamp
		if row and row[0] and "Barcode" in row[0] and len(text_array)==0:
			bs=row[0].split()
			bs1=bs[0].split("V-")
			ts=bs1[1]+' '+bs[1]
		if row and row[0] and "Results" in row[0]:
			sp=row[0].split("\xc2\xa0")
			if "F" in sp[1]:
				fn[scount]="F"
				scount=scount+1
			if "LUM" in sp[1]:
				fn[scount]="L"
			
# Convert the text array to a numpy array
fancy=np.array(text_array)
# Calculate some dimensions and reshape the array
rep=barcode_count/results_count #83
print("Number of measurements: %s" % rep)
rows, cols = fancy.shape # 83*3*8,13
shp=rows/barcode_count #8
# Define a lambda function to handle the range function edge case
f=lambda n: 0 if n == 0 else 1
# Loop through each set of measures(fluorescence, luminescence and OD)
for k in range(results_count):
	barc =np.empty((1,rep),dtype='U100')
	freshaped1=np.empty((shp*(cols-1),rep),dtype='U100')
	well=np.empty((shp*(cols-1),1),dtype='U100')
	pts=datetime.strptime(ts,"%Y-%m-%d %H:%M:%S")
	for i in range(rep*k,rep*(k+1)):
		if i!=rep*k:
			pts = pts + timedelta(minutes=interval)
		fts = datetime.strftime(pts, "%-m/%-d/%y %-H:%M")
		# Add the timestamp to the barcode array
		barc[0,i-k*rep]=fts
		# Loop through each well and add the data to the appropriate arrays
		for j in range(shp*(cols-1)):
			rem = (j) % shp
			res = ((j) // shp) +1	
			freshaped1[j,i-k*rep]=fancy[i*shp+rem,res]
			well[j,0]=fancy[rem,0]+str(res).zfill(2)
	barc=np.insert(barc,0,"Time",axis=1)
# fill the array with empty strings
	freshaped1 = np.hstack([well,freshaped1])
	freshaped1 = np.vstack([barc,freshaped1])
# Print the text array
	if k>0:
		ofilename ="20200530_erika_taga_00001_reader_plate_1_kinetic_supp_3_high_{}.csv".format(fn[k-1])
	else:
		ofilename="20200530_erika_taga_00001_reader_plate_1_kinetic_supp_abs.csv" #OD
	np.savetxt(ofilename, freshaped1, delimiter=',', fmt='%s') #luminiscence IT CAN BE FLUORESCENCE
	del freshaped1
	del barc
	del well
print("The formatted csv files are saved.")
