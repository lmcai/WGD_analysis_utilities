#get smoothing value from cv output, and replace the text with divergence time estimation
import fnmatch
import os

#get output filename
filenames=list()
for file in os.listdir('.'):
	 if fnmatch.fnmatch(file, '*.out'):filenames.append(file)
	 
#extract optimum smoothing value
for fn in filenames:
	with open(fn, 'r') as f:
		smoothing='na'
		for line in f:
			if 'Optimum' in line:
				smoothing=line.split()[2]
		f.close()
	
	if smoothing !='na':
		newline='set penalty=add smoothing='+smoothing+';\nDIVTIME method=pl algorithm=tn;\nshowage shownamed=yes;\n'	
		f1=open(fn.split('.')[0]+'.nex','r')
		f2=open(fn.split('.')[0]+'.2.nex','w')
		filedata = f1.readlines()
		for line in filedata:
			f2.write(newline if line.startswith('DIVTIME method=pl') else line)
			#replace old cv command with divtime command
		f1.close()
		f2.close()
	else:
		print fn
