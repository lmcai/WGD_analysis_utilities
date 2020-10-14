import sys

x=open(sys.argv[1]).readlines()
y=open(sys.argv[1]+'.top','a')
cur_rec=''
i=0
for l in x:
	if l.split()[0]!=cur_rec:
		cur_rec=l.split()[0]
		if l.split()[0]!=l.split()[1]:
			i=1
			y.write(l)
	else:
		if i<int(sys.argv[2]):
			if l.split()[0]!=l.split()[1]:
				y.write(l)
				i=i+1


y.close()