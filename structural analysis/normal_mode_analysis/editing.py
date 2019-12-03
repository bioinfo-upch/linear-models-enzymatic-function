x = open("output.ready.txt","r")
outfile = open("eigenvector.txt","w")

for i in x:
	i = str(i)
	i = i.replace('\t',' ')
	i = i.replace(' 0.','  0.')
	i = i.replace('CA  ','CA ')
	outfile.write(i)
outfile.close()

