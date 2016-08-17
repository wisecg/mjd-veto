f = open('P3KJR_MasterRunList.txt', 'w')
with open('p3kjr_runRanges.txt') as file:

	for line in file:

		col = line.split()
		r1 = int(col[0])
		r2 = int(col[1])
		for i in range(r1,r2):
			f.write(str(i)+"\n")