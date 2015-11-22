'''
#Generate test data

import numpy as np

vps=["vp01","vp02","vp03","vp04"]

for vp in vps:
	noise_array=np.zeros((512,61))

	for row in range(512):
		if row%4 == 0:
			print '0'
			noise_array[row,0:30]=np.random.normal(0,1,30)
			noise_array[row,30:61]=np.random.normal(10,1,31)
		elif row%4 == 1:
			print '1'
			noise_array[row,0:15]=np.random.normal(5,1,15)
			noise_array[row,15:40]=np.random.normal(1,1,25)
			noise_array[row,40:61]=np.random.normal(10,1,21)
		elif row%4 == 2:
			print '2'
			noise_array[row,0:15]=np.random.normal(10,1,15)
			noise_array[row,15:40]=np.random.normal(10,1,25)
			noise_array[row,40:61]=np.random.normal(0,1,21)
		elif row%4 == 3:
			print '3'
			noise_array[row,0:5]=np.random.normal(0,1,5)
			noise_array[row,5:50]=np.random.normal(10,1,45)
			noise_array[row,50:61]=np.random.normal(20,1,11)

			
	np.savetxt('{0}_O1.txt' .format(vp), noise_array)
'''


#Generate pm mean models based test data

import numpy as np

groups=["HEA", "PAT"]
vps=["vp01","vp02","vp03","vp04"]
conds=["Rest", "Test"]
runs=list(range(3))

four_maps=np.loadtxt('E://Programming//DEVELOPMENT//keypy_public//keypy//tests//6_sortmaps//data//mean_models_milz_etal_2015.asc')


for group in groups:
    for vp in vps:
        for cond in conds:
            for run in runs:
                noise_array=np.zeros((512,61))

                for row in range(512):
                    if row%4 == 0:
                        print('0')
                        noise_array[row,:]=four_maps[0]+np.random.normal(0,0.01,61)
                    elif row%4 == 1:
                        print('1')
                        noise_array[row,:]=four_maps[1]+np.random.normal(0,0.01,61)
                    elif row%4 == 2:
                        print('2')
                        noise_array[row,:]=four_maps[2]+np.random.normal(0,0.01,61)
                    elif row%4 == 3:
                        print('3')
                        noise_array[row,:]=four_maps[3]+np.random.normal(0,0.01,61)

                np.savetxt('{0}_{1}_{2}_{3}.txt' .format(group, vp, cond, run), noise_array)