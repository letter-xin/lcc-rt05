import numpy as np
#import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# read data from file
def read_file(file_name):
	file = open(file_name)
	X = []
	while 1:
		line = file.readline()
		if not line:
			break
		data = line[:-1].split(',')
		X.append(np.array(data).astype(np.float))
	return X



X1 = read_file("testdata/features")
#X2 = read_file("testdata/features_noise")


X1 = np.transpose(X1)
#X2 = np.transpose(X2)



#spatial grid
gas_l = np.arange(2300,2400,0.001*1)
#0.48214527964592
    
for i in range(10):
		plt.clf()
		plt.figure(1)
		plt.ylabel('Transmissivity [-] ',fontname="serif")
		plt.xlabel('Wavenumber [cm$^{-1}$] ',fontname="serif")
#		plt.axis([0, 32, 0.0001, 100])
#		plt.plot(gas_l, np.multiply(np.exp(res[:,i]),w), 'bo',gas_l, np.multiply(np.exp(testy[:,i]),w), 'k-',markersize=5)
#		plt.plot(gas_l, X1[:,i], 'b-', gas_l, X2[:,i], 'k-',linewidth=0.5, markersize=1)
		plt.plot(gas_l, X1[:,i], 'b-', linewidth=0.5, markersize=1)
#		plt.yscale('log')
		plt.tight_layout()
		name=str(i)
		plt.savefig('testplot/'+name +'.pdf')

#plt.show()
