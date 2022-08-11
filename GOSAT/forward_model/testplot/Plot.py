import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
matplotlib.use('Agg')
from numpy import zeros
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
############################################################
from sklearn.model_selection import train_test_split
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import LeavePOut
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import StratifiedKFold
############################################################

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

def read_file1(file_name):
	file = open(file_name)
	X = []
	while 1:
		line = file.readline()
		if not line:
			break
		data = line[:-1].split('\t')
		X.append(np.array(data).astype(np.float))
	return X

#####################################
Y_temp_sat = np.asarray(read_file1("SAT_0907260335_1   9  17_(-25.7361 145.7644)_46.0922_29.0366.dat"))
Y_fit = np.asarray(read_file("labels_FIT_090726.dat"))

Y_sat = []
for i in range(Y_temp_sat.shape[1]):
	if 6180 <= Y_temp_sat[0,i] <= 6280.1:
		Y_sat.append(np.array(Y_temp_sat[1:,i]).astype(np.float))

Y_sat = np.asarray(Y_sat).T

wavenumber = np.linspace(6180, 6280, 502)
for i in range(Y_fit.shape[0]):
	plt.figure()
	plt.subplot(211)
	plt.plot(wavenumber,Y_sat[i,:],'k-',linewidth = 0.5, label="GOSAT", markersize=3)
	plt.plot(wavenumber,Y_fit[i,],'r-',linewidth = 0.5, label="ForwardModel", markersize=3)
	plt.legend()
	plt.subplot(212)
	plt.ylim(max(Y_fit[i,:])*0.1,-max(Y_fit[i,:])*0.1)
	plt.plot(wavenumber,Y_fit[i,:]-Y_sat[i,:], 'b', linewidth = 0.5, label = "difference",markersize=3)
	plt.legend()
	name = str(i)
	plt.savefig('Spectrum_090726_' + name + '.png',dpi=300)



