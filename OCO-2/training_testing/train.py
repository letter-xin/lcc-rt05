from tokenize import Double
from turtle import shape
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
import tensorflow.keras.models
import tensorflow.keras
from tensorflow.keras import regularizers
import matplotlib
from scipy.stats import gaussian_kde
from tensorflow.python.keras import callbacks
from tensorflow.python.training.tracking.util import Checkpoint
matplotlib.use('Agg')
import tensorflow as tf
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import os
import time
os.environ["CUDA_VISIBLE_DEVICES"]="0"
scaler = StandardScaler()

config = tf.compat.v1.ConfigProto() 
config.gpu_options.per_process_gpu_memory_fraction = 0.95 # 占用GPU50%的显存 
session = tf.compat.v1.Session(config=config)

def r2(y_true, y_pred):
    """
    Calculation of the unadjusted r-squared, goodness of fit metric
    """
    sse  = np.square( y_pred - y_true ).sum()
    sst  = np.square( y_true - y_true.mean() ).sum()
    return 1 - sse/sst

def soundingNum(soundingID):
	id = np.double(soundingID)
	return np.floor((id+0.01)%10)

def read_file(file_name):
	file = open(file_name)
	X = []
	k = 0
	while 1:
		line = file.readline()
		k = k+1
		if k==420000:
			break
		if not line:
			break
		data = line[:-1].split(',')
		X.append(np.array(data).astype(np.float64))
	X = np.array(X)
	return X

t1 = time.time()
X = np.asarray(read_file("trainingdata/lt/filter/spec/specAsia16_18_flag.dat"))
y = np.asarray(read_file("trainingdata/lt/filter/train/trainAsia16_18_flag.dat"))
# X_1 = np.asarray(read_file("trainingdata/lt/filter/spec/specAsia16_18.dat"))
# y_1 = np.asarray(read_file("trainingdata/lt/filter/train/trainAsia16_18.dat"))
# X_2 = np.asarray(read_file("trainingdata/lt/filter/spec/specAsia20.dat"))
# y_2 = np.asarray(read_file("trainingdata/lt/filter/train/trainAsia20.dat"))
# X = np.append(X_1,X_2,axis=0)
# y = np.append(y_1,y_2,axis=0)

#X = X[:,list_num]
Geo = y[:,4:9]
y = y[:,9:11]
#y1 = y1*10000
# y2 = y[:,30]
# y2 = y2/100000
# y = np.delete(y,[0,1,2,3,4,5,6,7,8,9],axis=1)
#y = np.transpose(np.vstack((y1,y1)))
# for i in range(X.shape[0]):
# 	X[i,1] = soundingNum(X[i,0])
X = np.delete(X,[0],axis=1)
X = np.append(X,Geo,axis=1)

t2 = time.time()
print('data loading time ='+str(t2-t1))

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)


scaler.fit(X)  
[ux,sx]=[scaler.mean_,scaler.var_]
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

scaler.fit(y)  
[uy,sy]=[scaler.mean_,scaler.var_]
y_train = scaler.transform(y_train)
y_test = scaler.transform(y_test)

t3 = time.time()
print('data normalization time ='+str(t3-t2))

np.savetxt("../output/lt/flag/ux.dat",ux, fmt="%15.5e",delimiter=',')
np.savetxt("../output/lt/flag/sx.dat",sx, fmt="%15.5e",delimiter=',')
np.savetxt("../output/lt/flag/uy.dat",uy, fmt="%15.5e",delimiter=',')
np.savetxt("../output/lt/flag/sy.dat",sy, fmt="%15.5e",delimiter=',')
layer=[500,200,200,200,250,50]

n_features = X.shape[1]
filepath="../output/lt/flag/train.h5"

model = Sequential()
model.add(Dense(layer[0], activation='relu', kernel_initializer='he_uniform', input_shape=(n_features,)))
model.add(Dense(layer[1], activation='relu', kernel_initializer='he_uniform'))
model.add(Dense(layer[2], activation='relu', kernel_initializer='he_uniform'))
model.add(Dense(layer[3], activation='relu', kernel_initializer='he_uniform'))
# model.add(Dense(layer[4], activation='relu', kernel_initializer='he_uniform'))
# model.add(Dense(layer[5], activation='relu', kernel_initializer='he_uniform'))
model.add(Dense(y.shape[1]))
#model.load_weights(filepath)
#par_model = tensorflow.keras.utils.multi_gpu_model(model,gpus=2)
model.compile(optimizer='adam', loss='mse',metrics=['accuracy'])
checkpoint=callbacks.ModelCheckpoint(filepath,monitor='val_loss',verbose=0,save_best_only=True,mode='min')
callbacks_list=[checkpoint]
aa=model.fit(X_train, y_train,epochs=500, batch_size=64,validation_data=(X_test, y_test),callbacks=callbacks_list, verbose=1)
#model.fit(X_train, y_train,epochs=1000, batch_size=128,validation_split=0.1,callbacks=callbacks_list, verbose=1)
# loss = aa.history['loss']
# epochs = range(1, len(loss) + 1)
# plt.title('Loss')
# plt.ylim(0,2e-5)
# plt.xlim(1500,5000)
# plt.plot(epochs, loss, 'blue', label='loss')
# plt.legend()
# plt.savefig('output/loss.png')
res = model.predict(X_test)
res = res * (sy**(0.5)) + uy
y_test = y_test * (sy**(0.5)) + uy




LBL = y_test[:,0]
MLP = res[:,0]
xy = np.vstack([MLP, LBL])
z = gaussian_kde(xy)(xy)
z=(z-np.min(z))/(np.max(z)-np.min(z))
idx = z.argsort()
MMLP, LLBL, z = MLP[idx], LBL[idx], z[idx]
xyline = np.linspace(0.95*np.min(LBL),1.05*np.max(LBL),61)
fig = plt.figure(figsize=(8,7))
figg,ax = plt.subplots()
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.ylabel('Predicted[ppm]',fontname="serif",fontsize=14)
plt.xlabel('L2[ppm]',fontname="serif",fontsize=14)
plt.axis([395,420,395,420])
text1 = 'N:'+str(len(LBL)) 
text2 = 'R$^2$:'+str("%.3f" % r2(LBL,MLP))
text3 = 'ME: '+str("%.3f" % np.mean(LBL-MLP))+' ppm'
text4 = 'RMSE: '+str("%.3f" % np.sqrt(mean_squared_error(LBL, MLP)))+' ppm'
plt.plot(xyline, xyline, 'r-',label='Ideal$\pm$1%', linewidth=2)
plt.fill_between(xyline, xyline*1.01, xyline*0.99,
    alpha=0.5, edgecolor='0.4', facecolor='0.4',
    linewidth=1, linestyle='--', antialiased=True)
plt.scatter(LLBL, MMLP, c=z,s=7, alpha=1.0)#facecolors='none',
plt.colorbar(label='Number density')
plt.tick_params(labelsize=14)
plt.legend(fontsize=14)
plt.ticklabel_format(style='plain', scilimits=(0, 0), axis='both')
plt.text(0.02,0.85,text1,ha='left',va='top',transform=ax.transAxes,fontname="serif",fontsize=14)
plt.text(0.02,0.80,text2,ha='left',va='top',transform=ax.transAxes,fontname="serif",fontsize=14)
plt.text(0.02,0.75,text3,ha='left',va='top',transform=ax.transAxes,fontname="serif",fontsize=14)
plt.text(0.02,0.70,text4,ha='left',va='top',transform=ax.transAxes,fontname="serif",fontsize=14)
plt.tight_layout()
plt.savefig("../plot/lt/flag/scatterx.png",dpi=300)


LBL = y_test[:,1]/1000
MLP = res[:,1]/1000
xy = np.vstack([MLP, LBL])
z = gaussian_kde(xy)(xy)
z=(z-np.min(z))/(np.max(z)-np.min(z))
idx = z.argsort()
MMLP, LLBL, z = MLP[idx], LBL[idx], z[idx]
xyline = np.linspace(0.95*np.min(LBL),1.05*np.max(LBL),61)
fig = plt.figure(figsize=(8,7))
figg,ax = plt.subplots()
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.ylabel('Predicted[hPa]',fontname="serif",fontsize=14)
plt.xlabel('L2[hPa]',fontname="serif",fontsize=14)
plt.axis([0.99*np.min(LBL), 1.01*np.max(LBL), 0.99*np.min(LBL), 1.01*np.max(LBL)])
text1 = 'N:'+str(len(LBL)) 
text2 = 'R$^2$:'+str("%.3f" % r2(LBL,MLP))
text3 = 'ME: '+str("%.3f" % np.mean(LBL-MLP))+' hPa'
text4 = 'RMSE: '+str("%.3f" % np.sqrt(mean_squared_error(LBL, MLP)))+' hPa'
plt.plot(xyline, xyline, 'r-',label='Ideal$\pm$1%', linewidth=2)
plt.fill_between(xyline, xyline*1.01, xyline*0.99,
    alpha=0.5, edgecolor='0.4', facecolor='0.4',
    linewidth=1, linestyle='--', antialiased=True)
plt.scatter(LLBL, MMLP, c=z,s=7, alpha=1.0)#facecolors='none',
plt.colorbar(label='Number density')
plt.tick_params(labelsize=14)
plt.legend(fontsize=14)
plt.ticklabel_format(style='plain', scilimits=(0, 0), axis='both')
plt.text(0.02,0.85,text1,ha='left',va='top',transform=ax.transAxes,fontname="serif",fontsize=14)
plt.text(0.02,0.80,text2,ha='left',va='top',transform=ax.transAxes,fontname="serif",fontsize=14)
plt.text(0.02,0.75,text3,ha='left',va='top',transform=ax.transAxes,fontname="serif",fontsize=14)
plt.text(0.02,0.70,text4,ha='left',va='top',transform=ax.transAxes,fontname="serif",fontsize=14)
plt.tight_layout()
plt.savefig("../plot/lt/flag/scatterp.png",dpi=300)

