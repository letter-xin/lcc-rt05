import h5py
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
import tensorflow.keras.models
from scipy.stats import gaussian_kde
import tensorflow.keras
import tensorflow as tf
import matplotlib
from tensorflow.python.keras import callbacks
from tensorflow.python.training.tracking.util import Checkpoint
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.basemap import Basemap
import os
import time
#os.environ["CUDA_VISIBLE_DEVICES"]="0,1"
scaler = StandardScaler()

# config = tf.compat.v1.ConfigProto() 
# config.gpu_options.per_process_gpu_memory_fraction = 0.9 # 占用GPU50%的显存 
# session = tf.compat.v1.Session(config=config)

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

def ReLU(x):
    return np.maximum(x,0)

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
X = np.asarray(read_file("trainingdata/lt/filter/spec/specAsia19_flag.dat"))
y = np.asarray(read_file("trainingdata/lt/filter/train/trainAsia19_flag.dat"))
# date = y[:,0]
# lati = y[:,2]
# long = y[:,3]


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

#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.05, random_state=1)
ux = np.loadtxt("../output/lt/flag/ux.dat")
sx = np.loadtxt("../output/lt/flag/sx.dat")
uy = np.loadtxt("../output/lt/flag/uy.dat")
sy = np.loadtxt("../output/lt/flag/sy.dat")

X = (X - ux)/(sx**0.5)
X[:,790] = np.zeros(X.shape[0])
X[:,1150] = np.zeros(X.shape[0])
X[:,1900] = np.zeros(X.shape[0])
y = (y - uy)/(sy**0.5)

#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)
#y = y_test
t3 = time.time()
print('data normalization time ='+str(t3-t2))


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
#par_model = tensorflow.keras.utils.multi_gpu_model(model,gpus=2)
model.compile(optimizer='adam', loss='mse',metrics=['accuracy'])
checkpoint=callbacks.ModelCheckpoint(filepath,monitor='val_loss',verbose=1,save_best_only=True,mode='min')
callbacks_list=[checkpoint]
model.load_weights(filepath)
# model.compile(optimizer='adam', loss='mse',metrics=['accuracy'])
# checkpoint=callbacks.ModelCheckpoint(filepath,monitor='loss',verbose=1,save_best_only=True,mode='min')
# callbacks_list=[checkpoint]
# aa=model.fit(X, y,epochs=0, batch_size=256,callbacks=callbacks_list, verbose=1)
res = model.predict(X)

res = res * (sy**(0.5)) + uy
y = y * (sy**(0.5)) + uy

t4 = time.time()
print('data retrieval time ='+str(t4-t3))



LBL = y[:,0]
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
plt.axis([390,425,390,425])
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
plt.savefig("../plot/lt/flag/scatterx19.png",dpi=300)


LBL = y[:,1]/1000
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
plt.axis([0.97*np.min(LBL), 1.03*np.max(LBL), 0.97*np.min(LBL), 1.03*np.max(LBL)])
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
plt.savefig("../plot/lt/flag/scatterp19.png",dpi=300)










# d = res - y





# ##----------------map------------
# for month in range(4):
# 	begindate = 2019000100000000+(month*3+3)*10000000000
# 	enddate = 2019003100000000+(month*3+3)*10000000000
	
# 	detlaxco2 = []
# 	deltap = []
# 	latitude = []
# 	longitude = []
# 	for i in range(y.shape[0]):
# 		if (date[i]>begindate) and (date[i]<enddate):
# 			latitude=np.hstack((latitude,lati[i]))
# 			longitude=np.hstack((longitude,long[i]))
# 			detlaxco2 = np.hstack((detlaxco2,d[i,0]))
# 			deltap = np.hstack((deltap,d[i,1]))
# 	print(str(begindate)+'   '+str(enddate))
# 	print(latitude)
# 	print(longitude)
# 	print(detlaxco2)
# 	print('----------------------------')
	
# 	units='XCO$_{2}$[ppm]'
# 	fig = plt.figure(figsize=(8, 6), dpi=1000)
# 	m = Basemap(projection='cyl', resolution='h', llcrnrlat=20, urcrnrlat=45,llcrnrlon=110, urcrnrlon=150)
# 	#shp_info = m.readshapefile("CHN_adm_shp/CHN_adm1",'states',drawbounds=False)
# 	m.fillcontinents(color='white', lake_color='lightskyblue')
# 	m.drawmapboundary(fill_color='skyblue')
# 	m.drawmeridians(np.arange(110, 150, 10), labels=[1,0,0,1])
# 	m.drawparallels(np.arange(20, 45.001, 5), labels=[1,0,0,1])
# 	m.scatter(longitude, latitude,c=detlaxco2*1000000, s=0.2,zorder=10,cmap=plt.cm.jet,)
# 	cb = m.colorbar(location="bottom",pad='7%')
# 	cb.set_clim(-20,20)
# 	cb.set_label(units)
# 	plt.gcf()
# 	plt.title('2019 MLP-OCO Asia area')
# 	#plt.tight_layout()
# 	pngfile = "../plot/Time/Time19_asia_"+str(month)+"_x.png"
# 	fig.savefig(pngfile)

# 	units='Pressure[bar]'
# 	fig = plt.figure(figsize=(8, 6), dpi=1000)
# 	ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
# 	m = Basemap(projection='cyl', resolution='h', llcrnrlat=20, urcrnrlat=45,llcrnrlon=110, urcrnrlon=150)
# 	#shp_info = m.readshapefile("CHN_adm_shp/CHN_adm1",'states',drawbounds=False)
# 	m.fillcontinents(color='white', lake_color='lightskyblue')
# 	m.drawmapboundary(fill_color='skyblue')
# 	m.drawmeridians(np.arange(110, 150, 10), labels=[1,0,0,1])
# 	m.drawparallels(np.arange(20, 45.001, 5), labels=[1,0,0,1])
# 	m.scatter(longitude, latitude,c=deltap/100000, s=0.2,zorder=10,cmap=plt.cm.jet)#,cmap=plt.cm.jet,edgecolors=None, linewidth=0)
# 	cb = m.colorbar(location="bottom",pad='7%')
# 	cb.set_clim(-0.05,0.05)
# 	cb.set_label(units)
# 	plt.gcf()
# 	plt.title('2019 MLP-OCO Asia area')
# 	#plt.tight_layout()
# 	pngfile = "../plot/Time/Time19_asia_"+str(month)+"_p.png"
# 	fig.savefig(pngfile)









