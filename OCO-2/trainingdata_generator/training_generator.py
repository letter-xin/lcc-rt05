import numpy as np
import h5py
import os
import pandas as pd
import time
import re

def FindM(table,target,last_row):
    row =[]
    col =[]
    for i in range(last_row,table.shape[0]):
        for j in range(table.shape[1]):
            if target == table[i][j]:
                row = i
                col = j
                break
        if (row == []) == False:
            break
    return row, col

def FindV(table,target,last_row):
    row =[]
    #print(last_row)
    #print(len(table))
    new_table = table[last_row:,]
    # for i in range(last_row,len(table)):
    #     # print(i)
    #     if target == table[i]:
    #         row = i
    #         break
    #     if (row == []) == False:
    #         break
    # k = last_row
    # while 1:
    #     k = k+1
    #     if k > len(table):
    #         break
    #     line = table[k]
    #     if line == target:
    #         row = k
    #         break
    #     if int(table[k]) > int(target)+0.5:
    #         break
    ss = pd.DataFrame(np.arange(len(new_table)),index=new_table)
    ss.columns=['ss']
    try:
        row = last_row+ss.loc[target].ss
    except:
        row = []
    return row

def geoGET(l2,num:np.int32):
    #get some geometry informations and specialize the list
    AU = 1.49597871 * (10 ** 11)
    time93 = l2['RetrievalHeader']['retrieval_time_tai93'][num]
    latitude = l2['RetrievalGeometry']['retrieval_latitude'][num]
    longitude = l2['RetrievalGeometry']['retrieval_longitude'][num]
    distance = l2['RetrievalGeometry']['retrieval_solar_distance'][num]/AU
    zenith = l2['RetrievalGeometry']['retrieval_zenith'][num]/180*np.pi
    azimuth = l2['RetrievalGeometry']['retrieval_azimuth'][num]/180*np.pi
    solar_zenith = l2['RetrievalGeometry']['retrieval_solar_zenith'][num]/180*np.pi
    solar_azimuth = l2['RetrievalGeometry']['retrieval_solar_azimuth'][num]/180*np.pi
    return [time93,latitude,longitude,distance,zenith,azimuth,solar_zenith,solar_azimuth]

def stateVectorGET(l2,num:np.int32):
    #get some state-vector informations and specialize the list
    surface_pressure = l2['RetrievalResults']['surface_pressure_fph'][num]
    #print(len(co2_profile))
    return [surface_pressure]

def findSoundingID(l1,sounding_id,last_raw):
    sd_id_list = l1['SoundingGeometry']['sounding_id']
    row, col = FindM(sd_id_list,sounding_id,last_raw)
    return row, col

def specGET(l1,sound_id,last_raw):
    #get some state-vector informations and specialize the list
    row, col = findSoundingID(l1,sound_id,last_raw)
    #continuum_o2 = l1['SoundingMeasurements']['rad_continuum_o2'][row][col]
    #continuum_wco2 = l1['SoundingMeasurements']['rad_continuum_weak_co2'][row][col]
    #continuum_sco2 = l1['SoundingMeasurements']['rad_continuum_strong_co2'][row][col]
    spec_o2 = l1['SoundingMeasurements']['radiance_o2'][row][col][99:899]
    spec_o2 = spec_o2/spec_o2[790]
    spec_wco2 = l1['SoundingMeasurements']['radiance_weak_co2'][row][col][99:899]
    spec_wco2 = spec_wco2/spec_wco2[350]
    spec_sco2 = l1['SoundingMeasurements']['radiance_strong_co2'][row][col][99:899]
    spec_sco2 = spec_sco2/spec_sco2[300]
    return [sound_id]+list(spec_o2)+list(spec_wco2)+list(spec_sco2) , row

def trainData(readfile_l2,readfile_l1,readfile_l2_lt,savefile_sv,savefile_spec,land_part):

    #generate training dataset and save as a h5.file 
    #define land_part as read boundary

    savepath= open(savefile_sv,"a+")
    savepath_spec= open(savefile_spec,"a+")
    l1 = h5py.File(readfile_l1,mode='r')
    l2 = h5py.File(readfile_l2,mode='r')
    l2_lt = h5py.File(readfile_l2_lt,mode='r')

    #read sounding-id and xco2 and land_fraction

    sd_id_list = l2['RetrievalHeader']['sounding_id']
    xco2_list = l2['RetrievalResults']['xco2']
    xco2_qual_flag = l2_lt['xco2_quality_flag']
    xco2_lt_list = l2_lt['xco2']
    sd_id_lt_list = l2_lt['sounding_id']
    latitude_list = l2['RetrievalGeometry']['retrieval_latitude']
    longitude_list = l2['RetrievalGeometry']['retrieval_longitude']
    land_frac_list = l2['RetrievalGeometry']['retrieval_land_fraction']
    cloud_list = l2['PreprocessingResults']['cloud_flag_idp']
    k = 0
    j_past = 0
    row_past = 0
    for i in range(sd_id_list.shape[0]):
        if (xco2_list[i] < 0) or (land_frac_list[i] <= land_part):
            continue
        if (latitude_list[i] < 20) or (latitude_list[i] > 45):
            continue
        if (longitude_list[i] < 110) or (longitude_list[i] > 150):
            continue
        if cloud_list[i] <1.8:
            continue
        k = k+1
        print(str(k)+' and '+str(sd_id_list[i]))
        sd_id = sd_id_list[i]
        geo_info = geoGET(l2,num = i)
        j_new = FindV(sd_id_lt_list,sd_id_list[i],j_past)
        spec , row_new = specGET(l1,sound_id = sd_id,last_raw=row_past)
        if (j_new == []) or (row_new ==[]):
            print(j_new)
            print(row_new)
            continue
        if xco2_qual_flag[j_new]>0.1:
            continue
        j_past = j_new
        row_past = row_new
        xco2 = xco2_lt_list[j_new]
        sv = stateVectorGET(l2,num = i)
        #combine all need
        output = [sd_id]+geo_info+[xco2]+sv
        output = ','.join(str(j) for j in output)
        savepath.write(str(output)+'\n')
        spec = ','.join(str(j) for j in spec)
        savepath_spec.write(str(spec)+'\n')
    l1.close()
    l2.close()
    savepath.close()
    savepath_spec.close()


def Generator(savefile_label,savefile_spec,land_part):
    #define land_part as read boundary
    l2_files = os.listdir('../../L2_ND/')
    l2_lt_files = os.listdir('../../L2_Lt/')
    l1_files = os.listdir('../../L1B_ND/')
    k = 0
    for l2_file in l2_files:
        l2_name_path = os.path.join('',l2_file)
        model_name = 'ND'
        if model_name in l2_name_path:
            keywork_location = re.search(model_name,l2_name_path)
            year_num=keywork_location.span()[1]+8
            year = l2_name_path[year_num:year_num+2]
            datemark = l2_name_path[year_num-1:year_num+7]
            #print(datemark)
            if year in '16_18_17':
                sries = l2_name_path[10:26]
                l2_name_path = '../../L2_ND/'+l2_name_path
                for l1_file in l1_files:
                    l1_name_path = os.path.join('',l1_file) 
                    if sries in l1_name_path:
                        #print('fileNo='+sries)
                        l1_name_path = '../../L1B_ND/'+l1_name_path
                        for l2_lt_file in l2_lt_files:
                            l2_lt_name_path = os.path.join('',l2_lt_file)
                            if datemark in l2_lt_name_path:
                                l2_lt_name_path = '../../L2_Lt/'+l2_lt_name_path
                                k = k + 1
                                print('No.'+str(k)+' : '+l1_name_path+' and '+l2_name_path+' and '+l2_lt_name_path)
                                trainData(l2_name_path,l1_name_path,l2_lt_name_path,savefile_label,savefile_spec,land_part)
                                break
                        break
    File_num = k
    return File_num

if __name__ == '__main__':
    t1 = time.time()
    num = Generator('filter/train/trainAsia16_18_flag.dat','filter/spec/specAsia16_18_flag.dat',land_part=99.99)
    t2 = time.time()
    print('Cost time = '+str(t2-t1)+'s, File_Number = '+str(num))