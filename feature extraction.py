import os
import numpy as np
import radiomics
import pandas as pd
import six
import SimpleITK as sitk

def catch_features(imagePath, maskPath):
    '''
    结合脑膜瘤患者的mri扫描文件和mask文件（VOI）利用pyradiomics函数包提取影像图特征
    输入参数：
    imagePath：mri扫描文件路径
    maskPath：mask文件路径
    返回值：
    feature：提取到的特征值
    name：特征名
    '''
    if imagePath is None or maskPath is None:  # 出现问题，PyRadiomics 也会记录错误日志
        raise Exception('获取测试用例时出错！')  # 抛出异常以防止后续单元格运行（在“运行全部”情况下）
    settings = {}
    settings['binWidth'] = 25  # 直方图箱宽度，离散化
    settings['sigma'] = [3, 5] #指定LoG和小波滤波器的标准差
    settings['Interpolator'] = sitk.sitkBSpline #图像插值器类型
    settings['resampledPixelSpacing'] = [1, 1, 1]  # 重采样像素间距
    settings['normalize'] = True #归一化
    extractor =  radiomics.featureextractor.RadiomicsFeatureExtractor(**settings)
    #print('提取参数：\n\t', extractor.settings)

    extractor.enableImageTypeByName('LoG')
    extractor.enableImageTypeByName('Wavelet')
    extractor.enableAllFeatures()
    extractor.enableFeaturesByName(firstorder=['Energy', 'TotalEnergy', 'Entropy', 'Minimum', '10Percentile', '90Percentile', 'Maximum', 'Mean', 'Median', 'InterquartileRange', 'Range', 'MeanAbsoluteDeviation', 'RobustMeanAbsoluteDeviation', 'RootMeanSquared', 'StandardDeviation', 'Skewness', 'Kurtosis', 'Variance', 'Uniformity'])
    extractor.enableFeaturesByName(shape=['VoxelVolume', 'MeshVolume', 'SurfaceArea', 'SurfaceVolumeRatio', 'Compactness1', 'Compactness2', 'Sphericity', 'SphericalDisproportion','Maximum3DDiameter','Maximum2DDiameterSlice','Maximum2DDiameterColumn','Maximum2DDiameterRow', 'MajorAxisLength', 'MinorAxisLength', 'LeastAxisLength', 'Elongation', 'Flatness'])
    # 将一阶特征和形状特征中的默认禁用的特征都手动启用，为了之后特征筛选
    #print('启用的过滤器：\n\t', extractor.enabledImagetypes)
    feature_cur = []
    feature_name = []
    result = extractor.execute(imagePath, maskPath)

    for key, value in six.iteritems(result):
        #print('\t', key, ':', value)
        feature_name.append(key)
        feature_cur.append(value)
    #print(len(feature_cur[37:]))
    name = feature_name[37:]
    name = np.array(name)
    for i in range(len(feature_cur[37:])):
        #if type(feature_cur[i+22]) != type(feature_cur[30]):
        feature_cur[i+37] = float(feature_cur[i+37])
    return feature_cur[37:],name



def batch_features_catch(file_path):
    '''
    批量提取mri影像特征
    输入参数：
    file_path：mri影像文件夹路径
    返回值：
    final_feature：影像特征矩阵
    '''
    t1_fileNames = os.listdir(file_path+'\\t1')
    t1c_fileNames = os.listdir(file_path+'\\t1c')
    t2_fileNames = os.listdir(file_path+'\\t2')
    t2f_fileNames = os.listdir(file_path+'\\t2f')
    mask_fileNames =  os.listdir(file_path+'\\mask')

    patience = '-'.join(mask_fileNames[0].split('-')[0:3]) #获取文件名中患者信息用于匹配其他文件夹名
    seq = ['t1.nii.gz','t1c.nii.gz','t2.nii.gz','t2f.nii.gz'] 
    files = [t1_fileNames,t1c_fileNames,t2_fileNames,t2f_fileNames]
    mask_Dir = file_path+'\\mask\\'  + mask_fileNames[0]

    file_name = patience + '-' + seq[0] #将患者信息与序列名组合形成完整文件名
    if file_name in files[0]: #判断文件名是否存在，因为有些患者没有四种扫描序列影像
        img_seq = seq[0].split('.')[0].strip() #获取当前序列类型名
        seq_dir = file_path+'\\' +  img_seq + '\\' + file_name 
        print(seq_dir,mask_Dir)
        #进行特征提取
        feature,name=catch_features(seq_dir,mask_Dir)
        feature=pd.DataFrame([feature])
        feature.columns = name
        label = [{'patience_information':patience,
                    'sequence_type': img_seq}]
        label=pd.DataFrame(label)
        final_feature = pd.concat([label,feature],axis=1)
        
    for k in range(1,4):
            file_name = patience + '-' + seq[k] #将患者信息与序列名组合形成完整文件名
            if file_name in files[k]: #判断文件名是否存在，因为有些患者没有四种扫描序列影像
                img_seq = seq[k].split('.')[0].strip() #获取当前序列类型名
                seq_dir = file_path+'\\' +  img_seq + '\\' + file_name 
                print(seq_dir,mask_Dir)
                #进行特征提取
                feature,name=catch_features(seq_dir,mask_Dir)
                feature=pd.DataFrame([feature])
                feature.columns = name
                label = [{'patience_information':patience,
                        'sequence_type': img_seq}]
                label=pd.DataFrame(label)
                patience_feature = pd.concat([label,feature],axis=1)
                final_feature = pd.concat([final_feature,patience_feature])

    for i in mask_fileNames[1:]:
        print(i)
        patience = '-'.join(i.split('-')[0:3]) #获取文件名中患者信息用于匹配其他文件夹名
        seq = ['t1.nii.gz','t1c.nii.gz','t2.nii.gz','t2f.nii.gz'] 
        files = [t1_fileNames,t1c_fileNames,t2_fileNames,t2f_fileNames]
        mask_Dir = file_path+'\\mask\\'  + i 
        for k in range(4):
            file_name = patience + '-' + seq[k] #将患者信息与序列名组合形成完整文件名
            if file_name in files[k]: #判断文件名是否存在，因为有些患者没有四种扫描序列影像
                img_seq = seq[k].split('.')[0].strip() #获取当前序列类型名
                seq_dir = file_path+'\\' +  img_seq + '\\' + file_name 
                print(seq_dir,mask_Dir)
                #进行特征提取
                feature,name=catch_features(seq_dir,mask_Dir)
                feature=pd.DataFrame([feature])
                feature.columns = name
                label = [{'patience_information':patience,
                        'sequence_type': img_seq}]
                label=pd.DataFrame(label)
                patience_feature = pd.concat([label,feature],axis=1)
                final_feature = pd.concat([final_feature,patience_feature])
    return(final_feature)


