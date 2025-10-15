#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 22:09:09 2025

@author: lacter
"""

import pandas as pd
import time
from progress.bar import Bar
import pyopenms
import math
import numpy as np
import bisect
import re
import os
from scipy.optimize import linear_sum_assignment
from scipy.spatial import KDTree
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder, StandardScaler
from scipy.stats import chi2
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from statsmodels.nonparametric.smoothers_lowess import lowess
#from sklearn.linear_model import ElasticNetCV
from sklearn.cross_decomposition import PLSRegression
#from scipy.cluster.hierarchy import linkage, leaves_list
#import plotly.figure_factory as ff
#from plotly.subplots import make_subplots
#from sklearn.model_selection import train_test_split
#from sklearn.metrics import accuracy_score
#from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import warnings
import glob
#import json
import sqlite3
'''---------'''
import sys
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QGridLayout, QLabel, QTableWidget, QHeaderView,
    QPushButton, QLineEdit, QProgressBar, QDesktopWidget, QFileDialog,
    QTableWidgetItem, QApplication, QComboBox, QDialog, QMessageBox
)
import copy
from PyQt5.QtGui import QFont,QIcon
from PyQt5.QtCore import QThread, pyqtSignal, QTimer
#from PyQt5.QtCore import Qt
from PyQt5 import QtWidgets
#from PyQt5 import QtCore
#import pathos
import multiprocessing as mp
import queue
import pickle
import plotly.graph_objects as go
import plotly.express as px
#import plotly.io as pio
import networkx as nx
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

class MRMProcessorUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        
        # 设置标题
        self.resize(700,900)
        self.centerWindow()
        
        # 设置参数字典
        self.Data_params={}
        self.Calculate_params={'MS1_Tor':0.000010,'smooth':5,'min_Int':2000,'Points':10}
        self.Manual = False
    
    def centerWindow(self):
        screen = QDesktopWidget().screenGeometry()
        size  = self.geometry()
        LeftValue  = int((screen.width()-size.width())/2)
        TopValue = int((screen.height()-size.height())/2)
        self.move(LeftValue,TopValue)
    
    def set_initial_equal_width(self):
        """设置初始均分宽度"""
        table_width = self.TextBrowser_SampleSelect.width()
        if table_width > 0:
            column_width = table_width // 2
            self.TextBrowser_SampleSelect.setColumnWidth(0, column_width)
            self.TextBrowser_SampleSelect.setColumnWidth(1, column_width)
        else:
            # 如果宽度为0，重试一次
            QTimer.singleShot(50, self.set_initial_equal_width)
            
    def CloseEvent(self,event):
        os._exit(0)
        
    def initUI(self):
        globallayout = QVBoxLayout()
        # Data import
        Data_import_Widget = QWidget()
        Data_import_Layout = QGridLayout()
        # Union
        self.Label_SampleSelect = QLabel('Select data')
        self.Label_Target_Select = QLabel('Select Target List')
        self.Label_Params_Select = QLabel('Set params')
        
        self.TextBrowser_SampleSelect = QTableWidget()
        self.TextBrowser_SampleSelect.setColumnCount(2)
        self.TextBrowser_SampleSelect.setHorizontalHeaderLabels(['Name','Manual List'])
        # 设置表头行为
        header = self.TextBrowser_SampleSelect.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Interactive)  # 允许用户调整
        header.setStretchLastSection(False)  # 禁用最后一列特殊拉伸
        
        # 使用定时器在界面显示后设置初始宽度
        QTimer.singleShot(100, self.set_initial_equal_width)
        #self.TextBrowser_SampleSelect.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        #self.TextBrowser_SampleSelect.horizontalHeader().setStretchLastSection(True)
        self.PushButton_SampleSelect = QPushButton('Select')
        self.PushButton_Target_Select = QPushButton('Select')
        self.TextBrowser_Target_Select = QLabel('*')
        self.PushButton_SampleSelect.setToolTip('Select sample *.mzML documents')
        self.PushButton_Target_Select.setToolTip('Select Target List *.xlsx document')
        # Slot
        self.PushButton_SampleSelect.clicked.connect(self.SelectSampleFile)
        self.PushButton_Target_Select.clicked.connect(self.SelectTargetFile)
        
        Data_import_Layout.addWidget(self.Label_SampleSelect,1,0)
        Data_import_Layout.addWidget(self.TextBrowser_SampleSelect, 2, 0) # row 1， column 0
        Data_import_Layout.addWidget(self.PushButton_SampleSelect, 2, 1) # row 1， column 1
        Data_import_Layout.addWidget(self.Label_Target_Select,5,0)
        Data_import_Layout.addWidget(self.TextBrowser_Target_Select, 6, 0)
        Data_import_Layout.addWidget(self.PushButton_Target_Select, 6, 1) # 行，列，行高，列宽
        Data_import_Widget.setLayout(Data_import_Layout)
        
        # Params setting
        Params_setting_Widget = QWidget()
        Params_setting_Layout = QGridLayout()
        # Union
        self.Lable_MS_Tor = QLabel('MS1 Tolerance(Da)')
        self.Lable_IS_RT_Tor = QLabel('IS RT Tolerance(min)')
        self.Lable_RT_Tor = QLabel('RT Tolerance(min)')
        self.Lable_Int_min = QLabel('Min intensity')
        self.Lable_Point = QLabel('Min scan points')
        self.Lable_SingleData = QLabel('Report Type')
        self.LineEdit_MS_Tor = QLineEdit('0.35')
        self.LineEdit_IS_RT_Tor = QLineEdit('0.5')
        self.LineEdit_RT_Tor = QLineEdit('0.1')
        self.LineEdit_Int_min = QLineEdit('2000')
        self.LineEdit_Point = QLineEdit('10')
        self.Combox_Type = QComboBox()
        self.Combox_Type.addItems(['Detailed Report','Summary Report'])
        self.Combox_Type.setCurrentText('Summary Report')
        
        Params_setting_Layout.addWidget(self.Label_Params_Select,0,0)
        Params_setting_Layout.addWidget(self.Lable_MS_Tor,1,0)
        Params_setting_Layout.addWidget(self.LineEdit_MS_Tor,1,1)
        Params_setting_Layout.addWidget(self.Lable_IS_RT_Tor,1,2)
        Params_setting_Layout.addWidget(self.LineEdit_IS_RT_Tor,1,3)
        Params_setting_Layout.addWidget(self.Lable_RT_Tor,1,4)
        Params_setting_Layout.addWidget(self.LineEdit_RT_Tor,1,5)
        Params_setting_Layout.addWidget(self.Lable_Int_min,2,0)
        Params_setting_Layout.addWidget(self.LineEdit_Int_min,2,1)
        Params_setting_Layout.addWidget(self.Lable_Point,2,2)
        Params_setting_Layout.addWidget(self.LineEdit_Point,2,3)
        Params_setting_Layout.addWidget(self.Lable_SingleData,2,4)
        Params_setting_Layout.addWidget(self.Combox_Type,2,5)
        Params_setting_Widget.setLayout(Params_setting_Layout)
        
        # Run button
        Run_button_Widget = QWidget()
        Run_button_Layout = QGridLayout()
        self.PushButton_Run = QPushButton('Run')
        self.Label_Run1 = QLabel('')
        self.Label_Run2 = QLabel('')
        self.Label_Run3 = QLabel('')
        Run_button_Layout.addWidget(self.Label_Run1,0,0)
        Run_button_Layout.addWidget(self.Label_Run2,0,1)
        Run_button_Layout.addWidget(self.Label_Run3,0,2)
        Run_button_Layout.addWidget(self.PushButton_Run,0,3)
        Run_button_Widget.setLayout(Run_button_Layout)
        self.PushButton_Run.clicked.connect(self.Run)
        
        # statusbar 
        StatusBar_Widget = QWidget()
        StatusBar_Layout = QGridLayout()
        # Union
        self.Label_process = QLabel('Processing bar')
        self.Label_process_sub = QLabel('Step')
        self.process_bar = QProgressBar()
        self.process_bar.setStyleSheet("QProgressBar { border: 2px solid grey; border-radius: 5px; color: rgb(20,20,20);  background-color: #FFFFFF; text-align: center;}QProgressBar::chunk {background-color: rgb(100,200,200); border-radius: 10px; margin: 0.1px;  width: 1px;}")
        font = QFont()
        font.setBold(True)
        font.setWeight(30)
        self.process_bar.setFont(font)
        self.process_bar.setMaximum(100)
        self.process_bar.setMinimum(0)
        self.process_bar.setValue(0)
        StatusBar_Layout.addWidget(self.Label_process,0,0)
        StatusBar_Layout.addWidget(self.Label_process_sub,1,0)
        StatusBar_Layout.addWidget(self.process_bar,1,1)
        StatusBar_Widget.setLayout(StatusBar_Layout)
        
        # Global Layout setting
        globallayout.addWidget(Data_import_Widget)
        globallayout.addWidget(Params_setting_Widget)
        globallayout.addWidget(Run_button_Widget)
        globallayout.addWidget(StatusBar_Widget)
        self.setLayout(globallayout)
        self.setWindowTitle('MRM Processor')
        
        
    def SelectSampleFile(self):
        FileName,FileType = QFileDialog.getOpenFileNames(self,"Select files",os.getcwd(),"mzML Files (*.mzML)")
        if len(FileName) > 0:
            for i in FileName:
                name_begin = i.rfind('/')
                if name_begin == -1:
                    name_begin = i.rfind('\\')
                name_end = len(i)
                LineEdit_Manual= QLineEdit('')
                self.Data_params[i[name_begin+1:name_end-5]] = {'Name':i[name_begin+1:name_end-5],'Path':i,'LineEdit_Manual':LineEdit_Manual}
                self.TextBrowser_SampleSelect.insertRow(self.TextBrowser_SampleSelect.rowCount())
                self.TextBrowser_SampleSelect.setItem(self.TextBrowser_SampleSelect.rowCount()-1,0,QTableWidgetItem(i[name_begin+1:name_end-5]))
                self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,1,self.Data_params[i[name_begin+1:name_end-5]]['LineEdit_Manual'])
                QApplication.processEvents()
        
    def SelectTargetFile(self):
        FileName,FileType = QFileDialog.getOpenFileName(self,"选取文件",os.getcwd(),"Excel Files(*.xlsx)")
        name_begin = FileName.rfind('/')
        if name_begin == -1:
            name_begin = FileName.rfind('\\')
        name_end = len(FileName)
        self.TextBrowser_Target_Select.setText(FileName[name_begin+1:name_end])
        self.TargetPath = FileName
    
    def update_progress(self, value):
        self.process_bar.setValue(value)
        QApplication.processEvents()

    def on_multiprocess_done(self, result_dict):
        self.pool_result = result_dict
        self.All_Data_Area_1 = pd.read_excel(self.TargetPath)
        self.All_Data_Area_1 = self.All_Data_Area_1[self.All_Data_Area_1['IS'].notna()]
        self.All_Data_Area_1 = self.All_Data_Area_1.loc[:, ['Identity keys']]
        self.All_Data_Area_1.reset_index(drop=True,inplace=True)
        self.All_Data_Area_2 = self.All_Data_Area_1.copy()
        self.All_Data_Area_IS_1 = self.All_Data_Area_1.copy()
        self.All_Data_Area_IS_2 = self.All_Data_Area_1.copy()
        self.All_Data_Time = self.All_Data_Area_1.copy()
        for i in self.Data_params.keys():
            self.Data_params[i]['TargetList'] = self.pool_result[i][0].copy()
            self.Data_params[i]['Calibrant'] = self.pool_result[i][1].copy()
            self.All_Data_Area_1[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_1']
            self.All_Data_Area_2[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_2']
            self.All_Data_Area_IS_1[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_1/IS']
            self.All_Data_Area_IS_2[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_2/IS']
            self.All_Data_Time[i] = self.Data_params[i]['TargetList'].loc[:, 'RT_1']
        localtime = time.localtime()
        suffix = f"{localtime.tm_mon}{localtime.tm_mday}{localtime.tm_hour}{localtime.tm_min}"
        basepath = self.filepath_title
        self.All_Data_Area_1.to_excel(basepath + f'MRM_Area_1_{suffix}.xlsx', index=False)
        self.All_Data_Area_2.to_excel(basepath + f'MRM_Area_2_{suffix}.xlsx', index=False)
        self.All_Data_Area_IS_1.to_excel(basepath + f'MRM_Area_IS_1_{suffix}.xlsx', index=False)
        self.All_Data_Area_IS_2.to_excel(basepath + f'MRM_Area_IS_2_{suffix}.xlsx', index=False)
        self.All_Data_Time.to_excel(basepath + f'MRM_Time_{suffix}.xlsx', index=False)
        self.Label_process_sub.setText("Finish")
        self.process_bar.setValue(100)
        self.PushButton_Run.setEnabled(True)
        
    def Run(self):
        self.PushButton_Run.setEnabled(False)
        if self.Manual == False:
            self.process_bar.setValue(0)
            for i in self.Data_params.keys():
                name_begin = self.Data_params[i]['Path'].rfind('/')
                if name_begin > 0:      
                    self.filepath_title = self.Data_params[i]['Path'][0:name_begin]+'/'
                else:
                    name_begin = self.Data_params[i]['Path'].rfind('\\')
                    self.filepath_title = self.Data_params[i]['Path'][0:name_begin]+'\\'
                break
            ''' ---Bench Process MRM Data ---'''
            self.Label_process_sub.setText('Processing')
            self.process_bar.setMaximum(len(self.Data_params.keys()))
            self.params = {
                'MS_Tor': self.LineEdit_MS_Tor.text(),
                'IS_RT_Tor': self.LineEdit_IS_RT_Tor.text(),
                'RT_Tor': self.LineEdit_RT_Tor.text(),
                'Int_min': self.LineEdit_Int_min.text(),
                'Point': self.LineEdit_Point.text(),
                'TargetPath': self.TargetPath,
                'Type': self.Combox_Type.currentText(),
            }
            self.worker = MRMWorker(self.Data_params, self.params)
            self.worker.progress.connect(self.update_progress)
            self.worker.finished.connect(self.on_multiprocess_done)
            self.worker.start()
            self.Manual =True
        else:
            for i in self.Data_params.keys():
                if self.Data_params[i]['LineEdit_Manual'].text() != '':
                    temp_MRM = MRMProcess(self.Data_params[i]['Path'])
                    temp_MRM.set_param('MS1_Tor',float(self.LineEdit_MS_Tor.text()))
                    temp_MRM.set_param('IS_RT_Tor',float(self.LineEdit_IS_RT_Tor.text()*60))
                    temp_MRM.set_param('RT_Tor',float(self.LineEdit_RT_Tor.text()*60))
                    temp_MRM.set_param('min_Int',int(self.LineEdit_Int_min.text()))
                    temp_MRM.set_param('Points',int(self.LineEdit_Point.text()))
                    temp_MRM.load_TargetList(TargetPath=self.Data_params[i]['LineEdit_Manual'].text(),FilePath='')
                    temp_MRM.detect_Peak_MRM(FilePath='')
                    self.Data_params[i]['TargetList'] = temp_MRM.TargetList.copy()
                    self.Data_params[i]['Calibrant'] = temp_MRM.Calibrant.copy()
                    self.All_Data_Area_1[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_1']
                    self.All_Data_Area_2[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_2']
                    self.All_Data_Area_IS_1[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_1/IS']
                    self.All_Data_Area_IS_2[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_2/IS']
                    self.All_Data_Time[i] = self.Data_params[i]['TargetList'].loc[:, 'RT_1']
            localtime = time.localtime()
            suffix = f"{localtime.tm_mon}{localtime.tm_mday}{localtime.tm_hour}{localtime.tm_min}"
            basepath = self.filepath_title
            self.All_Data_Area_1.to_excel(basepath + f'MRM_Area_1_{suffix}.xlsx', index=False)
            self.All_Data_Area_2.to_excel(basepath + f'MRM_Area_2_{suffix}.xlsx', index=False)
            self.All_Data_Area_IS_1.to_excel(basepath + f'MRM_Area_IS_1_{suffix}.xlsx', index=False)
            self.All_Data_Area_IS_2.to_excel(basepath + f'MRM_Area_IS_2_{suffix}.xlsx', index=False)
            self.All_Data_Time.to_excel(basepath + f'MRM_Time_{suffix}.xlsx', index=False)
            self.Label_process_sub.setText("Finish")
            self.process_bar.setValue(100)

class MRMWorker(QThread):
    progress = pyqtSignal(int)       # 当前进度
    finished = pyqtSignal(dict)      # 所有结果完成

    def __init__(self, Data_params, params, parent=None):
        super().__init__(parent)
        self.Data_params = Data_params
        self.params = params

    def run(self):
        total_tasks = len(self.Data_params)
        result_dict = {}

        # 使用 ProcessPoolExecutor 提交任务
        with ProcessPoolExecutor(max_workers=max(1, os.cpu_count()-1)) as executor:
            future_to_key = {}
            for key, info in self.Data_params.items():
                future = executor.submit(
                    pool_MRM,
                    info['Path'],
                    float(self.params['MS_Tor']),
                    float(self.params['IS_RT_Tor']) * 60,
                    float(self.params['RT_Tor']) * 60,
                    int(self.params['Int_min']),
                    int(self.params['Point']),
                    self.params['TargetPath'],
                    self.params['Type']
                )
                future_to_key[future] = key

            tasks_done = 0
            # as_completed 会按完成顺序迭代返回结果
            for future in as_completed(future_to_key):
                key = future_to_key[future]
                try:
                    res = future.result()  # 返回 tuple of DataFrames
                    # 用 copy() 避免 Manager/dict pickle 问题
                    result_dict[key] = (res[0].copy(), res[1].copy())
                except Exception as e:
                    print(f"Error in {key}: {e}")

                tasks_done += 1
                percent = int(tasks_done / total_tasks * 100)
                self.progress.emit(percent)

        # 所有任务完成，发射 finished 信号
        self.finished.emit(result_dict)

def pool_MRM(FilePath,MS1_Tor,IS_RT_Tor,RT_Tor,min_Int,Points,List_Path,Detailed):
    if Detailed == 'Detailed Report':
        folder_path = Path(FilePath[0:-5])
        folder_path.mkdir(exist_ok=True)
        Output_Path = FilePath[0:-5]
    else:
        Output_Path = ''
    temp_MRM = MRMProcess(FilePath)
    temp_MRM.set_param('MS1_Tor',MS1_Tor)
    temp_MRM.set_param('IS_RT_Tor',IS_RT_Tor)
    temp_MRM.set_param('RT_Tor',RT_Tor)
    temp_MRM.set_param('min_Int',min_Int)
    temp_MRM.set_param('Points',Points)
    temp_MRM.load_TargetList(TargetPath=List_Path,FilePath=Output_Path)
    temp_MRM.detect_Peak_MRM(FilePath=Output_Path)
    return copy.deepcopy(temp_MRM.TargetList),copy.deepcopy(temp_MRM.Calibrant)


class MRMProcess(object):   
    def __init__(self,DataPath):
        self.OriginData = pyopenms.MSExperiment()
        self.file_path = DataPath
        ''' 储存文件pyopenms.MzMLFile().store("filtered.mzML", exp) '''
        pyopenms.MzMLFile().load(self.file_path,self.OriginData)
        self.OriginData.sortSpectra(True)
        self.__param = {'MS1_Tor':0.35,'RT_Tor':6,'IS_RT_Tor':30,'Target_RT_Tor':18,'min_Int':2000,'min_IS_Int':10000,'min_MZ':100,'min_RT':60,'min_RT_width':6,'max_Noise':2000,
                        'RI_Tor':0.02,'Deconvolution':False,'FeatureDetectPlot':3,'MergeRule':'Union',
                        'UpDown_gap':10,'saveAutoList':False,'smooth':5,'Points':8,'DeconvolutionSimilarityScore':0.98,
                        'assign MS2':True,'FlowRT':2}
        temp_DF_List = []
        for i in self.OriginData.getChromatograms():
            RT,Int = i.get_peaks()
            temp_DF = pd.DataFrame({'Pre_MZ':i.getPrecursor().getMZ(),'Pro_MZ':i.getProduct().getMZ(), 'Int':[Int], 'RTList':[RT]})
            temp_DF_List.append(temp_DF)
        self.Auto_Peak = pd.concat(temp_DF_List,ignore_index=True)
        self.__param['Flow_RT'] = int(self.__param['min_RT_width']/(1+4*sum(np.diff(self.Auto_Peak.loc[0,'RTList']))/(len(np.diff(self.Auto_Peak.loc[0,'RTList']))+1)))
        
    def load_TargetList(self,TargetPath,FilePath=''):
        self.TargetList = pd.read_excel(TargetPath, engine='openpyxl')
        if 'IS' in self.TargetList.keys():
            print('IS')
            self.Calibrant = self.TargetList[self.TargetList['IS'].isna()]
            self.Calibrant.reset_index(drop=True,inplace=True)
            self.TargetList = self.TargetList[self.TargetList['IS'].notna()]
            self.TargetList.reset_index(drop=True,inplace=True)
            self.Calibrant['RT_1'] = ''
            self.Calibrant['Area_1'] = ''
            self.Calibrant['Hight_1'] = ''
            self.Calibrant['Noise_1'] = ''
            self.Calibrant['Edge_left_1'] = ''
            self.Calibrant['Edge_right_1'] = ''
            self.Calibrant['Noise_1'] = ''
            self.Calibrant['RT_2'] = ''
            self.Calibrant['Area_2'] = ''
            self.Calibrant['Hight_2'] = ''
            self.Calibrant['Noise_2'] = ''
            self.Calibrant['Edge_left_2'] = ''
            self.Calibrant['Edge_right_2'] = ''
            self.Calibrant['RT_Calibrant'] = ''
            bar = Bar('Pick Calibrants', max=len(self.Calibrant))
            for i in range(len(self.Calibrant)):  
                bar.next()
                MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.Calibrant.loc[i,'1_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.Calibrant.loc[i,'1_Q3'])<=self.get_param('MS1_Tor'))]
                if pd.notna(self.Calibrant.loc[i,'2_Q1']):
                    MRM_2 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.Calibrant.loc[i,'2_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.Calibrant.loc[i,'2_Q3'])<=self.get_param('MS1_Tor'))]
                else:
                    MRM_2 = MRM_1
                if len(MRM_1)>0 and len(MRM_2)>0:
                    Origin_RT_List_1 = []
                    Origin_RT_List_2 = []
                    Origin_Int_List_1 = []
                    Origin_Int_List_2 = []
                    for i_MRM in MRM_1.index:
                        if len(Origin_Int_List_1) == 0:
                            Origin_RT_List_1 = list(MRM_1.loc[i_MRM,'RTList'])
                            Origin_Int_List_1 = list(MRM_1.loc[i_MRM,'Int'])
                        else:
                            for ii_MRM in range(len(MRM_1.loc[i_MRM,'RTList'])):
                                if MRM_1.loc[i_MRM,'RTList'][ii_MRM] not in Origin_RT_List_1:
                                    Origin_Int_List_1.insert(bisect.bisect_left(Origin_RT_List_1,MRM_1.loc[i_MRM,'RTList'][ii_MRM]), MRM_1.loc[i_MRM,'Int'][ii_MRM])
                                    Origin_RT_List_1.insert(bisect.bisect_left(Origin_RT_List_1,MRM_1.loc[i_MRM,'RTList'][ii_MRM]), MRM_1.loc[i_MRM,'RTList'][ii_MRM])
                                else:
                                    Origin_Int_List_1[np.where(Origin_RT_List_1==MRM_1.loc[i_MRM,'RTList'][ii_MRM])[0][0]] = max(MRM_1.loc[i_MRM,'Int'][ii_MRM],Origin_Int_List_1[np.where(Origin_RT_List_1==MRM_1.loc[i_MRM,'RTList'][ii_MRM])[0][0]])
                    for i_MRM in MRM_2.index:
                        if len(Origin_Int_List_2) == 0:
                            Origin_RT_List_2 = list(MRM_2.loc[i_MRM,'RTList'])
                            Origin_Int_List_2 = list(MRM_2.loc[i_MRM,'Int'])
                        else:
                            for ii_MRM in range(len(MRM_2.loc[i_MRM,'RTList'])):
                                if MRM_2.loc[i_MRM,'RTList'][ii_MRM] not in Origin_RT_List_2:
                                    Origin_Int_List_2.insert(bisect.bisect_left(Origin_RT_List_2,MRM_2.loc[i_MRM,'RTList'][ii_MRM]), MRM_2.loc[i_MRM,'Int'][ii_MRM])
                                    Origin_RT_List_2.insert(bisect.bisect_left(Origin_RT_List_2,MRM_2.loc[i_MRM,'RTList'][ii_MRM]), MRM_2.loc[i_MRM,'RTList'][ii_MRM])
                                else:
                                    Origin_Int_List_2[np.where(Origin_RT_List_2==MRM_2.loc[i_MRM,'RTList'][ii_MRM])[0][0]] = max(MRM_2.loc[i_MRM,'Int'][ii_MRM],Origin_Int_List_2[np.where(Origin_RT_List_2==MRM_2.loc[i_MRM,'RTList'][ii_MRM])[0][0]])
                    if len(FilePath)>0:
                        fig_1 = go.Figure()
                        fig_1.add_trace(
                            go.Scatter(
                                x=[x/60 for x in Origin_RT_List_1],
                                y=Origin_Int_List_1,
                                mode='lines',
                                name='Raw Data',
                                line={'width':1},
                                ))
                        if pd.notna(self.Calibrant.loc[i,'2_Q1']):
                            fig_2 = go.Figure()
                            fig_2.add_trace(
                                go.Scatter(
                                    x=[x/60 for x in Origin_RT_List_2],
                                    y=Origin_Int_List_2,
                                    mode='lines',
                                    name='Raw Data',
                                    line={'width':1},
                                    ))
                    Noise_1 = self.Calculate_Noise(Origin_RT_List_1,Origin_Int_List_1)
                    Peak_Result_1 = self.MRMPeakDetecter(Origin_RT_List_1,Origin_Int_List_1)
                    (Peak_Begin_1,Peak_End_1,Peak_Top_1,Area_List_1,Int_List_1) = Peak_Result_1
                    Index_1 = [x for x in range(len(Int_List_1)) if Int_List_1[x]>=self.get_param('min_IS_Int')]
                    Peak_Begin_1 = [Peak_Begin_1[x] for x in Index_1]
                    Peak_End_1 = [Peak_End_1[x] for x in Index_1]
                    Peak_Top_1 = [Peak_Top_1[x] for x in Index_1]
                    Area_List_1 = [Area_List_1[x] for x in Index_1]
                    Int_List_1 = [Int_List_1[x] for x in Index_1]
                    RT_List_1 = [Origin_RT_List_1[x] for x in Peak_Top_1]
                    RT_L_1 = [Origin_RT_List_1[x] for x in Peak_Begin_1]
                    RT_R_1 = [Origin_RT_List_1[x] for x in Peak_End_1]
                    Noise_2 = self.Calculate_Noise(Origin_RT_List_2,Origin_Int_List_2)
                    Peak_Result_2 = self.MRMPeakDetecter(Origin_RT_List_2,Origin_Int_List_2)
                    (Peak_Begin_2,Peak_End_2,Peak_Top_2,Area_List_2,Int_List_2) = Peak_Result_2
                    Index_2 = [x for x in range(len(Int_List_2)) if Int_List_2[x]>=self.get_param('min_IS_Int')]
                    Peak_Begin_2 = [Peak_Begin_2[x] for x in Index_2]
                    Peak_End_2 = [Peak_End_2[x] for x in Index_2]
                    Peak_Top_2 = [Peak_Top_2[x] for x in Index_2]
                    Area_List_2 = [Area_List_2[x] for x in Index_2]
                    Int_List_2 = [Int_List_2[x] for x in Index_2]
                    RT_List_2 = [Origin_RT_List_2[x] for x in Peak_Top_2]
                    RT_L_2 = [Origin_RT_List_2[x] for x in Peak_Begin_2]
                    RT_R_2 = [Origin_RT_List_2[x] for x in Peak_End_2]
                    Score_matrix= np.zeros([len(RT_List_1),len(RT_List_2)])
                    for i_1 in range(len(RT_List_1)):
                        for i_2 in range(len(RT_List_2)):
                            if abs(RT_List_1[i_1]-RT_List_2[i_2])<self.get_param('RT_Tor'):
                                Score_matrix[i_1,i_2] = 1-abs(RT_List_1[i_1]-RT_List_2[i_2])/self.get_param('RT_Tor')
                    Sm_1,Sm_2 = linear_sum_assignment(Score_matrix,True)
                    Sm_Filter = [(x,y) for x,y in zip(Sm_1,Sm_2) if Score_matrix[x,y]>0]
                    if len(Sm_Filter) > 0:
                        #Diff_RT = [0.5*abs(self.Calibrant.loc[i,'RT']*120-RT_List_1[x[0]]-RT_List_2[x[1]]) for x in Sm_Filter]
                        Score = [1-(0.5*abs(self.Calibrant.loc[i,'RT']*120-RT_List_1[x[0]]-RT_List_2[x[1]]))/(self.Calibrant.loc[i,'RT']*60)+0.5*(Int_List_1[x[0]]/max(Int_List_1)+Int_List_2[x[1]]/max(Int_List_2)) if 0.5*abs(self.Calibrant.loc[i,'RT']*120-RT_List_1[x[0]]-RT_List_2[x[1]])<self.get_param('IS_RT_Tor') else 0 for x in Sm_Filter]
                        Sm_Place = [x for x in range(len(Score)) if Score[x]==max(Score) and max(Score)>0]
                        if len(Sm_Place) > 0:
                            Sm_Place = Sm_Place[0]
                            self.Calibrant.loc[i,'RT_1'] = np.around(RT_List_1[Sm_Filter[Sm_Place][0]]/60,2)
                            self.Calibrant.loc[i,'Area_1'] = int(Area_List_1[Sm_Filter[Sm_Place][0]])
                            self.Calibrant.loc[i,'Hight_1'] = int(Int_List_1[Sm_Filter[Sm_Place][0]])
                            self.Calibrant.loc[i,'Noise_1'] = int(Noise_1)
                            self.Calibrant.loc[i,'Edge_left_1'] = np.around(RT_L_1[Sm_Filter[Sm_Place][0]]/60,2)
                            self.Calibrant.loc[i,'Edge_right_1'] = np.around(RT_R_1[Sm_Filter[Sm_Place][0]]/60,2)
                            if pd.notna(self.Calibrant.loc[i,'2_Q1']):
                                self.Calibrant.loc[i,'RT_2'] = np.around(RT_List_2[Sm_Filter[Sm_Place][1]]/60,2)
                                self.Calibrant.loc[i,'Area_2'] = int(Area_List_2[Sm_Filter[Sm_Place][1]])
                                self.Calibrant.loc[i,'Hight_2'] = int(Int_List_2[Sm_Filter[Sm_Place][1]])
                                self.Calibrant.loc[i,'Noise_2'] = int(Noise_2)
                                self.Calibrant.loc[i,'Edge_left_2'] = np.around(RT_L_2[Sm_Filter[Sm_Place][1]]/60,2)
                                self.Calibrant.loc[i,'Edge_right_2'] = np.around(RT_R_2[Sm_Filter[Sm_Place][1]]/60,2)
                            self.Calibrant.loc[i,'RT_Calibrant'] = np.around(0.5*(RT_List_1[Sm_Filter[Sm_Place][0]]/60+RT_List_2[Sm_Filter[Sm_Place][1]]/60),2)
                            if len(FilePath)>0:
                                fig_1.add_trace(
                                    go.Scatter(
                                        x=[x/60 for x in Origin_RT_List_1][Peak_Begin_1[Sm_Filter[Sm_Place][0]]:Peak_End_1[Sm_Filter[Sm_Place][0]]+1],
                                        y=Origin_Int_List_1[Peak_Begin_1[Sm_Filter[Sm_Place][0]]:Peak_End_1[Sm_Filter[Sm_Place][0]]+1],
                                        mode='lines',
                                        name=self.TargetList.loc[i,'Identity keys'],
                                        line={'width':1},
                                        fill='tozeroy',
                                        ))
                                
                                if pd.notna(self.Calibrant.loc[i,'2_Q1']):
                                    fig_2.add_trace(
                                        go.Scatter(
                                            x=[x/60 for x in Origin_RT_List_2][Peak_Begin_2[Sm_Filter[Sm_Place][1]]:Peak_End_2[Sm_Filter[Sm_Place][1]]+1],
                                            y=Origin_Int_List_2[Peak_Begin_2[Sm_Filter[Sm_Place][1]]:Peak_End_2[Sm_Filter[Sm_Place][1]]+1],
                                            mode='lines',
                                            name=self.TargetList.loc[i,'Identity keys'],
                                            line={'width':1},
                                            fill='tozeroy',
                                            ))
                    if len(FilePath)>0:
                        fig_1.update_layout(
                            width=600,
                            height=500,
                            titlefont={'size':10},
                            title='Q1:'+str(self.Calibrant.loc[i,'1_Q1'])+' Q3:'+str(self.Calibrant.loc[i,'1_Q3']),
                            plot_bgcolor='rgba(0,0,0,0)',
                            xaxis={'title':{'font':{'size':10},'standoff':0},
                                   'linecolor':'black',
                                   'tickfont':{'size':10,'color':'black'},
                                   'ticks':'outside',
                                   'ticklen':2,
                                   'tickformat':'0.01f', # 统一小数点
                                   },
                            yaxis={'title':{'font':{'size':10},'standoff':0},
                                   'linecolor':'black',
                                   'tickfont':{'size':10,'color':'black'},
                                   #'dtick':10000,
                                   'ticks':'outside',
                                   'ticklen':2,
                                   #'range':(0,20001),
                                   'exponentformat':'e', # 科学记数法
                                   'tickformat':'0.01E', # 统一小数点
                                   },
                            legend_title_text='Gradeint :',
                            legend_traceorder='reversed',
                            legend={
                                'font': {'size':10},
                                'orientation':'v',  # 改为垂直方向
                                'yanchor':'top',    # 锚点改为顶部
                                'y':1,             # 顶部对齐
                                'xanchor':'left',   # 左侧对齐
                                'x':1.02,          # 放在图表右侧外部
                                'bgcolor':'rgba(0,0,0,0)',
                                #'bordercolor':'black',
                                'borderwidth':0,
                                'itemwidth':30,     # 图例项宽度
                            },
                            )
                        fig_1.layout.font.family = 'Helvetica'
                        fig_1.update_layout(showlegend=True)
                        fig_1.update_xaxes(hoverformat='.2f')
                        fig_1.write_html(FilePath+'/'+self.Calibrant.loc[i,'Identity keys']+'-Q1-'+str(self.Calibrant.loc[i,'1_Q1'])+' Q3-'+str(self.Calibrant.loc[i,'1_Q3'])+'.html',config={'responsive': False})
                        if pd.notna(self.Calibrant.loc[i,'2_Q1']):
                            fig_2.update_layout(
                                width=600,
                                height=500,
                                titlefont={'size':10},
                                title='Q1:'+str(self.Calibrant.loc[i,'2_Q1'])+' Q3:'+str(self.Calibrant.loc[i,'2_Q3']),
                                plot_bgcolor='rgba(0,0,0,0)',
                                xaxis={'title':{'font':{'size':10},'standoff':0},
                                       'linecolor':'black',
                                       'tickfont':{'size':10,'color':'black'},
                                       'ticks':'outside',
                                       'ticklen':2,
                                       'tickformat':'0.01f', # 统一小数点
                                       },
                                yaxis={'title':{'font':{'size':10},'standoff':0},
                                       'linecolor':'black',
                                       'tickfont':{'size':10,'color':'black'},
                                       #'dtick':10000,
                                       'ticks':'outside',
                                       'ticklen':2,
                                       #'range':(0,20001),
                                       'exponentformat':'e', # 科学记数法
                                       'tickformat':'0.01E', # 统一小数点
                                       },
                                legend_title_text='Gradeint :',
                                legend_traceorder='reversed',
                                legend={
                                    'font': {'size':10},
                                    'orientation':'v',  # 改为垂直方向
                                    'yanchor':'top',    # 锚点改为顶部
                                    'y':1,             # 顶部对齐
                                    'xanchor':'left',   # 左侧对齐
                                    'x':1.02,          # 放在图表右侧外部
                                    'bgcolor':'rgba(0,0,0,0)',
                                    #'bordercolor':'black',
                                    'borderwidth':0,
                                    'itemwidth':30,     # 图例项宽度
                                },
                                )
                            fig_2.layout.font.family = 'Helvetica'
                            fig_2.update_layout(showlegend=True)
                            fig_2.update_xaxes(hoverformat='.2f')
                            fig_2.write_html(FilePath+'/'+self.Calibrant.loc[i,'Identity keys']+'-Q1-'+str(self.Calibrant.loc[i,'2_Q1'])+' Q3-'+str(self.Calibrant.loc[i,'2_Q3'])+'.html',config={'responsive': False})
                        self.Calibrant.to_excel(FilePath+'/Calibrant-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
                else:
                    Noise_1 = self.Calculate_Noise(MRM_1.iloc[0,:]['RTList'],MRM_1.iloc[0,:]['Int'])
                    Peak_Result_1 = self.MRMPeakDetecter(MRM_1.iloc[0,:]['RTList'],MRM_1.iloc[0,:]['Int'])
                    (Peak_Begin_1,Peak_End_1,Peak_Top_1,Area_List_1,Int_List_1) = Peak_Result_1
                    Index_1 = [x for x in range(len(Int_List_1)) if Int_List_1[x]>=self.get_param('min_Int')]
                    Peak_Begin_1 = [Peak_Begin_1[x] for x in Index_1]
                    Peak_End_1 = [Peak_End_1[x] for x in Index_1]
                    Peak_Top_1 = [Peak_Top_1[x] for x in Index_1]
                    Area_List_1 = [Area_List_1[x] for x in Index_1]
                    Int_List_1 = [Int_List_1[x] for x in Index_1]
                    RT_List_1 = [MRM_1.iloc[0,:]['RTList'][x] for x in Peak_Top_1]
                    RT_L_1 = [MRM_1.iloc[0,:]['RTList'][x] for x in Peak_Begin_1]
                    RT_R_1 = [MRM_1.iloc[0,:]['RTList'][x] for x in Peak_End_1]
                    Score = [1-(abs(self.Calibrant.loc[i,'RT']*60-RT_List_1[x]))/(self.Calibrant.loc[i,'RT']*60)+Int_List_1[x]/max(Int_List_1) if abs(self.Calibrant.loc[i,'RT']*60-RT_List_1[x])<self.get_param('IS_RT_Tor') else 0 for x in range(len(RT_List_1))]
                    Sm_Place = [x for x in range(len(Score)) if Score[x]==max(Score) and max(Score)>0]
                    if len(Sm_Place) > 0:
                        Sm_Place = Sm_Place[0]
                        self.Calibrant.loc[i,'RT_1'] = np.around(RT_List_1[range(len(RT_List_1))[Sm_Place]]/60,2)
                        self.Calibrant.loc[i,'Area_1'] = int(Area_List_1[range(len(RT_List_1))[Sm_Place]])
                        self.Calibrant.loc[i,'Hight_1'] = int(Int_List_1[range(len(RT_List_1))[Sm_Place]])
                        self.Calibrant.loc[i,'Noise_1'] = int(Noise_1)
                        self.Calibrant.loc[i,'Edge_left_1'] = np.around(RT_L_1[range(len(RT_List_1))[Sm_Place]]/60,2)
                        self.Calibrant.loc[i,'Edge_right_1'] = np.around(RT_R_1[range(len(RT_List_1))[Sm_Place]]/60,2)
                        self.Calibrant.loc[i,'RT_Calibrant'] = np.around(RT_List_1[range(len(RT_List_1))[Sm_Place]]/60,2)
                        if len(FilePath)>0:
                            fig_1 = go.Figure()
                            fig_1.add_trace(
                                go.Scatter(
                                    x=list(MRM_1.loc[:,'RTList']/60)[0],
                                    y=list(MRM_1.loc[:,'Int'])[0],
                                    mode='lines',
                                    name='Raw Data',
                                    line={'width':1},))
                            fig_1.add_trace(
                                go.Scatter(
                                    x=list(MRM_1.loc[:,'RTList']/60)[0][Peak_Begin_1[range(len(RT_List_1))[Sm_Place]]:Peak_End_1[range(len(RT_List_1))[Sm_Place]]+1],
                                    y=list(MRM_1.loc[:,'Int'])[0][Peak_Begin_1[range(len(RT_List_1))[Sm_Place]]:Peak_End_1[range(len(RT_List_1))[Sm_Place]]+1],
                                    mode='lines',
                                    name=self.TargetList.loc[i,'Identity keys'],
                                    line={'width':1},
                                    fill='tozeroy',)
                                )
                            fig_1.update_layout(
                                width=600,
                                height=500,
                                titlefont={'size':10},
                                title='Q1:'+str(self.Calibrant.loc[i,'1_Q1'])+' Q3:'+str(self.Calibrant.loc[i,'1_Q3']),
                                plot_bgcolor='rgba(0,0,0,0)',
                                xaxis={'title':{'font':{'size':10},'standoff':0},
                                       'linecolor':'black',
                                       'tickfont':{'size':10,'color':'black'},
                                       'ticks':'outside',
                                       'ticklen':2,
                                       'tickformat':'0.01f', # 统一小数点
                                       },
                                yaxis={'title':{'font':{'size':10},'standoff':0},
                                       'linecolor':'black',
                                       'tickfont':{'size':10,'color':'black'},
                                       #'dtick':10000,
                                       'ticks':'outside',
                                       'ticklen':2,
                                       #'range':(0,20001),
                                       'exponentformat':'e', # 科学记数法
                                       'tickformat':'0.01E', # 统一小数点
                                       },
                                legend_title_text='Gradeint :',
                                legend_traceorder='reversed',
                                legend={
                                    'font': {'size':10},
                                    'orientation':'v',  # 改为垂直方向
                                    'yanchor':'top',    # 锚点改为顶部
                                    'y':1,             # 顶部对齐
                                    'xanchor':'left',   # 左侧对齐
                                    'x':1.02,          # 放在图表右侧外部
                                    'bgcolor':'rgba(0,0,0,0)',
                                    #'bordercolor':'black',
                                    'borderwidth':0,
                                    'itemwidth':30,     # 图例项宽度
                                },
                                )
                            fig_1.layout.font.family = 'Helvetica'
                            fig_1.update_layout(showlegend=True)
                            fig_1.update_xaxes(hoverformat='.2f')
                            fig_1.write_html(FilePath+'/'+self.Calibrant.loc[i,'Identity keys']+'-Q1-'+str(self.Calibrant.loc[i,'1_Q1'])+' Q3-'+str(self.Calibrant.loc[i,'1_Q3'])+'.html',config={'responsive': False})
                if len(FilePath)>0:
                    self.Calibrant.to_excel(FilePath+'/Calibrant-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
        if 'IS' in self.TargetList.keys() and 'Edge_left_1' not in self.TargetList.keys():
            self.Calibrant = self.Calibrant[self.Calibrant['RT_Calibrant'] != '']
            self.Calibrant.reset_index(drop=True,inplace=True)
            self.TargetList.loc[:,'RT From List'] = self.TargetList.loc[:,'RT']
            if len(self.Calibrant) > 0:
                for i in range(len(self.TargetList)):  
                    Diff_RT = self.Calibrant['RT'] - self.TargetList.loc[i,'RT']
                    RT_Right = bisect.bisect(Diff_RT,0)
                    RT_Left = RT_Right-1
                    if RT_Left < 0:
                        RT_Fix = self.TargetList.loc[i,'RT']*self.Calibrant.loc[RT_Right,'RT_Calibrant']/self.Calibrant.loc[RT_Right,'RT']
                    elif RT_Right == len(Diff_RT):
                        RT_Fix = self.TargetList.loc[i,'RT']*self.Calibrant.loc[RT_Right-1,'RT_Calibrant']/self.Calibrant.loc[RT_Right-1,'RT']
                    else:
                        RT_Fix = 0.5*(self.TargetList.loc[i,'RT']*self.Calibrant.loc[RT_Right,'RT_Calibrant']/self.Calibrant.loc[RT_Right,'RT']+self.TargetList.loc[i,'RT']*self.Calibrant.loc[RT_Left,'RT_Calibrant']/self.Calibrant.loc[RT_Left,'RT'])
                    self.TargetList.loc[i,'RT'] = RT_Fix
        MRM_List = MRM_List = [(x,y) for x,y in zip(self.TargetList['1_Q1'],self.TargetList['1_Q3'])]
        MRM_List = list(set(MRM_List))
        TargetList = []
        for i in MRM_List:
            temp_TargetList = self.TargetList[(self.TargetList['1_Q1']==i[0])&(self.TargetList['1_Q3']==i[1])]
            temp_TL_set = pd.DataFrame({'Identity keys':[list(temp_TargetList.loc[temp_TargetList.index,'Identity keys'])],'RT':[list(temp_TargetList.loc[temp_TargetList.index,'RT'])],'1_Q1':i[0],'1_Q3':i[1],'2_Q1':temp_TargetList.loc[temp_TargetList.index[0],'2_Q1'],'2_Q3':temp_TargetList.loc[temp_TargetList.index[0],'2_Q3'],'TL_Index':[list(temp_TargetList.index)]})
            TargetList.append(temp_TL_set)
        self.TargetList_set=pd.concat(TargetList,ignore_index=True)

    def MRMPeakDetecter(self,Auto_RT_List,Auto_Int_List):
        # 单个色谱峰识别模块，由下方另一函数调用
        def GaussSmooth(x):
            if len(x)==5:
                op = x[0]*0.07+x[1]*0.23+x[2]*0.4+x[3]*0.23+x[4]*0.07
            elif len(x)==3:
                op = x[0]*0.17 +x[1]*0.66 +x[2]*0.17
            else:
                op = sum(x)/len(x)
            return op
        Diff_List = np.diff(Auto_Int_List[1:len(Auto_Int_List)-2])
        Diff_List = np.array(list(map(lambda x:GaussSmooth(Diff_List[x-2:x+2+1]) if x in range(2,len(Diff_List)-2) else Diff_List[x],range(len(Diff_List)))))
        FD_List = np.array(list(map(lambda x: MRMProcess.TFFD(x, Auto_Int_List, Auto_RT_List), range(2, len(Auto_Int_List)-2))))
        FD_List = np.array(list(map(lambda x:GaussSmooth(FD_List[x-2:x+2+1]) if x in range(2,len(FD_List)-2) else FD_List[x],range(len(FD_List)))))
        SD_List = np.array(list(map(lambda x: MRMProcess.TFSD(x, Auto_Int_List, Auto_RT_List), range(2, len(Auto_Int_List)-2))))
        SD_List = np.array(list(map(lambda x:GaussSmooth(SD_List[x-2:x+2+1]) if x in range(2,len(SD_List)-2) else SD_List[x],range(len(SD_List)))))
        ABS_FD_List = abs(np.array(FD_List))
        FD_Median = MRMProcess.FD_Line(ABS_FD_List)
        FD_Median_N = FD_Median*(-1)
        Diff_Median = MRMProcess.Diff_Line(Diff_List)
        SD_Median = MRMProcess.SD_Line(SD_List)
        FD_P = list(filter(lambda x: FD_List[x] > FD_Median, range(1, len(FD_List))))
        FD_N = list(filter(lambda x: FD_List[x] < FD_Median_N, range(1, len(FD_List))))
        FD_Change = list(filter(lambda x: FD_List[x-1] > 0 and FD_List[x] < 0, range(1, len(FD_List))))
        Diff_P = list(filter(lambda x: Diff_List[x] > Diff_Median, range(len(Diff_List))))
        Diff_N = list(filter(lambda x: Diff_List[x] < (Diff_Median*-1), range(len(Diff_List))))
        SD_Place = list(filter(lambda x: SD_List[x] < SD_Median and SD_List[x] < 0 and (SD_List[x-1]>SD_List[x] or SD_List[x+1]>SD_List[x]), range(1,len(SD_List)-1)))
        if len(Diff_P) == 0 or len(FD_P) == 0 or len(Diff_N) == 0 or len(FD_N) == 0 or len(SD_Place) == 0:
            return [],[],[],[],[]
        else:
            # 引入人机交互，调整参数？
            Begin_Place,Complet_BP = MRMProcess.Find_FContinuous(FD_List, FD_P, Diff_P, FDMode='P',mergeRule=self.get_param('MergeRule'),FC_Number=self.get_param('FeatureDetectPlot'))
            Begin_Place += 2
            Complet_BP += 2
            End_Place,Complet_EP = MRMProcess.Find_FContinuous(FD_List, FD_N, Diff_N, FDMode='N',mergeRule=self.get_param('MergeRule'),FC_Number=self.get_param('FeatureDetectPlot'))
            End_Place += 2
            Complet_EP += 2
            Peak_Place = MRMProcess.Find_SDChange(FD_Change, SD_Place)  # 峰顶位置
            Peak_Place += 2
        if len(Begin_Place) >= 1 and len(End_Place) >= 1 and Begin_Place[0] < End_Place[-1]:
            Peak_Begin = []
            Peak_End = []
            for i_FD_Seq in range(len(Begin_Place)):
                Peak_Begin.append(Begin_Place[i_FD_Seq])
                for ii_FD_Seq in range(len(End_Place)):
                    if End_Place[ii_FD_Seq]>=Complet_BP[i_FD_Seq]:
                        if i_FD_Seq != len(Begin_Place)-1 and End_Place[ii_FD_Seq]<=Begin_Place[i_FD_Seq+1]+self.get_param('FeatureDetectPlot'):
                            Peak_End.append(End_Place[ii_FD_Seq])
                            break
                        elif i_FD_Seq != len(Begin_Place)-1 and End_Place[ii_FD_Seq]>Begin_Place[i_FD_Seq+1]+self.get_param('FeatureDetectPlot'):
                            Peak_End.append(-1)
                            break
                        else:
                            Peak_End.append(End_Place[ii_FD_Seq])
                            break
            Peak_filter = list(filter(lambda x:Peak_End[x]>0,range(len(Peak_End))))           
            Peak_Begin = list(map(lambda x:Peak_Begin[x],Peak_filter))
            Peak_End = list(map(lambda x:Peak_End[x],Peak_filter))       
            # 删除上升和下降之间差距过大的峰，例如冲顶峰和平峰以及噪声
            gap_del = []
            for i_com in range(len(Peak_Begin)):
                BP_gap = np.where(Begin_Place==Peak_Begin[i_com])[0][0]
                EP_gap = np.where(End_Place==Peak_End[i_com])[0][0]
                if Complet_EP[EP_gap]-Complet_BP[BP_gap]>=self.get_param('UpDown_gap'):
                    gap_del.append(i_com)
            Peak_Begin = np.array(Peak_Begin)
            Peak_Begin = list(np.delete(Peak_Begin,gap_del))
            Peak_End = np.array(Peak_End)
            Peak_End = list(np.delete(Peak_End,gap_del))
            if len(Peak_Begin) != len(Peak_End):
                raise ValueError('Peak len not match')
            Peak_Top = []
            for i_top in range(len(Peak_Begin)-1, -1, -1):
                temp_top = np.where((Peak_Begin[i_top] < Peak_Place) & (Peak_End[i_top] > Peak_Place))[0]
                if len(temp_top) >= 1:
                    top_ran = np.array(range(len(temp_top)))
                    Auto_Int_List = np.array(Auto_Int_List)
                    top_ran = list(filter(lambda x: Auto_Int_List[Peak_Place[temp_top[x]]] == np.max(Auto_Int_List[Peak_Place[temp_top]]), top_ran))
                    Peak_Top.append(Peak_Place[temp_top[top_ran[0]]])
                else:
                    del Peak_Begin[i_top]
                    del Peak_End[i_top]
            Peak_Top = list(map(lambda x: Peak_Top[-x], range(1, len(Peak_Top)+1)))
            unfit = []
            for i_top in range(len(Peak_Begin)):
                if Peak_Begin[i_top]!=0:
                    if Peak_Begin[i_top] < self.get_param('Flow_RT'):
                        begin_ran = np.array(range(self.get_param('Flow_RT')))
                    else:
                        begin_ran = np.array(range(Peak_Begin[i_top]-self.get_param('Flow_RT')+1, Peak_Begin[i_top]+1))
                    temp_be = list(filter(lambda x: Auto_Int_List[x] == np.min(Auto_Int_List[begin_ran]), begin_ran))
                    Peak_Begin[i_top] = temp_be[-1]
                if Peak_End[i_top]!=len(Auto_Int_List)-1:
                    if Peak_End[i_top]+self.get_param('Flow_RT') >= len(Auto_Int_List):
                        end_ran = np.array(range(len(Auto_Int_List)-self.get_param('Flow_RT'), len(Auto_Int_List)))
                    else:
                        end_ran = np.array(range(Peak_End[i_top], Peak_End[i_top]+self.get_param('Flow_RT')))
                    temp_be = list(filter(lambda x: Auto_Int_List[x] == np.min(Auto_Int_List[end_ran]), end_ran))
                    Peak_End[i_top] = temp_be[0]
                top_ran = np.array(range(Peak_Top[i_top]-1, Peak_Top[i_top]+2))
                temp_be = list(filter(lambda x: Auto_Int_List[x] == np.max(Auto_Int_List[top_ran]), top_ran))
                Peak_Top[i_top] = temp_be[0]
                if Peak_Begin[i_top]>=Peak_End[i_top]:
                    unfit.append(i_top)
            if len(unfit)>0:
                Peak_Begin = np.array(Peak_Begin)
                Peak_End = np.array(Peak_End)
                Peak_Top = np.array(Peak_Top)
                Peak_Begin = list(np.delete(Peak_Begin,unfit))
                Peak_End = list(np.delete(Peak_End,unfit))
                Peak_Top = list(np.delete(Peak_Top,unfit))
            for i_edge in range(len(Peak_Begin)):
                if Peak_Begin[i_edge]>=self.get_param('FeatureDetectPlot'):
                    Begin_Filter = np.where(Auto_Int_List[Peak_Begin[i_edge]-self.get_param('FeatureDetectPlot'):Peak_Begin[i_edge]+2*self.get_param('FeatureDetectPlot')+1]==min(Auto_Int_List[Peak_Begin[i_edge]-self.get_param('FeatureDetectPlot'):Peak_Begin[i_edge]+2*self.get_param('FeatureDetectPlot')+1]))[0][-1]
                    if Peak_Begin[i_edge]-self.get_param('FeatureDetectPlot')+Begin_Filter < Peak_Top[i_edge]:
                        Peak_Begin[i_edge] = Peak_Begin[i_edge]-self.get_param('FeatureDetectPlot')+Begin_Filter
                if Peak_End[i_edge]<len(Auto_Int_List)-self.get_param('FeatureDetectPlot'):
                    End_Filter = np.where(Auto_Int_List[Peak_End[i_edge]-2*self.get_param('FeatureDetectPlot'):Peak_End[i_edge]+self.get_param('FeatureDetectPlot')+1]==min(Auto_Int_List[Peak_End[i_edge]-2*self.get_param('FeatureDetectPlot'):Peak_End[i_edge]+self.get_param('FeatureDetectPlot')+1]))[0][0]
                    if Peak_End[i_edge]+End_Filter-2*self.get_param('FeatureDetectPlot')>Peak_Top[i_edge]:
                        Peak_End[i_edge] = Peak_End[i_edge]+End_Filter-2*self.get_param('FeatureDetectPlot')
            Area_List = []
            Int_List = []
            Final = []
            for i_edge in range(len(Peak_Begin)):
                temp_Area = 0
                for ii_edge in range(Peak_Begin[i_edge],Peak_End[i_edge]):
                    temp_Area += 0.5*(Auto_RT_List[ii_edge+1]-Auto_RT_List[ii_edge])*(Auto_Int_List[ii_edge+1]+Auto_Int_List[ii_edge])
                Area_List.append(temp_Area)
                Int_List.append(max(np.array(Auto_Int_List)[Peak_Begin[i_edge]:Peak_End[i_edge]+1]))
                if max(np.array(Auto_Int_List)[Peak_Begin[i_edge]:Peak_End[i_edge]+1]) > self.get_param('min_Int'):
                    Final.append(i_edge)
            Peak_Begin = [Peak_Begin[x] for x in Final]
            Peak_End = [Peak_End[x] for x in Final]
            Peak_Top = [Peak_Top[x] for x in Final]
            Area_List = [Area_List[x] for x in Final]
            Int_List = [Int_List[x] for x in Final]
            return Peak_Begin,Peak_End,Peak_Top,Area_List,Int_List
        else:
            return [],[],[],[],[]
                
    def Calculate_Noise(self,RT_List,Int_List):
        def GaussSmooth(x):
            if len(x)==5:
                op = x[0]*0.07+x[1]*0.23+x[2]*0.4+x[3]*0.23+x[4]*0.07
            elif len(x)==3:
                op = x[0]*0.17 +x[1]*0.66 +x[2]*0.17
            else:
                op = sum(x)/len(x)
            return op
        #self.Final_Peak_Detect['S/N']=self.Final_Peak_Detect['AverageMZ'].apply(lambda x:0)
        Int_List = np.array(Int_List)
        Auto_Int_List = np.array(list(map(lambda x:GaussSmooth(Int_List[x-2:x+3]) if 2<=x<=len(Int_List)-2 else Int_List[x],range(len(Int_List)))))
        if np.log10(max(Int_List)/self.get_param('min_Int'))>=2:
            noise = sum(abs(Int_List[(Int_List>0)&(Int_List<max(Int_List)*0.1)]-Auto_Int_List[(Int_List>0)&(Int_List<max(Int_List)*0.1)]))/(len(Int_List[(Int_List>0)&(Int_List<max(Int_List)*0.1)])+1)
        else:
            noise = sum(abs(Int_List[Int_List>0]-Auto_Int_List[Int_List>0]))/(len(Int_List[Int_List>0])+1)
        return noise
    
    def CosineSimilarity(ExestList, NewList, LOG=False):
        if  len(ExestList) == 0 :
            return 0
        else:
            if len(ExestList) != len(NewList):
                print("CosineSimilarity error", len(ExestList), len(NewList))
                return 0
            if LOG == True:
                temp_ExestList = []
                for i in range(len(ExestList)):
                    if ExestList[i] != 0:
                        temp_ExestList.append(math.log(ExestList[i]))
                    else:
                        temp_ExestList.append(0)
                ExestList = temp_ExestList
                temp_NewList = []
                for ii in range(len(NewList)):
                    if NewList[ii] != 0:
                        temp_NewList.append(math.log(NewList[ii]))
                    else:
                        temp_NewList.append(0)
                NewList = temp_NewList
            ExestList = np.array(ExestList)
            NewList = np.array(NewList)
            CosineSimilarityValue = ExestList.dot(NewList)/(np.linalg.norm(ExestList) * np.linalg.norm(NewList))
            return CosineSimilarityValue
    
    def detect_Peak_MRM(self,FilePath=''):
        def gaussian(x,a=1,b=10,c=2):
            y=a*math.exp((-1*(x-b)**2)/(2*c**2))
            return y
        def vectorized_ROI(temp_low, temp_high, keys, values):
            mask = (keys >= temp_low) & (keys <= temp_high)
            selected_values = values[mask]
            return [item for sublist in selected_values for item in sublist]
        # 色谱峰识别
        if 'Edge_left_1' not in self.TargetList.keys():
            self.TargetList['IS_Area'] = ''
            self.TargetList['RT_1'] = ''
            self.TargetList['Area_1'] = ''
            self.TargetList['Area_1/IS'] = ''
            self.TargetList['Hight_1'] = ''
            self.TargetList['Noise_1'] = ''
            self.TargetList['Edge_left_1'] = ''
            self.TargetList['Edge_right_1'] = ''
            self.TargetList['RT_2'] = ''
            self.TargetList['Area_2'] = ''
            self.TargetList['Area_2/IS'] = ''
            self.TargetList['Hight_2'] = ''
            self.TargetList['Noise_2'] = ''
            self.TargetList['Edge_left_2'] = ''
            self.TargetList['Edge_right_2'] = ''
            self.TargetList['Similarity'] = ''
            bar = Bar('Processing', max=len(self.TargetList_set))
            for i in range(len(self.TargetList_set)):  
                bar.next()
                MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList_set.loc[i,'1_Q1'])<=0.35)&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList_set.loc[i,'1_Q3'])<=0.35)]
                if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                    MRM_2 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList_set.loc[i,'2_Q1'])<=0.35)&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList_set.loc[i,'2_Q3'])<=0.35)]
                else:
                    MRM_2 = MRM_1
                if len(MRM_1)>0 and len(MRM_2)>0:
                    Origin_RT_List_1 = []
                    Origin_RT_List_2 = []
                    Origin_Int_List_1 = []
                    Origin_Int_List_2 = []
                    for i_MRM in MRM_1.index:
                        if len(Origin_Int_List_1) == 0:
                            Origin_RT_List_1 = list(MRM_1.loc[i_MRM,'RTList'])
                            Origin_Int_List_1 = list(MRM_1.loc[i_MRM,'Int'])
                        else:
                            for ii_MRM in range(len(MRM_1.loc[i_MRM,'RTList'])):
                                if MRM_1.loc[i_MRM,'RTList'][ii_MRM] not in Origin_RT_List_1:
                                    Origin_Int_List_1.insert(bisect.bisect_left(Origin_RT_List_1,MRM_1.loc[i_MRM,'RTList'][ii_MRM]), MRM_1.loc[i_MRM,'Int'][ii_MRM])
                                    Origin_RT_List_1.insert(bisect.bisect_left(Origin_RT_List_1,MRM_1.loc[i_MRM,'RTList'][ii_MRM]), MRM_1.loc[i_MRM,'RTList'][ii_MRM])
                                else:
                                    Origin_Int_List_1[np.where(Origin_RT_List_1==MRM_1.loc[i_MRM,'RTList'][ii_MRM])[0][0]] = max(MRM_1.loc[i_MRM,'Int'][ii_MRM],Origin_Int_List_1[np.where(Origin_RT_List_1==MRM_1.loc[i_MRM,'RTList'][ii_MRM])[0][0]])
                    for i_MRM in MRM_2.index:
                        if len(Origin_Int_List_2) == 0:
                            Origin_RT_List_2 = list(MRM_2.loc[i_MRM,'RTList'])
                            Origin_Int_List_2 = list(MRM_2.loc[i_MRM,'Int'])
                        else:
                            for ii_MRM in range(len(MRM_2.loc[i_MRM,'RTList'])):
                                if MRM_2.loc[i_MRM,'RTList'][ii_MRM] not in Origin_RT_List_2:
                                    Origin_Int_List_2.insert(bisect.bisect_left(Origin_RT_List_2,MRM_2.loc[i_MRM,'RTList'][ii_MRM]), MRM_2.loc[i_MRM,'Int'][ii_MRM])
                                    Origin_RT_List_2.insert(bisect.bisect_left(Origin_RT_List_2,MRM_2.loc[i_MRM,'RTList'][ii_MRM]), MRM_2.loc[i_MRM,'RTList'][ii_MRM])
                                else:
                                    Origin_Int_List_2[np.where(Origin_RT_List_2==MRM_2.loc[i_MRM,'RTList'][ii_MRM])[0][0]] = max(MRM_2.loc[i_MRM,'Int'][ii_MRM],Origin_Int_List_2[np.where(Origin_RT_List_2==MRM_2.loc[i_MRM,'RTList'][ii_MRM])[0][0]])
                    Noise_1 = self.Calculate_Noise(Origin_RT_List_1,Origin_Int_List_1)
                    Peak_Result_1 = self.MRMPeakDetecter(Origin_RT_List_1,Origin_Int_List_1)
                    (Peak_Begin_1,Peak_End_1,Peak_Top_1,Area_List_1,Int_List_1) = Peak_Result_1
                    if len(Peak_Begin_1) > 0:
                        Ex_Peak_Begin_1 = [0]+Peak_End_1
                        Ex_Peak_End_1 = Peak_Begin_1+[len(Origin_Int_List_1)]
                        for v in range(len(Ex_Peak_Begin_1)):
                            if Ex_Peak_Begin_1[v]-Ex_Peak_End_1[v] >= self.get_param('Points') and max(Origin_Int_List_1[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]) >= self.get_param('min_Int'):
                                temp_RT = Origin_RT_List_1[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]
                                temp_Int = Origin_Int_List_1[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]
                                Peak_Result_1 = self.MRMPeakDetecter(temp_RT,temp_Int)
                                (Ex_Peak_Begin_1,Ex_Peak_End_1,Ex_Peak_Top_1,Ex_Area_List_1,Ex_Int_List_1) = Peak_Result_1
                                Ex_Time_Begin_1 = [temp_RT[x] for x in Ex_Peak_Begin_1]
                                Ex_Time_End_1 = [temp_RT[x] for x in Ex_Peak_End_1]
                                Ex_Time_Top_1 = [temp_RT[x] for x in Ex_Peak_Top_1]
                                Ex_Peak_Begin_1 = [x for x in range(len(Origin_RT_List_1)) if Origin_RT_List_1[x] in Ex_Time_Begin_1]
                                Ex_Peak_End_1 = [x for x in range(len(Origin_RT_List_1)) if Origin_RT_List_1[x] in Ex_Time_End_1]
                                Ex_Peak_Top_1 = [x for x in range(len(Origin_RT_List_1)) if Origin_RT_List_1[x] in Ex_Time_Top_1]
                                Peak_Begin_1 += Ex_Peak_Begin_1
                                Peak_End_1 += Ex_Peak_End_1
                                Peak_Top_1 += Ex_Peak_Top_1
                                Area_List_1 += Ex_Area_List_1
                                Int_List_1 += Ex_Int_List_1
                    Index_1 = [x for x in range(len(Int_List_1)) if Int_List_1[x]>=self.get_param('min_Int')]
                    Peak_Begin_1 = [Peak_Begin_1[x] for x in Index_1]
                    Peak_End_1 = [Peak_End_1[x] for x in Index_1]
                    Peak_Top_1 = [Peak_Top_1[x] for x in Index_1]
                    Area_List_1 = [Area_List_1[x] for x in Index_1]
                    Int_List_1 = [Int_List_1[x] for x in Index_1]
                    RT_List_1 = [Origin_RT_List_1[x] for x in Peak_Top_1]
                    RT_L_1 = [Origin_RT_List_1[x] for x in Peak_Begin_1]
                    RT_R_1 = [Origin_RT_List_1[x] for x in Peak_End_1]
                    Noise_2 = self.Calculate_Noise(Origin_RT_List_2,Origin_Int_List_2)
                    Peak_Result_2 = self.MRMPeakDetecter(Origin_RT_List_2,Origin_Int_List_2)
                    (Peak_Begin_2,Peak_End_2,Peak_Top_2,Area_List_2,Int_List_2) = Peak_Result_2
                    if len(Peak_Begin_2) > 0:
                        Ex_Peak_Begin_2 = [0]+Peak_End_2
                        Ex_Peak_End_2 = Peak_Begin_2+[len(Origin_Int_List_2)]
                        for v in range(len(Ex_Peak_Begin_2)):
                            if Ex_Peak_Begin_2[v]-Ex_Peak_End_2[v] >= self.get_param('Points') and max(Origin_Int_List_2[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]) >= self.get_param('min_Int'):
                                temp_RT = Origin_RT_List_2[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]
                                temp_Int = Origin_Int_List_2[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]
                                Peak_Result_2 = self.MRMPeakDetecter(temp_RT,temp_Int)
                                (Ex_Peak_Begin_2,Ex_Peak_End_2,Ex_Peak_Top_2,Ex_Area_List_2,Ex_Int_List_2) = Peak_Result_2
                                Ex_Time_Begin_2 = [temp_RT[x] for x in Ex_Peak_Begin_2]
                                Ex_Time_End_2 = [temp_RT[x] for x in Ex_Peak_End_2]
                                Ex_Time_Top_2 = [temp_RT[x] for x in Ex_Peak_Top_2]
                                Ex_Peak_Begin_2 = [x for x in range(len(Origin_RT_List_2)) if Origin_RT_List_2[x] in Ex_Time_Begin_2]
                                Ex_Peak_End_2 = [x for x in range(len(Origin_RT_List_2)) if Origin_RT_List_2[x] in Ex_Time_End_2]
                                Ex_Peak_Top_2 = [x for x in range(len(Origin_RT_List_2)) if Origin_RT_List_2[x] in Ex_Time_Top_2]
                                Peak_Begin_2 += Ex_Peak_Begin_2
                                Peak_End_2 += Ex_Peak_End_2
                                Peak_Top_2 += Ex_Peak_Top_2
                                Area_List_2 += Ex_Area_List_2
                                Int_List_2 += Ex_Int_List_2
                    Index_2 = [x for x in range(len(Int_List_2)) if Int_List_2[x]>=self.get_param('min_Int')]
                    Peak_Begin_2 = [Peak_Begin_2[x] for x in Index_2]
                    Peak_End_2 = [Peak_End_2[x] for x in Index_2]
                    Peak_Top_2 = [Peak_Top_2[x] for x in Index_2]
                    Area_List_2 = [Area_List_2[x] for x in Index_2]
                    Int_List_2 = [Int_List_2[x] for x in Index_2]
                    RT_List_2 = [Origin_RT_List_2[x] for x in Peak_Top_2]
                    RT_L_2 = [Origin_RT_List_2[x] for x in Peak_Begin_2]
                    RT_R_2 = [Origin_RT_List_2[x] for x in Peak_End_2]
                    Score_matrix= np.zeros([len(RT_List_1),len(RT_List_2)])
                    for i_1 in range(len(RT_List_1)):
                        for i_2 in range(len(RT_List_2)):
                            if abs(RT_List_1[i_1]-RT_List_2[i_2])<self.get_param('RT_Tor'):
                                Score_matrix[i_1,i_2] = 1-abs(RT_List_1[i_1]-RT_List_2[i_2])/self.get_param('RT_Tor')
                    Sm_1,Sm_2 = linear_sum_assignment(Score_matrix,True)
                    Sm_Filter = [(x,y) for x,y in zip(Sm_1,Sm_2) if Score_matrix[x,y]>0]
                    
                    Score_matrix= np.zeros([len(Sm_Filter),len(self.TargetList_set.loc[i,'RT'])])
                    for i_1 in range(len(Sm_Filter)):
                        for i_2 in range(len(self.TargetList_set.loc[i,'RT'])):
                            if abs(RT_List_1[Sm_Filter[i_1][0]]-self.TargetList_set.loc[i,'RT'][i_2]*60)< self.get_param('Target_RT_Tor') and abs(RT_List_2[Sm_Filter[i_1][1]]-self.TargetList_set.loc[i,'RT'][i_2]*60)< self.get_param('Target_RT_Tor'):
                                Score_matrix[i_1,i_2] = 2-abs(RT_List_1[Sm_Filter[i_1][0]]-self.TargetList_set.loc[i,'RT'][i_2]*60)/ self.get_param('Target_RT_Tor') - abs(RT_List_2[Sm_Filter[i_1][1]]-self.TargetList_set.loc[i,'RT'][i_2]*60)/ self.get_param('Target_RT_Tor') + 0.5*(Int_List_1[Sm_Filter[i_1][0]]/max(Int_List_1) + Int_List_2[Sm_Filter[i_1][1]]/max(Int_List_2))
                    Sm_1,Sm_2 = linear_sum_assignment(Score_matrix,True)         
                    Sm_Assign = [(x,y) for x,y in zip(Sm_1,Sm_2) if Score_matrix[x,y]>0]
                    for i_1,i_2 in Sm_Assign:
                        IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS']]['Area_1'])
                        if len(IS_Area)>0:
                            IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS']]['Area_1'])[0]
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS_Area'] = int(IS_Area)
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_1/IS'] = np.around(Area_List_1[Sm_Filter[i_1][0]]/IS_Area,4)
                        else:
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_1/IS'] = 'None'
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS_Area'] = 'None'
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'RT_1'] = np.around(RT_List_1[Sm_Filter[i_1][0]]/60,2)
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_1'] = int(Area_List_1[Sm_Filter[i_1][0]])
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Hight_1'] = int(Int_List_1[Sm_Filter[i_1][0]])
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Noise_1'] = int(Noise_1)
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_left_1'] = np.around(RT_L_1[Sm_Filter[i_1][0]]/60,2)
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_right_1'] = np.around(RT_R_1[Sm_Filter[i_1][0]]/60,2)
                        if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'RT_2'] = np.around(RT_List_2[Sm_Filter[i_1][1]]/60,2)
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_2'] = int(Area_List_2[Sm_Filter[i_1][1]])
                            IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS']]['Area_1'])
                            if len(IS_Area)>0:
                                IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS']]['Area_1'])[0]
                                self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_2/IS'] = np.around(Area_List_2[Sm_Filter[i_1][1]]/IS_Area,4)
                            else:
                                self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_2/IS'] = 'None'
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Hight_2'] = int(Int_List_2[Sm_Filter[i_1][1]])
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Noise_2'] = int(Noise_2)
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_left_2'] = np.around(RT_L_2[Sm_Filter[i_1][1]]/60,2)
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_right_2'] = np.around(RT_R_2[Sm_Filter[i_1][1]]/60,2)
                        x_1=list(MRM_1.loc[:,'RTList']/60)[0][Peak_Begin_1[Sm_Filter[i_1][0]]:Peak_End_1[Sm_Filter[i_1][0]]+1]
                        y_1=list(MRM_1.loc[:,'Int'])[0][Peak_Begin_1[Sm_Filter[i_1][0]]:Peak_End_1[Sm_Filter[i_1][0]]+1]
                        x_2=list(MRM_2.loc[:,'RTList']/60)[0][Peak_Begin_2[Sm_Filter[i_1][1]]:Peak_End_2[Sm_Filter[i_1][1]]+1]
                        y_2=list(MRM_2.loc[:,'Int'])[0][Peak_Begin_2[Sm_Filter[i_1][1]]:Peak_End_2[Sm_Filter[i_1][1]]+1]
                        y_1 = [y_1[x] for x in range(len(x_1)) if x_1[x] in x_2]
                        y_2 = [y_2[x] for x in range(len(x_2)) if x_2[x] in x_1]
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Similarity'] = np.around(MRMProcess.CosineSimilarity(y_1, y_2, LOG=False),3)
                    if len(FilePath)>0:
                        fig_1 = go.Figure()
                        fig_1.add_trace(
                            go.Scatter(
                                x=[x/60 for x in Origin_RT_List_1],
                                y=Origin_Int_List_1,
                                mode='lines',
                                name='Raw Data',
                                line={'width':1},
                                ))
                        fig_2 = go.Figure()
                        fig_2.add_trace(
                            go.Scatter(
                                x=[x/60 for x in Origin_RT_List_2],
                                y=Origin_Int_List_2,
                                mode='lines',
                                name='RawData',
                                line={'width':1},
                                ))   
                        for i_1,i_2 in Sm_Assign:
                            fig_1.add_trace(
                                go.Scatter(
                                    x=[x/60 for x in Origin_RT_List_1][Peak_Begin_1[Sm_Filter[i_1][0]]:Peak_End_1[Sm_Filter[i_1][0]]+1],
                                    y=Origin_Int_List_1[Peak_Begin_1[Sm_Filter[i_1][0]]:Peak_End_1[Sm_Filter[i_1][0]]+1],
                                    mode='lines',
                                    name=self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Identity keys'],
                                    line={'width':1},
                                    fill='tozeroy',
                                    ))
                            fig_2.add_trace(
                                go.Scatter(
                                    x=[x/60 for x in Origin_RT_List_2][Peak_Begin_2[Sm_Filter[i_1][1]]:Peak_End_2[Sm_Filter[i_1][1]]+1],
                                    y=Origin_Int_List_2[Peak_Begin_2[Sm_Filter[i_1][1]]:Peak_End_2[Sm_Filter[i_1][1]]+1],
                                    mode='lines',
                                    name=self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Identity keys'],
                                    line={'width':1},
                                    fill='tozeroy',
                                    ))
                        fig_1.update_layout(
                            width=600,
                            height=500,
                            titlefont={'size':10},
                            title='Q1:'+str(self.TargetList_set.loc[i,'1_Q1'])+' Q3:'+str(self.TargetList_set.loc[i,'1_Q3']),
                            plot_bgcolor='rgba(0,0,0,0)',
                            xaxis={'title':{'font':{'size':10},'standoff':0},
                                   'linecolor':'black',
                                   'tickfont':{'size':10,'color':'black'},
                                   'ticks':'outside',
                                   'ticklen':2,
                                   'tickformat':'0.01f', # 统一小数点
                                   },
                            yaxis={'title':{'font':{'size':10},'standoff':0},
                                   'linecolor':'black',
                                   'tickfont':{'size':10,'color':'black'},
                                   #'dtick':10000,
                                   'ticks':'outside',
                                   'ticklen':2,
                                   #'range':(0,20001),
                                   'exponentformat':'e', # 科学记数法
                                   'tickformat':'0.01E', # 统一小数点
                                   },
                            legend_title_text='Gradeint :',
                            legend_traceorder='reversed',
                            legend={
                                'font': {'size':10},
                                'orientation':'v',  # 改为垂直方向
                                'yanchor':'top',    # 锚点改为顶部
                                'y':1,             # 顶部对齐
                                'xanchor':'left',   # 左侧对齐
                                'x':1.02,          # 放在图表右侧外部
                                'bgcolor':'rgba(0,0,0,0)',
                                #'bordercolor':'black',
                                'borderwidth':0,
                                'itemwidth':30,     # 图例项宽度
                            },
                            )
                        fig_1.layout.font.family = 'Helvetica'
                        fig_1.update_layout(showlegend=True)
                        fig_1.update_xaxes(hoverformat='.2f')
                        fig_2.update_layout(
                            width=600,
                            height=500,
                            titlefont={'size':10},
                            title='Q1:'+str(self.TargetList_set.loc[i,'2_Q1'])+' Q3:'+str(self.TargetList_set.loc[i,'2_Q3']),
                            plot_bgcolor='rgba(0,0,0,0)',
                            xaxis={'title':{'font':{'size':10},'standoff':0},
                                   'linecolor':'black',
                                   'tickfont':{'size':10,'color':'black'},
                                   'ticks':'outside',
                                   'ticklen':2,
                                   'tickformat':'0.1f', # 统一小数点
                                   },
                            yaxis={'title':{'font':{'size':10},'standoff':0},
                                   'linecolor':'black',
                                   'tickfont':{'size':10,'color':'black'},
                                   #'dtick':10000,
                                   'ticks':'outside',
                                   'ticklen':2,
                                   #'range':(0,20001),
                                   'exponentformat':'e', # 科学记数法
                                   'tickformat':'0.01E', # 统一小数点
                                   },
                            legend_title_text='Gradeint :',
                            legend_traceorder='reversed',
                            legend={
                                'font': {'size':10},
                                'orientation':'v',  # 改为垂直方向
                                'yanchor':'top',    # 锚点改为顶部
                                'y':1,             # 顶部对齐
                                'xanchor':'left',   # 左侧对齐
                                'x':1.02,          # 放在图表右侧外部
                                'bgcolor':'rgba(0,0,0,0)',
                                #'bordercolor':'black',
                                'borderwidth':0,
                                'itemwidth':30,     # 图例项宽度
                            },
                            )
                        fig_2.layout.font.family = 'Helvetica'
                        fig_2.update_layout(showlegend=True)
                        fig_2.update_xaxes(hoverformat='.2f')
                        fig_1.write_html(FilePath+'/'+'Q1-'+str(self.TargetList_set.loc[i,'1_Q1'])+' Q3-'+str(self.TargetList_set.loc[i,'1_Q3'])+'.html',config={'responsive': False})
                        if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                            fig_2.write_html(FilePath+'/'+'Q1-'+str(self.TargetList_set.loc[i,'2_Q1'])+' Q3-'+str(self.TargetList_set.loc[i,'2_Q3'])+'.html',config={'responsive': False})
            bar.finish()
            if len(FilePath)>0:
                self.TargetList.to_excel(FilePath+'/MRM-Results-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
        else:
            bar = Bar('Processing', max=len(self.TargetList))
            for i in range(len(self.TargetList)):  
                if (type(self.TargetList.loc[i,'Area_1']) == str or np.isnan(self.TargetList.loc[i,'Area_1']) == True) and np.isnan(self.TargetList.loc[i,'Edge_left_1']) == False:
                    self.TargetList.loc[i,'Area_1'] = ''
                    self.TargetList.loc[i,'Hight_1'] = ''
                    self.TargetList.loc[i,'Noise_1'] = ''
                    self.TargetList.loc[i,'Area_2'] = ''
                    self.TargetList.loc[i,'Hight_2'] = ''
                    self.TargetList.loc[i,'Noise_2'] = ''
                    MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList.loc[i,'1_Q1'])<=0.35)&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList.loc[i,'1_Q3'])<=0.35)]
                    if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                        MRM_2 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList_set.loc[i,'2_Q1'])<=0.35)&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList_set.loc[i,'2_Q3'])<=0.35)]
                    else:
                        MRM_2 = MRM_1
                    if len(MRM_1)>0 and len(MRM_2)>0:
                        Noise_1 = self.Calculate_Noise(MRM_1.iloc[0,:]['RTList'],MRM_1.iloc[0,:]['Int'])
                        Noise_2 = self.Calculate_Noise(MRM_2.iloc[0,:]['RTList'],MRM_2.iloc[0,:]['Int'])
                        Index_1 = range(bisect.bisect_left(MRM_1.iloc[0,:]['RTList'],self.TargetList.loc[i,'Edge_left_1']*60),bisect.bisect_right(MRM_1.iloc[0,:]['RTList'],self.TargetList.loc[i,'Edge_right_1']*60))
                        Index_2 = range(bisect.bisect_left(MRM_2.iloc[0,:]['RTList'],self.TargetList.loc[i,'Edge_left_2']*60),bisect.bisect_right(MRM_2.iloc[0,:]['RTList'],self.TargetList.loc[i,'Edge_right_2']*60))
                        Area_1 = 0
                        for ii in Index_1[1:]:
                            Area_1 += 0.5*(MRM_1.iloc[0,:]['Int'][ii-1]+MRM_1.iloc[0,:]['Int'][ii])*(MRM_1.iloc[0,:]['RTList'][ii]+MRM_1.iloc[0,:]['RTList'][ii-1])
                        Area_2 = 0
                        for ii in Index_2[1:]:
                            Area_2 += 0.5*(MRM_2.iloc[0,:]['Int'][ii-1]+MRM_2.iloc[0,:]['Int'][ii])*(MRM_2.iloc[0,:]['RTList'][ii]+MRM_2.iloc[0,:]['RTList'][ii-1])
                        Int_1 = max([MRM_1.iloc[0,:]['Int'][x] for x in Index_1])
                        Int_2 = max([MRM_2.iloc[0,:]['Int'][x] for x in Index_2])
                        RT_1 = [x for x in Index_1 if MRM_1.iloc[0,:]['Int'][x] == Int_1]
                        RT_1 = np.around(MRM_1.iloc[0,:]['RTList'][RT_1[0]]/60,2)
                        RT_2 = [x for x in Index_2 if MRM_2.iloc[0,:]['Int'][x] == Int_2]
                        RT_2 = np.around(MRM_2.iloc[0,:]['RTList'][RT_2[0]]/60,2)
                        IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[i,'IS']]['Area_1'])
                        if len(IS_Area)>0:
                            IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[i,'IS']]['Area_1'])[0]
                            self.TargetList.loc[i,'IS_Area'] = int(IS_Area)
                            self.TargetList.loc[i,'Area_1/IS'] = np.around(Area_1/IS_Area,4)
                        else:
                            self.TargetList.loc[i,'Area_1/IS'] = 'None'
                            self.TargetList.loc[i,'IS_Area'] = 'None'
                        self.TargetList.loc[i,'RT_1'] = RT_1
                        self.TargetList.loc[i,'Area_1'] = int(Area_1)
                        self.TargetList.loc[i,'Hight_1'] = int(Int_1)
                        self.TargetList.loc[i,'Noise_1'] = int(Noise_1)
                        if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                            self.TargetList.loc[i,'RT_2'] = RT_2
                            self.TargetList.loc[i,'Area_2'] = int(Area_2)
                            IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[i,'IS']]['Area_1'])
                            if len(IS_Area)>0:
                                IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[i,'IS']]['Area_1'])[0]
                                self.TargetList.loc[i,'Area_2/IS'] = np.around(Area_1/IS_Area,4)
                            else:
                                self.TargetList.loc[i,'Area_2/IS'] = 'None'
                            self.TargetList.loc[i,'Hight_2'] = int(Int_2)
                            self.TargetList.loc[i,'Noise_2'] = int(Noise_2)
                        if len(FilePath)>0:
                            fig_1 = go.Figure()
                            fig_1.add_trace(
                                go.Scatter(
                                    x=list(MRM_1.loc[:,'RTList']/60)[0],
                                    y=list(MRM_1.loc[:,'Int'])[0],
                                    mode='lines',
                                    name='Raw Data',
                                    line={'width':1},
                                    ))
                            fig_2 = go.Figure()
                            fig_2.add_trace(
                                go.Scatter(
                                    x=list(MRM_2.loc[:,'RTList']/60)[0],
                                    y=list(MRM_2.loc[:,'Int'])[0],
                                    mode='lines',
                                    name='RawData',
                                    line={'width':1},
                                    ))
                            fig_1.add_trace(
                                go.Scatter(
                                    x=[MRM_1.iloc[0,:]['RTList'][x]/60 for x in Index_1],
                                    y=[MRM_1.iloc[0,:]['Int'][x] for x in Index_1],
                                    mode='lines',
                                    name=self.TargetList.loc[i,'Identity keys'],
                                    line={'width':1},
                                    fill='tozeroy',
                                    ))
                            fig_2.add_trace(
                                go.Scatter(
                                    x=[MRM_2.iloc[0,:]['RTList'][x]/60 for x in Index_2],
                                    y=[MRM_2.iloc[0,:]['Int'][x] for x in Index_2],
                                    mode='lines',
                                    name=self.TargetList.loc[i,'Identity keys'],
                                    line={'width':1},
                                    fill='tozeroy',
                                    ))
                            fig_1.update_layout(
                                width=600,
                                height=500,
                                titlefont={'size':10},
                                title='Q1:'+str(self.TargetList.loc[i,'1_Q1'])+' Q3:'+str(self.TargetList.loc[i,'1_Q3']),
                                plot_bgcolor='rgba(0,0,0,0)',
                                xaxis={'title':{'font':{'size':10},'standoff':0},
                                       'linecolor':'black',
                                       'tickfont':{'size':10,'color':'black'},
                                       'ticks':'outside',
                                       'ticklen':2,
                                       'tickformat':'0.01f', # 统一小数点
                                       },
                                yaxis={'title':{'font':{'size':10},'standoff':0},
                                       'linecolor':'black',
                                       'tickfont':{'size':10,'color':'black'},
                                       #'dtick':10000,
                                       'ticks':'outside',
                                       'ticklen':2,
                                       #'range':(0,20001),
                                       'exponentformat':'e', # 科学记数法
                                       'tickformat':'0.01E', # 统一小数点
                                       },
                                legend_title_text='Gradeint :',
                                legend_traceorder='reversed',
                                legend={
                                    'font': {'size':10},
                                    'orientation':'v',  # 改为垂直方向
                                    'yanchor':'top',    # 锚点改为顶部
                                    'y':1,             # 顶部对齐
                                    'xanchor':'left',   # 左侧对齐
                                    'x':1.02,          # 放在图表右侧外部
                                    'bgcolor':'rgba(0,0,0,0)',
                                    #'bordercolor':'black',
                                    'borderwidth':0,
                                    'itemwidth':30,     # 图例项宽度
                                },
                                )
                            fig_1.layout.font.family = 'Helvetica'
                            fig_1.update_layout(showlegend=True)
                            fig_1.update_xaxes(hoverformat='.2f')  # x轴悬浮值显示两位小数
                            fig_2.update_layout(
                                width=600,
                                height=500,
                                titlefont={'size':10},
                                title='Q1:'+str(self.TargetList.loc[i,'2_Q1'])+' Q3:'+str(self.TargetList.loc[i,'2_Q3']),
                                plot_bgcolor='rgba(0,0,0,0)',
                                xaxis={'title':{'font':{'size':10},'standoff':0},
                                       'linecolor':'black',
                                       'tickfont':{'size':10,'color':'black'},
                                       'ticks':'outside',
                                       'ticklen':2,
                                       'tickformat':'0.01f', # 统一小数点
                                       },
                                yaxis={'title':{'font':{'size':10},'standoff':0},
                                       'linecolor':'black',
                                       'tickfont':{'size':10,'color':'black'},
                                       #'dtick':10000,
                                       'ticks':'outside',
                                       'ticklen':2,
                                       #'range':(0,20001),
                                       'exponentformat':'e', # 科学记数法
                                       'tickformat':'0.01E', # 统一小数点
                                       },
                                legend_title_text='Gradeint :',
                                legend_traceorder='reversed',
                                legend={
                                    'font': {'size':10},
                                    'orientation':'v',  # 改为垂直方向
                                    'yanchor':'top',    # 锚点改为顶部
                                    'y':1,             # 顶部对齐
                                    'xanchor':'left',   # 左侧对齐
                                    'x':1.02,          # 放在图表右侧外部
                                    'bgcolor':'rgba(0,0,0,0)',
                                    #'bordercolor':'black',
                                    'borderwidth':0,
                                    'itemwidth':30,     # 图例项宽度
                                },
                                )
                            fig_2.layout.font.family = 'Helvetica'
                            fig_2.update_layout(showlegend=True)
                            fig_2.update_xaxes(hoverformat='.2f')
                            fig_1.write_html(FilePath+'/1-'+self.TargetList.loc[i,'Identity keys']+'.html',config={'responsive': False})
                            fig_2.write_html(FilePath+'/2-'+self.TargetList.loc[i,'Identity keys']+'.html',config={'responsive': False})        
            bar.finish()
            if len(FilePath)>0:
                self.TargetList.to_excel(FilePath+'/MRM-Manual Results-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)   

    def add_0(x):
        x.append(0)
        return x
    def add_Blank(self,DataPath):
        self.BlankData = pyopenms.MSExperiment()
        if DataPath.endswith('mzML'):
            pyopenms.MzMLFile().load(DataPath,self.BlankData)
        elif DataPath.endswith('mzXML'):
            pyopenms.MzXMLFile().load(DataPath,self.BlankData)
        self.BlankData.sortSpectra(True)
        self.Blank_RT_List = np.array([])
        self.Blank_MZ_List = []
        self.Blank_Int_List = []
        for i in self.BlankData:
            if i.getMSLevel()==1:
                MZ_temp, Int_temp = i.get_peaks()
                self.Blank_MZ_List.append(np.around(MZ_temp,5))
                self.Blank_Int_List.append(np.around(Int_temp,0))
                self.Blank_RT_List = np.append(self.Blank_RT_List,i.getRT())
        self.Blank_RT_List = np.around(self.Blank_RT_List,3)
    def add_RT(x, RT):
        x.append(RT)
        return x
    def add_param(self,name,value):
        self.__param[name]=value
        print('set',name,' = ',value)
    
        
    def set_param(self,name,value):
        if name in self.__param:
            self.__param[name]=value
        else:
            print('no such param')
    def get_param(self,key=''):
        if len(key) == 0:
            print(self.__param)
        else:
            return self.__param[key]
    def TFFD(x, Auto_Int_List, Auto_RT_List):
        FirstDerivative = (Auto_Int_List[x+1]*8+Auto_Int_List[x-2]-Auto_Int_List[x-1] * 8-Auto_Int_List[x+2])/((Auto_RT_List[x+2]-Auto_RT_List[x-2])*3)
        return FirstDerivative
    def TFSD(x, Auto_Int_List, Auto_RT_List):
        SecondDerivative = (Auto_Int_List[x+1]*16+Auto_Int_List[x-1]*16-Auto_Int_List[x] * 30-Auto_Int_List[x+2]-Auto_Int_List[x-2])/((Auto_RT_List[x+2]-Auto_RT_List[x-2])*3)
        return SecondDerivative
    def FD_Line(ABS_FD_List):
        ABS_FD_List = list(filter(lambda x: x < max(ABS_FD_List)*0.20 and x > 0, ABS_FD_List))
        FD_Median = np.median(ABS_FD_List)
        return FD_Median
    def SD_Line(SD_List):
        #SD_List = abs(SD_List)
        if len(SD_List) >= 1:
            SD_Median = np.median(SD_List[SD_List<0])
            return SD_Median
        else:
            return 0
    def Diff_Line(Diff_List,limit=0.1):
        Diff_List = list(map(lambda x: abs(x), Diff_List))
        Diff_List = list(filter(lambda x: x < max(Diff_List)*limit and x != 0, Diff_List))
        if len(Diff_List) >= 1:
            Diff_Median = np.median(Diff_List)
        else:
            Diff_Median = 0
        return Diff_Median 
    
    def Find_FContinuous(FD_List, FD_P, Diff_P, FDMode='P' ,FC_Number=2 ,mergeRule='Intersection'):
        if mergeRule == 'Intersection':
            First_List = list(filter(lambda x:x in FD_P,Diff_P))
            if len(First_List)>0:
                if FDMode == 'P':
                    Second_List = [First_List[0]]
                    Second_Score = [First_List[0]]
                    for i_FSB in range(1,len(First_List)):
                        if First_List[i_FSB]-First_List[i_FSB-1]<=FC_Number-1:
                            Second_Score[-1] = First_List[i_FSB]
                        else:
                            Second_List.append(First_List[i_FSB])
                            Second_Score.append(First_List[i_FSB])
                    Second_List = np.array(Second_List)
                    Second_Score = np.array(Second_Score)
                    temp_Place = np.where(Second_Score-Second_List >= FC_Number-1)[0]
                elif FDMode == 'N':
                    Second_List = [First_List[0]]
                    Second_Score = [First_List[0]]
                    for i_FSB in range(1,len(First_List)):
                        if First_List[i_FSB]-First_List[i_FSB-1]<=FC_Number-1:
                            Second_Score[-1] = First_List[i_FSB]
                        else:
                            Second_List.append(First_List[i_FSB])
                            Second_Score.append(First_List[i_FSB])
                    Second_List = np.array(Second_List)
                    Second_Score = np.array(Second_Score)
                    temp_Place = np.where(Second_Score-Second_List >= FC_Number-1)[0]
                    Second_List = Second_List[temp_Place]
                    Second_Score = Second_Score[temp_Place]
                    for i_FSB in range(len(Second_List)):
                        FD_List_short = FD_List[Second_List[i_FSB]:Second_Score[i_FSB]]
                        Int_min = min(FD_List_short)
                        count = 0
                        for ii_FSB in range(np.where(FD_List_short==Int_min)[0][0],len(FD_List_short)):
                            if FD_List_short[ii_FSB] >= Int_min * 0.05:
                                count += 1
                            if count >= 3 and ii_FSB >= FC_Number-1:
                                Second_Score[i_FSB] = Second_List[i_FSB]+ii_FSB
                                break
        elif mergeRule == 'Union':
            First_List = list(set(FD_P+Diff_P))
            First_List.sort()
            if len(First_List)>0:
                if FDMode == 'P':
                    Second_List = [First_List[0]]
                    Second_Score = [First_List[0]]
                    for i_FSB in range(1,len(First_List)):
                        if First_List[i_FSB]-First_List[i_FSB-1]==1:
                            Second_Score[-1] = First_List[i_FSB]
                        else:
                            Second_List.append(First_List[i_FSB])
                            Second_Score.append(First_List[i_FSB])
                    Second_List = np.array(Second_List)
                    Second_Score = np.array(Second_Score)
                    temp_Place = np.where(Second_Score-Second_List >= FC_Number-1)[0]
                elif FDMode == 'N':
                    Second_List = [First_List[0]]
                    Second_Score = [First_List[0]]
                    for i_FSB in range(1,len(First_List)):
                        if First_List[i_FSB]-First_List[i_FSB-1]==1:
                            Second_Score[-1] = First_List[i_FSB]
                        else:
                            Second_List.append(First_List[i_FSB])
                            Second_Score.append(First_List[i_FSB])
                    Second_List = np.array(Second_List)
                    Second_Score = np.array(Second_Score)
                    temp_Place = np.where(Second_Score-Second_List >= FC_Number-1)[0]
                    Second_List = Second_List[temp_Place]
                    Second_Score = Second_Score[temp_Place]
                    for i_FSB in range(len(Second_List)):
                        FD_List_short = FD_List[Second_List[i_FSB]:Second_Score[i_FSB]]
                        Int_min = min(FD_List_short)
                        count = 0
                        for ii_FSB in range(np.where(FD_List_short==Int_min)[0][0],len(FD_List_short)):
                            if FD_List_short[ii_FSB] >= Int_min * 0.05:
                                count += 1
                            if count >= 3 and ii_FSB >= FC_Number-1:
                                Second_Score[i_FSB] = Second_List[i_FSB]+ii_FSB
                                break
            if len(temp_Place) > 0:
                if FDMode == 'P':
                    return Second_List[temp_Place],Second_Score[temp_Place]
                elif FDMode == 'N':
                    return Second_Score,Second_List
            else:
                return np.array([]),np.array([])
        else:
            return np.array([]),np.array([])
    
    def Find_SDChange(FD_Change, SD_Place):
        if len(FD_Change) <= len(SD_Place):
            list_a = np.array(FD_Change)
            list_b = np.array(SD_Place)
        else:
            list_b = np.array(FD_Change)
            list_a = np.array(SD_Place)
        Change_Place = []
        for i in range(len(list_a)):
            temp = np.where(abs(list_a[i]-list_b) <= 2)[0]
            if len(temp) >= 1:
                Change_Place.append(i)
        if len(Change_Place) > 0:
            Change_Place = np.array(Change_Place)
            return list_a[Change_Place]
        else:
            return np.array([])

 


if __name__ == '__main__':
    mp.freeze_support()
    if not QtWidgets.QApplication.instance():
        app = QtWidgets.QApplication(sys.argv)
    else:
        app = QtWidgets.QApplication.instance()
    app.setWindowIcon(QIcon('./MultipleGradientProcessor.ico'))
    app.setQuitOnLastWindowClosed(True)
    main_MRM = MRMProcessorUI()
    main_MRM.show()
    sys.exit(app.exec_())