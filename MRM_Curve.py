#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 11:42:13 2025

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
from sklearn.linear_model import LinearRegression
import json
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
            column_width = table_width // 3
            self.TextBrowser_SampleSelect.setColumnWidth(0, column_width)
            self.TextBrowser_SampleSelect.setColumnWidth(1, column_width)
            self.TextBrowser_SampleSelect.setColumnWidth(2, column_width)
            #self.TextBrowser_SampleSelect.setColumnWidth(3, column_width)
        else:
            # 如果宽度为0，重试一次
            QTimer.singleShot(50, self.set_initial_equal_width)
            
    def closeEvent(self,event):
        try:
            if hasattr(self, "worker") and self.worker is not None and self.worker.isRunning():
                self.worker.request_stop()      # 你要在 MRMWorker 里实现
                self.worker.wait(3000)          # 最多等 3 秒（你可调大）
        except Exception as e:
            print("close cleanup error:", e)
        event.accept()
        QtWidgets.QApplication.quit()
        
    def initUI(self):
        globallayout = QVBoxLayout()
        # Data import
        Data_import_Widget = QWidget()
        Data_import_Layout = QGridLayout()
        # Union
        self.Label_SampleSelect = QLabel('Select data')
        self.Label_Target_Select = QLabel('Select Target List')
        self.Label_Manual_Select = QLabel('Select Manual List')
        self.Label_Params_Select = QLabel('Peak picking params')
        self.Label_Curve = QLabel('Curve params')
        # 设置表头行为
        self.TextBrowser_SampleSelect = QTableWidget()
        self.TextBrowser_SampleSelect.setColumnCount(3)
        self.TextBrowser_SampleSelect.setHorizontalHeaderLabels(['Name','Type','Concentration'])
        header = self.TextBrowser_SampleSelect.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Interactive)  # 允许用户调整
        header.setStretchLastSection(False)  # 禁用最后一列特殊拉伸
        self.set_initial_equal_width()
        # 使用定时器在界面显示后设置初始宽度
        QTimer.singleShot(100, self.set_initial_equal_width)
        #self.TextBrowser_SampleSelect.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        #self.TextBrowser_SampleSelect.horizontalHeader().setStretchLastSection(True)
        self.PushButton_SampleSelect = QPushButton('Select')
        self.PushButton_Target_Select = QPushButton('Select')
        self.TextBrowser_Target_Select = QLabel('')
        self.PushButton_SampleSelect.setToolTip('Select sample *.mzML documents')
        self.PushButton_Target_Select.setToolTip('Select target list *.xlsx document')
        self.TextBrowser_Manual_Select = QLabel('')
        self.PushButton_Manual_Select = QPushButton('Select')
        self.PushButton_Manual_Select.setToolTip('Select manual list *.xlsx document')
        self.PushButton_Manual_Reset = QPushButton('Reset')
        self.PushButton_Manual_Reset.setToolTip('Reset manual list')
        self.ManualListPath = ''
        # Slot
        self.PushButton_SampleSelect.clicked.connect(self.SelectSampleFile)
        self.PushButton_Target_Select.clicked.connect(self.SelectTargetFile)
        self.PushButton_Manual_Select.clicked.connect(self.SelectManualFile)
        self.PushButton_Manual_Reset.clicked.connect(self.ResetManualFile)
        
        Data_import_Layout.addWidget(self.Label_SampleSelect,1,0)
        Data_import_Layout.addWidget(self.TextBrowser_SampleSelect, 2, 0) # row 1， column 0
        Data_import_Layout.addWidget(self.PushButton_SampleSelect, 2, 2) # row 1， column 1
        Data_import_Layout.addWidget(self.Label_Target_Select,5,0)
        Data_import_Layout.addWidget(self.TextBrowser_Target_Select, 6, 0)
        Data_import_Layout.addWidget(self.PushButton_Target_Select, 6, 1) # 行，列，行高，列宽
        Data_import_Layout.addWidget(self.Label_Manual_Select,7,0)
        Data_import_Layout.addWidget(self.TextBrowser_Manual_Select,8,0)
        Data_import_Layout.addWidget(self.PushButton_Manual_Select,8,1)
        Data_import_Layout.addWidget(self.PushButton_Manual_Reset,8,2)
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
        self.Combox_Export = QComboBox()
        self.Combox_Export.addItems(['Detailed Report','Summary Report'])
        self.Combox_Export.setCurrentText('Summary Report')
        
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
        Params_setting_Layout.addWidget(self.Combox_Export,2,5)
        Params_setting_Widget.setLayout(Params_setting_Layout)
        
        # Curve setting
        Curve_setting_Widget = QWidget()
        Cueve_setting_Layout = QGridLayout()
        self.Lable_R2_Tor = QLabel('R²')
        self.Lable_RE_Tor = QLabel('Relative Error(%)')
        self.Lable_Curve_Weight = QLabel('Weighting factor')
        self.LineEdit_R2_Tor = QLineEdit('0.99')
        self.LineEdit_RE_Tor = QLineEdit('20')
        self.Combox_Weight = QComboBox()
        self.Combox_Weight.addItems(['1','1/x','1/x²'])
        self.Combox_Weight.setCurrentText('1')
        Cueve_setting_Layout.addWidget(self.Label_Curve,0,0)
        Cueve_setting_Layout.addWidget(self.Lable_R2_Tor,1,0)
        Cueve_setting_Layout.addWidget(self.LineEdit_R2_Tor,1,1)
        Cueve_setting_Layout.addWidget(self.Lable_RE_Tor,1,2)
        Cueve_setting_Layout.addWidget(self.LineEdit_RE_Tor,1,3)
        Cueve_setting_Layout.addWidget(self.Lable_Curve_Weight,1,4)
        Cueve_setting_Layout.addWidget(self.Combox_Weight,1,5)
        Curve_setting_Widget.setLayout(Cueve_setting_Layout)
        
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
        globallayout.addWidget(Curve_setting_Widget)
        globallayout.addWidget(Run_button_Widget)
        globallayout.addWidget(StatusBar_Widget)
        self.setLayout(globallayout)
        self.setWindowTitle('MRM Processor')
        
    def SelectSampleFile(self):
        FileName,FileType = QFileDialog.getOpenFileNames(self,"Select files",os.getcwd(),"mzML Files (*.mzML)")
        self.FileName = FileName
        if len(FileName) > 0:
            for i in FileName:
                Combox_Type = QComboBox()
                Combox_Type.addItems(['STD','Blank'])
                Combox_Type.setCurrentText('STD')
                Combox_Type.currentIndexChanged.connect(self.SampleTypeChange)
                LineEdit_Concentration = QLineEdit('')
                name_begin = i.rfind('/')
                if name_begin == -1:
                    name_begin = i.rfind('\\')
                name_end = len(i)
                self.Data_params[i[name_begin+1:name_end-5]] = {'Name':i[name_begin+1:name_end-5],'Path':i,'Type':'STD','Group':1,'Index':'','Combox':Combox_Type,'LineEdit_Concentration':LineEdit_Concentration,'rowCount':self.TextBrowser_SampleSelect.rowCount()}
                self.TextBrowser_SampleSelect.insertRow(self.TextBrowser_SampleSelect.rowCount())
                self.TextBrowser_SampleSelect.setItem(self.TextBrowser_SampleSelect.rowCount()-1,0,QTableWidgetItem(i[name_begin+1:name_end-5]))
                self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,1,self.Data_params[i[name_begin+1:name_end-5]]['Combox'])
                self.TextBrowser_SampleSelect.setCellWidget(self.TextBrowser_SampleSelect.rowCount()-1,2,self.Data_params[i[name_begin+1:name_end-5]]['LineEdit_Concentration'])
        
    def SelectTargetFile(self):
        FileName,FileType = QFileDialog.getOpenFileName(self,"Select",os.getcwd(),"Excel Files(*.xlsx)")
        name_begin = FileName.rfind('/')
        if name_begin == -1:
            name_begin = FileName.rfind('\\')
        name_end = len(FileName)
        self.TextBrowser_Target_Select.setText(FileName[name_begin+1:name_end])
        self.TargetPath = FileName
        
    def SelectManualFile(self):
        FileName,FileType = QFileDialog.getOpenFileName(self,"Select",os.getcwd(),"Excel Files(*.xlsx)")
        name_begin = FileName.rfind('/')
        if name_begin == -1:
            name_begin = FileName.rfind('\\')
        name_end = len(FileName)
        self.TextBrowser_Manual_Select.setText(FileName[name_begin+1:name_end])
        self.ManualListPath = FileName
    
    def ResetManualFile(self):
        self.TextBrowser_Manual_Select.setText('')
        self.ManualListPath = ''
    
    def SampleTypeChange(self):
        for i in self.Data_params.keys():
            if self.sender() == self.Data_params[i]['Combox']:
                self.Data_params[i]['Type'] = self.Data_params[i]['Combox'].currentText()
                if self.Data_params[i]['Combox'].currentText() == 'STD':
                    self.Data_params[i]['Type'] = 'STD'
                    self.Data_params[i]['LineEdit_Dilution'] = QLineEdit('')
                    self.Data_params[i]['LineEdit_Dilution'].setReadOnly(False)
                    self.TextBrowser_SampleSelect.setCellWidget(self.Data_params[i]['rowCount'],2,self.Data_params[i]['LineEdit_Dilution'])
                elif self.Data_params[i]['Combox'].currentText() == 'Blank':
                    self.Data_params[i]['Type'] = 'Blank'
                    self.Data_params[i]['LineEdit_Dilution'] = QLineEdit('Blank')
                    self.Data_params[i]['LineEdit_Dilution'].setReadOnly(True)
                    self.TextBrowser_SampleSelect.setCellWidget(self.Data_params[i]['rowCount'],2,self.Data_params[i]['LineEdit_Dilution'])
            QApplication.processEvents()
    
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
        self.All_Data_Manual = pd.read_excel(self.TargetPath)
        self.All_Data_Manual = self.All_Data_Manual.loc[:,['Identity keys','Polarity','IS','1_Q1','1_Q3']]
        for i in self.Data_params.keys():
            self.Data_params[i]['TargetList'] = self.pool_result[i][0].copy()
            self.Data_params[i]['Calibrant'] = self.pool_result[i][1].copy()
            self.All_Data_Area_1[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_1']
            self.All_Data_Area_2[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_2']
            self.All_Data_Area_IS_1[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_1/IS']
            self.All_Data_Area_IS_2[i] = self.Data_params[i]['TargetList'].loc[:, 'Area_2/IS']
            self.All_Data_Time[i] = self.Data_params[i]['TargetList'].loc[:, 'RT_1']
            self.All_Data_Manual[i+'-L'] = ''
            self.All_Data_Manual[i+'-R'] = ''
            for ii in range(len(self.Data_params[i]['TargetList'])):
                if type(self.Data_params[i]['TargetList'].loc[ii,'Edge_left_1']) != str and self.Data_params[i]['TargetList'].loc[ii,'Edge_left_1'] > 0:
                    Name_Match = [x for x in range(len(self.All_Data_Manual)) if self.All_Data_Manual.loc[x,'Identity keys'] == self.Data_params[i]['TargetList'].loc[ii,'Identity keys']][0]
                    self.All_Data_Manual.loc[Name_Match,i+'-L'] = self.Data_params[i]['TargetList'].loc[ii,'Edge_left_1']
                    self.All_Data_Manual.loc[Name_Match,i+'-R'] = self.Data_params[i]['TargetList'].loc[ii,'Edge_right_1']
            for ii in range(len(self.Data_params[i]['Calibrant'])):
                if type(self.Data_params[i]['Calibrant'].loc[ii,'Edge_left_1']) != str and self.Data_params[i]['Calibrant'].loc[ii,'Edge_left_1'] > 0:
                    Name_Match = [x for x in range(len(self.All_Data_Manual)) if self.All_Data_Manual.loc[x,'Identity keys'] == self.Data_params[i]['Calibrant'].loc[ii,'Identity keys']][0]
                    self.All_Data_Manual.loc[Name_Match,i+'-L'] = self.Data_params[i]['Calibrant'].loc[ii,'Edge_left_1']
                    self.All_Data_Manual.loc[Name_Match,i+'-R'] = self.Data_params[i]['Calibrant'].loc[ii,'Edge_right_1']
        localtime = time.localtime()
        suffix = f"{localtime.tm_mon}{localtime.tm_mday}{localtime.tm_hour}{localtime.tm_min}"
        basepath = self.filepath_title
        self.All_Data_Area_1.to_excel(basepath + f'MRM_Area_1_{suffix}.xlsx', index=False)
        self.All_Data_Area_IS_1.to_excel(basepath + f'MRM_Area_IS_1_{suffix}.xlsx', index=False)
        self.All_Data_Time.to_excel(basepath + f'MRM_Time_{suffix}.xlsx', index=False)
        if self.ManualListPath == '':
            self.All_Data_Area_2.to_excel(basepath + f'MRM_Area_2_{suffix}.xlsx', index=False)
            self.All_Data_Area_IS_2.to_excel(basepath + f'MRM_Area_IS_2_{suffix}.xlsx', index=False)
            self.All_Data_Manual.to_excel(basepath + f'MRM_RT_range_{suffix}.xlsx', index=False)
        Concentration = []
        STD_key = []
        Blank_key = []
        for i in self.Data_params.keys():
            if self.Data_params[i]['Type'] == 'STD':
                Concentration.append(float(self.Data_params[i]['LineEdit_Concentration'].text()))
                STD_key.append(i)
            else:
                Blank_key.append(i)
        Concentration_set = list(set(Concentration))
        Concentration_set.sort()
        results = {}
        results_DF = pd.DataFrame(columns=['Name','slope','intercept','r2','RE','X_L','X_R'])
        for i in self.All_Data_Area_IS_1.index:
            score = []
            i_L = Concentration_set[0]
            Warm_RE = ''
            for i_R in Concentration_set[1:]:
                X = []
                Y = []
                for ii in self.Data_params.keys():
                    if type(self.All_Data_Area_IS_1.loc[i,ii]) != str and self.All_Data_Area_IS_1.loc[i,ii]>0:
                        if i_L<=float(self.Data_params[ii]['LineEdit_Concentration'].text())<=i_R:
                            X.append(float(self.Data_params[ii]['LineEdit_Concentration'].text()))
                            Y.append(self.All_Data_Area_IS_1.loc[i,ii])
                X_set = list(set(X)) 
                X_set.sort()
                if len(X_set) >= 3:
                    Y_mean = []
                    SD = []
                    RSD = []
                    X_set_temp = []
                    for ii in X_set:
                        if len([Y[x] for x in range(len(X)) if X[x]==ii]) == len([x for x in X if x == ii]):
                            y = np.mean([Y[x] for x in range(len(X)) if X[x]==ii])
                            Y_mean.append(y)
                            y_SD = np.std([Y[x] for x in range(len(X)) if X[x]==ii])
                            SD.append(y_SD)
                            RSD.append(y_SD/y)
                            X_set_temp.append(ii)
                    X_set = [[x] for x in X_set_temp]
                    Y_mean = [[x] for x in Y_mean]
                    if self.Combox_Weight.currentText() == '1/x²':
                        w = 1/(np.array(X_set).ravel()**2)
                    elif self.Combox_Weight.currentText() == '1/x':
                        w = 1/(np.array(X_set).ravel())
                    else:
                        w = np.array([1]*len(X_set))
                    model = LinearRegression().fit(X_set, Y_mean, sample_weight=w)
                    slope, intercept, r2 = model.coef_[0], model.intercept_, model.score(X_set, Y_mean)
                    RE = [np.around(((Y_mean[x][0]-intercept[0])/slope[0]-X_set[x][0])/X_set[x][0],4) for x in range(len(X_set))]
                    if r2 >= float(self.LineEdit_R2_Tor.text()) and max(abs(np.array(RE))) <= float(self.LineEdit_RE_Tor.text())/100:
                        score.append((X_set[0],X_set[-1],X_set[-1][0]/X_set[0][0]))
            if len(score) > 0:
                X_all = []
                Y_all = []
                for ii in self.Data_params.keys():
                    if type(self.All_Data_Area_IS_1.loc[i,ii]) != str and self.All_Data_Area_IS_1.loc[i,ii]>0:
                        X_all.append(float(self.Data_params[ii]['LineEdit_Concentration'].text()))
                        Y_all.append(self.All_Data_Area_IS_1.loc[i,ii])
                X_all_set = list(set(X_all)) 
                X_all_set.sort()
                Y_all_mean = []
                SD = []
                RSD = []
                X_all_set_temp = []
                for ii in X_all_set:
                    if len([Y_all[x] for x in range(len(X_all)) if X_all[x]==ii]) == len([x for x in X_all if x == ii]):
                        y = np.mean([Y_all[x] for x in range(len(X_all)) if X_all[x]==ii])
                        Y_all_mean.append(y)
                        y_SD = np.std([Y_all[x] for x in range(len(X_all)) if X_all[x]==ii])
                        SD.append(y_SD)
                        RSD.append(y_SD/y)
                        X_all_set_temp.append(ii)
                X_all_set = X_all_set_temp
                score_len = [x[2] for x in score]
                score = [score[x] for x in range(len(score_len)) if score_len[x] == max(score_len)]
                score_min = [x[0] for x in score]
                score = [score[x] for x in range(len(score_min)) if score_min[x] == min(score_min)][0]
                X = []
                Y = []
                for ii in self.Data_params.keys():
                    if type(self.All_Data_Area_IS_1.loc[i,ii]) != str and self.All_Data_Area_IS_1.loc[i,ii]>0:
                        if score[0][0]<=float(self.Data_params[ii]['LineEdit_Concentration'].text())<=score[1][0]:
                            X.append(float(self.Data_params[ii]['LineEdit_Concentration'].text()))
                            Y.append(self.All_Data_Area_IS_1.loc[i,ii])
                X_set = list(set(X)) 
                X_set.sort()
                Y_all_mean_L = [Y_all_mean[x] for x in range(len(X_all_set)) if X_all_set[x] <= max(X_set)]
                X_all_set_L = [X_all_set[x] for x in range(len(X_all_set)) if X_all_set[x] <= max(X_set)]
                SD_L = [SD[x] for x in range(len(X_all_set)) if X_all_set[x] <= max(X_set)]
                Y_all_mean_R = [Y_all_mean[x] for x in range(len(X_all_set)) if X_all_set[x] > max(X_set)]
                X_all_set_R = [X_all_set[x] for x in range(len(X_all_set)) if X_all_set[x] > max(X_set)]
                SD_R = [SD[x] for x in range(len(X_all_set)) if X_all_set[x] > max(X_set)]
                X_all_set = [[x] for x in X_all_set]
                Y_all_mean = [[x] for x in Y_all_mean]
                Y_mean = []
                X_set_temp = []
                for ii in X_set:
                    if len([Y[x] for x in range(len(X)) if X[x]==ii]) == len([x for x in X if x == ii]):
                        y = np.mean([Y[x] for x in range(len(X)) if X[x]==ii])
                        Y_mean.append(y)
                        X_set_temp.append(ii)
                X_set = [[x] for x in X_set_temp]
                Y_mean = [[x] for x in Y_mean]
                if self.Combox_Weight.currentText() == '1/x²':
                    w = 1/(np.array(X_set).ravel()**2)
                elif self.Combox_Weight.currentText() == '1/x':
                    w = 1/(np.array(X_set).ravel())
                else:
                    w = np.array([1]*len(X_set))
                model = LinearRegression().fit(X_set, Y_mean, sample_weight=w)
                slope, intercept, r2 = model.coef_[0], model.intercept_, model.score(X_set, Y_mean)
                RE = [((Y_mean[x][0]-intercept[0])/slope[0]-X_set[x][0])/X_set[x][0] for x in range(len(X_set))]
                results[self.All_Data_Area_IS_1.loc[i,'Identity keys']] = {"slope": slope[0], "intercept": intercept[0], "r2": r2, 'RE':RE, 'X_set':X_set}
                temp_DF = pd.DataFrame({"Name":self.All_Data_Area_IS_1.loc[i,'Identity keys'],"slope": slope[0], "intercept": intercept[0], "r2": r2, 'RE':[RE], 'X_L':X_set[0],'X_R':X_set[-1]})
                results_DF = pd.concat([results_DF,temp_DF])
                fig = go.Figure()
                fig.add_trace(
                    go.Scatter(
                        x=[x[0] for x in X_all_set],
                        y=[x[0] for x in list(model.predict(X_all_set))],# 绑定额外数据
                        mode='lines',
                        name=self.All_Data_Area_IS_1.loc[i,'Identity keys'],
                        line={'width':2,'color':'rgba(247,174,177,1)'},
                        ))
                fig.add_trace(
                    go.Scatter(
                        x=X_all_set_L,
                        y=Y_all_mean_L,# 绑定额外数据
                        hovertemplate=(
                        "Concentration: %{x}<br>"
                        "Area/IS: %{y}<br>"
                        "<extra></extra>"),
                        mode='markers',
                        name=self.All_Data_Area_IS_1.loc[i,'Identity keys']+' dot',
                        legendgroup=self.All_Data_Area_IS_1.loc[i,'Identity keys']+' dot',
                        showlegend=True,
                        line={'width':2},
                        marker=dict(
                            size=5,
                            color='rgba(247,174,177,0.7)',
                            line=dict(width=1, color="black")
                                ),
                        error_y=dict(
                            type='data',
                            array=SD_L,      # 上下误差 = SD
                            visible=True)
                        ))
                fig.add_trace(
                    go.Scatter(
                        x=X_all_set_R,
                        y=Y_all_mean_R,# 绑定额外数据
                        hovertemplate=(
                        "Concentration: %{x}<br>"
                        "Area/IS: %{y}<br>"
                        "<extra></extra>"),
                        mode='markers',
                        name=self.All_Data_Area_IS_1.loc[i,'Identity keys']+' dot',
                        legendgroup=self.All_Data_Area_IS_1.loc[i,'Identity keys']+' dot',
                        showlegend=False,
                        line={'width':2},
                        marker=dict(
                            size=5,
                            color='rgba(31, 119, 180,0.7)',
                            line=dict(width=1, color="black")
                                ),
                        error_y=dict(
                            type='data',
                            array=SD_R,      # 上下误差 = SD
                            visible=True)
                        ))
                fig.update_layout(
                    legend={
                        'xanchor': 'right',
                        'x': 1.3,  # 改为 0.95，让图例在画布内部（95% 宽度处）
                        'y': 0.5,   # 垂直居中
                        'bgcolor': 'rgba(255,255,255,0.3)',  # 可选：添加半透明背景
                    },
                    #margin=dict(l=50, r=150, b=50, t=50, pad=4),
                    autosize=False,
                    width=700,
                    height=600,
                    title='Curve of '+self.All_Data_Area_IS_1.loc[i,'Identity keys'],
                    titlefont={'size':20},
                    plot_bgcolor='rgba(0,0,0,0)',
                    xaxis={'title':{'text':'Concentration','font':{'size':15},'standoff':0},
                           'linecolor':'black',
                           'tickfont':{'size':11},
                           'ticks':'outside',
                           'ticklen':2,
                           # 添加range设置，在两端留出空白
                           #'range':[-0.5, len(x_List)-0.5]  # 这里-0.5和+0.5表示两端各留出相当于半个类别的空白
                           },
                    yaxis={'title':{'text':'Peak area ratio','font':{'size':15},'standoff':0},
                           'linecolor':'black',
                           'tickfont':{'size':11},
                           'ticks':'outside',
                           'ticklen':2,
                           'side':'left',
                           #'exponentformat':'e', # 科学记数法
                           #'tickformat':'0.01E', # 统一小数点
                           #'range':(0,1500),
                           }
                    )
                fig.layout.font.family = 'Helvetica'
                fig.write_html(self.filepath_title+'Curve of '+self.All_Data_Area_IS_1.loc[i,'Identity keys']+' - '+Warm_RE+f"R²={r2:.4f}"+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.html',config={'responsive': False}) 
                try:
                    #if self.ManualListPath == '':
                        fig = go.Figure()
                        for ii in self.Data_params.keys():
                            if type(self.Data_params[ii]['TargetList'].loc[i,'Plot_RT_1']) != str:
                                fig.add_trace(
                                    go.Scatter(
                                        x = [x/60 for x in self.Data_params[ii]['TargetList'].loc[i,'Plot_RT_1']],
                                        y = self.Data_params[ii]['TargetList'].loc[i,'Plot_Int_1'],
                                        hovertemplate=(
                                        "RT: %{x}<br>"
                                        "Int: %{y}<br>"
                                        "<extra></extra>"),
                                        mode='lines',
                                        name=ii,
                                        line={'width':2,'color':'rgba(247,174,177,1)'},
                                        hoverinfo="skip",
                                        showlegend=False,
                                        legendgroup=ii,
                                        fill="tozeroy",
                                        fillcolor="rgba(31, 119, 180, 0.25)",
                                        ))
                                fig.add_trace(
                                    go.Scatter(
                                        x = [x/60 for x in self.Data_params[ii]['TargetList'].loc[i,'Plot_RT_All_1']],
                                        y = self.Data_params[ii]['TargetList'].loc[i,'Plot_Int_All_1'],
                                        hovertemplate=(
                                        "RT: %{x}<br>"
                                        "Int: %{y}<br>"
                                        "<extra></extra>"),
                                        mode='lines',
                                        name=ii,
                                        line={'width':2,'color':'rgba(247,174,177,1)'},
                                        showlegend=True,
                                        legendgroup=ii,
                                        ))
                        fig.update_layout(
                            #margin=dict(l=50, r=150, b=50, t=50, pad=4),
                            autosize=False,
                            width=700,
                            height=600,
                            title='Curve of '+self.All_Data_Area_IS_1.loc[i,'Identity keys'],
                            titlefont={'size':20},
                            plot_bgcolor='rgba(0,0,0,0)',
                            xaxis={'title':{'text':'Concentration','font':{'size':15},'standoff':0},
                                   'linecolor':'black',
                                   'tickfont':{'size':11},
                                   'ticks':'outside',
                                   'ticklen':2,
                                   # 添加range设置，在两端留出空白
                                   #'range':[-0.5, len(x_List)-0.5]  # 这里-0.5和+0.5表示两端各留出相当于半个类别的空白
                                   },
                            yaxis={'title':{'text':'Peak area ratio','font':{'size':15},'standoff':0},
                                   'linecolor':'black',
                                   'tickfont':{'size':11},
                                   'ticks':'outside',
                                   'ticklen':2,
                                   'side':'left',
                                   #'exponentformat':'e', # 科学记数法
                                   #'tickformat':'0.01E', # 统一小数点
                                   #'range':(0,1500),
                                   }
                            )
                        fig.layout.font.family = 'Helvetica'
                        fig.write_html(self.filepath_title+'Area of '+self.All_Data_Area_IS_1.loc[i,'Identity keys']+'-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.html',config={'responsive': False})
                except Exception:
                    print('Area Draw error')
            else:
                results[self.All_Data_Area_IS_1.loc[i,'Identity keys']] = {"slope": '', "intercept": '', "r2": '', 'RE':[''], 'X_set':['']}
                temp_DF = pd.DataFrame({"Identity keys":self.All_Data_Area_IS_1.loc[i,'Identity keys'],"slope": slope[0], "intercept": intercept[0], "r2": r2, 'RE':[RE], 'X_L':X_set[0],'X_R':X_set[-1]})
                results_DF = pd.concat([results_DF,temp_DF])
        with open(self.filepath_title+'Curve Data-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.json', 'w') as f:
            json.dump(results, f)
        results_DF['Regression equation'] = ''
        for i in range(len(results_DF)):
            if results_DF.loc[i,'slope'] != '':
                if results_DF.loc[i,'intercept'] >= 0:
                    results_DF.loc[i,'Regression equation'] = 'y = ' + str(np.around(results_DF.loc[i,'slope'],4)) + 'x + ' + str(np.around(results_DF.loc[i,'intercept'],4))
                else:
                    results_DF.loc[i,'Regression equation'] = 'y = ' + str(np.around(results_DF.loc[i,'slope'],4)) + 'x - ' + str(np.around(abs(results_DF.loc[i,'intercept']),4))
        results_DF['Regression equation'] = results_DF.apply(lambda row: 'y = '+str(np.around(row['slope'],4))+'')
        results_DF.to_excel(self.filepath_title+'Curve Data-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx')
        self.Label_process_sub.setText("Finish")
        self.process_bar.setValue(100)
        self.PushButton_Run.setEnabled(True)
        
    def Run(self):
        self.PushButton_Run.setEnabled(False)
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
            'Export': self.Combox_Export.currentText(),
            'ManualList':self.ManualListPath
        }
        self.worker = MRMWorker(self.Data_params, self.params)
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.on_multiprocess_done)
        self.worker.start()

class MRMWorker(QThread):
    progress = pyqtSignal(int)       # 当前进度
    finished = pyqtSignal(dict)      # 所有结果完成

    def __init__(self, Data_params, params, parent=None):
        super().__init__(parent)
        self.Data_params = Data_params
        self.params = params

    def request_stop(self):
        self._stop = True
        # 尽量取消还没开始的任务
        for f in self.futures:
            try:
                f.cancel()
            except Exception:
                pass
        # 让进程池尽快关闭（会取消未开始的 future）
        if self.executor is not None:
            try:
                self.executor.shutdown(wait=False, cancel_futures=True)
            except TypeError:
                # 旧 Python 版本没有 cancel_futures 参数
                self.executor.shutdown(wait=False)

    def run(self):
        total_tasks = len(self.Data_params)
        result_dict = {}
        # 使用 ProcessPoolExecutor 提交任务
        with ProcessPoolExecutor(max_workers=max(1, os.cpu_count()-1)) as executor:
            future_to_key = {}
            for key, info in self.Data_params.items():
                if self.params['ManualList'] == '':
                    future = executor.submit(
                        pool_MRM,
                        info['Path'],
                        float(self.params['MS_Tor']),
                        float(self.params['IS_RT_Tor']) * 60,
                        float(self.params['RT_Tor']) * 60,
                        int(self.params['Int_min']),
                        int(self.params['Point']),
                        self.params['TargetPath'],
                        self.params['Export'],
                        info['Type'],
                        False
                    )
                else:
                    future = executor.submit(
                        pool_MRM,
                        info['Path'],
                        float(self.params['MS_Tor']),
                        float(self.params['IS_RT_Tor']) * 60,
                        float(self.params['RT_Tor']) * 60,
                        int(self.params['Int_min']),
                        int(self.params['Point']),
                        self.params['ManualList'],
                        self.params['Export'],
                        info['Type'],
                        True
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

def pool_MRM(FilePath,MS1_Tor,IS_RT_Tor,RT_Tor,min_Int,Points,List_Path,Detailed,Sample_Type,Manual=False):
    if Detailed == 'Detailed Report':
        folder_path = Path(FilePath[0:-5])
        folder_path.mkdir(exist_ok=True)
        Output_Path = FilePath[0:-5]
    else:
        Output_Path = ''
    if Manual==False:
        if Sample_Type != 'Blank':
            temp_MRM = MRMProcess(FilePath)
            temp_MRM.set_param('MS1_Tor',MS1_Tor)
            temp_MRM.set_param('IS_RT_Tor',IS_RT_Tor)
            temp_MRM.set_param('Target_RT_Tor',RT_Tor)
            temp_MRM.set_param('min_Int',min_Int)
            temp_MRM.set_param('Points',Points)
            temp_MRM.load_TargetList(TargetPath=List_Path,DefFilePath=Output_Path)
            temp_MRM.detect_Peak_MRM(FilePath=Output_Path)
            return copy.deepcopy(temp_MRM.TargetList),copy.deepcopy(temp_MRM.Calibrant)
        else:
            temp_MRM = MRMProcess(FilePath)
            temp_MRM.set_param('MS1_Tor',MS1_Tor)
            temp_MRM.set_param('IS_RT_Tor',IS_RT_Tor)
            temp_MRM.set_param('Target_RT_Tor',RT_Tor)
            temp_MRM.set_param('min_Int',min_Int)
            temp_MRM.set_param('Points',Points)
            temp_MRM.load_TargetList(TargetPath=List_Path,DefFilePath=Output_Path)
            temp_MRM.detect_Blank(FilePath=Output_Path)
            return copy.deepcopy(temp_MRM.TargetList),copy.deepcopy(temp_MRM.Calibrant)
    else:
        temp_MRM = MRMProcess(FilePath)
        temp_MRM.set_param('MS1_Tor',MS1_Tor)
        temp_MRM.set_param('IS_RT_Tor',IS_RT_Tor)
        temp_MRM.set_param('Target_RT_Tor',RT_Tor)
        temp_MRM.set_param('min_Int',min_Int)
        temp_MRM.set_param('Points',Points)
        name_begin = FilePath.rfind('/')
        if name_begin == -1:
            name_begin = FilePath.rfind('\\')
        Name = FilePath[name_begin+1:len(FilePath)-5]
        temp_MRM.load_ManualList(List_Path,Name)
        return copy.deepcopy(temp_MRM.TargetList),copy.deepcopy(temp_MRM.Calibrant)

class MRMProcess(object):   
    def __init__(self,DataPath):
        self.OriginData = pyopenms.MSExperiment()
        self.file_path = DataPath
        ''' 储存文件pyopenms.MzMLFile().store("filtered.mzML", exp) '''
        pyopenms.MzMLFile().load(self.file_path,self.OriginData)
        self.OriginData.sortSpectra(True)
        self.__param = {'MS1_Tor':0.35,'RT_Tor':6,'IS_RT_Tor':30,'Target_RT_Tor':6,'min_Int':2000,'min_IS_Int':10000,'min_MZ':100,'min_RT':60,'min_RT_width':6,'max_Noise':2000,
                        'RI_Tor':0.02,'Deconvolution':False,'FeatureDetectPlot':2,'MergeRule':'Union',
                        'UpDown_gap':10,'saveAutoList':False,'smooth':5,'Points':8,'DeconvolutionSimilarityScore':0.98,
                        'assign MS2':True,'FlowRT':2}
        temp_DF_List = []
        for i in self.OriginData.getChromatograms():
            RT,Int = i.get_peaks()
            Polarity = 'Unknow'
            try:
                if i.getNativeID().startswith('-'):
                    Polarity = 'Negative'
                else:
                    Polarity = 'Positive'
            except Exception as e:
                print("No polarity",e)
            temp_DF = pd.DataFrame({'Polarity':Polarity,'Pre_MZ':i.getPrecursor().getMZ(),'Pro_MZ':i.getProduct().getMZ(), 'Int':[Int], 'RTList':[RT]})
            temp_DF_List.append(temp_DF)
        self.Auto_Peak = pd.concat(temp_DF_List,ignore_index=True)
        self.__param['Flow_RT'] = int(self.__param['min_RT_width']/(1+4*sum(np.diff(self.Auto_Peak.loc[0,'RTList']))/(len(np.diff(self.Auto_Peak.loc[0,'RTList']))+1)))
        
    def load_TargetList(self,TargetPath,DefFilePath=''):
        def generate_and_merge_intervals(elements, x):
            counts = [1]*len(elements)
            # 1. 生成初始区间
            intervals = []
            for element, count in zip(elements, counts):
                center = element
                intervals.append((center - x, center + x, count))
            
            # 2. 按区间起始点排序
            intervals.sort(key=lambda interval: interval[0])
            
            # 3. 合并重叠区间
            merged_intervals = []
            element_counts = []
            
            for interval in intervals:
                start, end, count = interval
                
                # 如果是第一个区间，直接添加
                if not merged_intervals:
                    merged_intervals.append((start, end))
                    element_counts.append(count)
                else:
                    # 检查当前区间是否与最后一个合并区间重叠
                    last_start, last_end = merged_intervals[-1]
                    
                    # 判断是否重叠：当前区间的起始点 <= 上一个区间的结束点
                    if start <= last_end:
                        # 合并区间，取更大的结束点
                        merged_intervals[-1] = (last_start, max(last_end, end))
                        # 合并计数
                        element_counts[-1] += count
                    else:
                        # 不重叠，添加新区间
                        merged_intervals.append((start, end))
                        element_counts.append(count)
            
            return merged_intervals, element_counts
        self.TargetList = pd.read_excel(TargetPath, engine='openpyxl')
        if 'IS' in self.TargetList.keys():
            self.Calibrant = self.TargetList[self.TargetList['IS'].isna()]
            self.Calibrant.reset_index(drop=True,inplace=True)
            MRM_List = [(x,y,z) for x,y,z in zip(self.Calibrant['1_Q1'],self.Calibrant['1_Q3'],self.Calibrant['Polarity'])]
            MRM_List = list(set(MRM_List))
            CalibrantList = []
            for i in MRM_List:
                temp_Calibrant = self.Calibrant[(abs(self.Calibrant['1_Q1']-i[0])<=self.get_param('MS1_Tor'))&(abs(self.Calibrant['1_Q3']-i[1])<=self.get_param('MS1_Tor'))&(self.Calibrant['Polarity']==i[2])]
                temp_CL_set = pd.DataFrame({'Identity keys':[list(temp_Calibrant.loc[temp_Calibrant.index,'Identity keys'])],'RT':[list(temp_Calibrant.loc[temp_Calibrant.index,'RT'])],'Polarity':i[2],'1_Q1':i[0],'1_Q3':i[1],'2_Q1':temp_Calibrant.loc[temp_Calibrant.index[0],'2_Q1'],'2_Q3':temp_Calibrant.loc[temp_Calibrant.index[0],'2_Q3'],'TL_Index':[list(temp_Calibrant.index)]})
                CalibrantList.append(temp_CL_set)
            self.Calibrant_set=pd.concat(CalibrantList,ignore_index=True)
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
            bar = Bar('Pick Calibrants', max=len(self.Calibrant_set))
            for i in range(len(self.Calibrant_set)):  
                bar.next()
                MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.Calibrant_set.loc[i,'1_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.Calibrant_set.loc[i,'1_Q3'])<=self.get_param('MS1_Tor'))&(self.Auto_Peak['Polarity']==self.Calibrant_set.loc[i,'Polarity'])]
                if pd.notna(self.Calibrant_set.loc[i,'2_Q1']):
                    MRM_2 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.Calibrant_set.loc[i,'2_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.Calibrant_set.loc[i,'2_Q3'])<=self.get_param('MS1_Tor'))&(self.Auto_Peak['Polarity']==self.Calibrant_set.loc[i,'Polarity'])]
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
                    if len(DefFilePath)>0:
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
                    ''' MRM-1 '''
                    Noise_1 = self.Calculate_Noise(Origin_RT_List_1,Origin_Int_List_1)
                    Baseline_1 = MRMProcess.calc_global_baseline_with_mask_index(Origin_RT_List_1,Origin_Int_List_1,[],[],Noise_1,q=0.20,thr_k=4.0,pad_s=1.5,pad_points=None,min_nonpeak_frac=0.05,max_iter=3,drop_zeros=True)
                    time_seq,Count_seq = generate_and_merge_intervals([x*60 for x in self.Calibrant_set.loc[i,'RT']],self.get_param('IS_RT_Tor')+6)
                    All_Peak_Begin_1 = []
                    All_Peak_End_1 = []
                    All_Peak_Top_1 = []
                    All_Area_List_1 = []
                    All_Int_List_1 = []
                    All_RT_List_1 = []
                    All_RT_L_1 = []
                    All_RT_R_1 = []
                    for i_seq in range(len(time_seq)):
                        Index_forDetect = [x for x in range(len(Origin_RT_List_1)) if time_seq[i_seq][0]<=Origin_RT_List_1[x]<=time_seq[i_seq][1]]
                        RT_List_forDetect = [Origin_RT_List_1[x] for x in range(len(Origin_RT_List_1)) if time_seq[i_seq][0]<=Origin_RT_List_1[x]<=time_seq[i_seq][1]]
                        Int_List_forDetect = [Origin_Int_List_1[x] for x in range(len(Origin_RT_List_1)) if time_seq[i_seq][0]<=Origin_RT_List_1[x]<=time_seq[i_seq][1]]
                        Peak_Result_1 = self.MRMPeakDetecter(RT_List_forDetect,Int_List_forDetect,Baseline_1)
                        (Peak_Begin_1,Peak_End_1,Peak_Top_1,Area_List_1,Int_List_1) = Peak_Result_1
                        if len(Peak_Begin_1) > 0 and len(Peak_Begin_1) < Count_seq[i_seq]:
                            Ex_Peak_Begin_1 = [0]+Peak_End_1
                            Ex_Peak_End_1 = Peak_Begin_1+[len(Int_List_forDetect)]
                            for v in range(len(Ex_Peak_Begin_1)):
                                if Ex_Peak_Begin_1[v]-Ex_Peak_End_1[v] >= self.get_param('Points') and max(Int_List_forDetect[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]) >= self.get_param('min_Int'):
                                    temp_RT = RT_List_forDetect[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]
                                    temp_Int = Int_List_forDetect[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]
                                    Peak_Result_1 = self.MRMPeakDetecter(temp_RT,temp_Int)
                                    (Ex_Peak_Begin_1,Ex_Peak_End_1,Ex_Peak_Top_1,Ex_Area_List_1,Ex_Int_List_1) = Peak_Result_1
                                    Ex_Time_Begin_1 = [temp_RT[x] for x in Ex_Peak_Begin_1]
                                    Ex_Time_End_1 = [temp_RT[x] for x in Ex_Peak_End_1]
                                    Ex_Time_Top_1 = [temp_RT[x] for x in Ex_Peak_Top_1]
                                    Ex_Peak_Begin_1 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Begin_1]
                                    Ex_Peak_End_1 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_End_1]
                                    Ex_Peak_Top_1 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Top_1]
                                    Peak_Begin_1 += Ex_Peak_Begin_1
                                    Peak_End_1 += Ex_Peak_End_1
                                    Peak_Top_1 += Ex_Peak_Top_1
                                    Area_List_1 += Ex_Area_List_1
                                    Int_List_1 += Ex_Int_List_1
                        Index_1 = [x for x in range(len(Int_List_1)) if Int_List_1[x]>=self.get_param('min_Int')]
                        Peak_Begin_1 = [Peak_Begin_1[x] for x in Index_1]
                        Peak_Begin_1 = [Index_forDetect[x] for x in Peak_Begin_1]
                        Peak_End_1 = [Peak_End_1[x] for x in Index_1]
                        Peak_End_1 = [Index_forDetect[x] for x in Peak_End_1]
                        Peak_Top_1 = [Peak_Top_1[x] for x in Index_1]
                        Peak_Top_1 = [Index_forDetect[x] for x in Peak_Top_1]
                        Area_List_1 = [Area_List_1[x] for x in Index_1]
                        Int_List_1 = [Int_List_1[x] for x in Index_1]
                        RT_List_1 = [Origin_RT_List_1[x] for x in Peak_Top_1]
                        RT_L_1 = [Origin_RT_List_1[x] for x in Peak_Begin_1]
                        RT_R_1 = [Origin_RT_List_1[x] for x in Peak_End_1]
                        All_Peak_Begin_1 = All_Peak_Begin_1 + Peak_Begin_1
                        All_Peak_End_1 = All_Peak_End_1 + Peak_End_1
                        All_Peak_Top_1 = All_Peak_Top_1 + Peak_Top_1
                        All_Area_List_1 = All_Area_List_1 + Area_List_1
                        All_Int_List_1 = All_Int_List_1 + Int_List_1
                        All_RT_List_1 = All_RT_List_1 + RT_List_1
                        All_RT_L_1 = All_RT_L_1 + RT_L_1
                        All_RT_R_1 = All_RT_R_1 + RT_R_1
                    Peak_Begin_1 = All_Peak_Begin_1
                    Peak_End_1 = All_Peak_End_1
                    Peak_Top_1 = All_Peak_Top_1
                    Area_List_1 = All_Area_List_1
                    Int_List_1 = All_Int_List_1
                    RT_List_1 = All_RT_List_1
                    RT_L_1 = All_RT_L_1
                    RT_R_1 = All_RT_R_1
                    ''' MRM-2 '''
                    Noise_2 = self.Calculate_Noise(Origin_RT_List_2,Origin_Int_List_2)
                    Baseline_2 = MRMProcess.calc_global_baseline_with_mask_index(Origin_RT_List_2,Origin_Int_List_2,[],[],Noise_2,q=0.20,thr_k=4.0,pad_s=1.5,pad_points=None,min_nonpeak_frac=0.05,max_iter=3,drop_zeros=True)
                    All_Peak_Begin_2 = []
                    All_Peak_End_2 = []
                    All_Peak_Top_2 = []
                    All_Area_List_2 = []
                    All_Int_List_2 = []
                    All_RT_List_2 = []
                    All_RT_L_2 = []
                    All_RT_R_2 = []
                    for i_seq in range(len(time_seq)):
                        Index_forDetect = [x for x in range(len(Origin_RT_List_2)) if time_seq[i_seq][0]<=Origin_RT_List_2[x]<=time_seq[i_seq][1]]
                        RT_List_forDetect = [Origin_RT_List_2[x] for x in range(len(Origin_RT_List_2)) if time_seq[i_seq][0]<=Origin_RT_List_2[x]<=time_seq[i_seq][1]]
                        Int_List_forDetect = [Origin_Int_List_2[x] for x in range(len(Origin_RT_List_2)) if time_seq[i_seq][0]<=Origin_RT_List_2[x]<=time_seq[i_seq][1]]
                        Peak_Result_2 = self.MRMPeakDetecter(RT_List_forDetect,Int_List_forDetect,Baseline_2)
                        (Peak_Begin_2,Peak_End_2,Peak_Top_2,Area_List_2,Int_List_2) = Peak_Result_2
                        if len(Peak_Begin_2) > 0 and len(Peak_Begin_2) < Count_seq[i_seq]:
                            Ex_Peak_Begin_2 = [0]+Peak_End_2
                            Ex_Peak_End_2 = Peak_Begin_2+[len(Int_List_forDetect)]
                            for v in range(len(Ex_Peak_Begin_2)):
                                if Ex_Peak_Begin_2[v]-Ex_Peak_End_2[v] >= self.get_param('Points') and max(Int_List_forDetect[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]) >= self.get_param('min_Int'):
                                    temp_RT = RT_List_forDetect[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]
                                    temp_Int = Int_List_forDetect[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]
                                    Peak_Result_2 = self.MRMPeakDetecter(temp_RT,temp_Int)
                                    (Ex_Peak_Begin_2,Ex_Peak_End_2,Ex_Peak_Top_2,Ex_Area_List_2,Ex_Int_List_2) = Peak_Result_2
                                    Ex_Time_Begin_2 = [temp_RT[x] for x in Ex_Peak_Begin_2]
                                    Ex_Time_End_2 = [temp_RT[x] for x in Ex_Peak_End_2]
                                    Ex_Time_Top_2 = [temp_RT[x] for x in Ex_Peak_Top_2]
                                    Ex_Peak_Begin_2 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Begin_2]
                                    Ex_Peak_End_2 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_End_2]
                                    Ex_Peak_Top_2 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Top_2]
                                    Peak_Begin_2 += Ex_Peak_Begin_2
                                    Peak_End_2 += Ex_Peak_End_2
                                    Peak_Top_2 += Ex_Peak_Top_2
                                    Area_List_2 += Ex_Area_List_2
                                    Int_List_2 += Ex_Int_List_2
                        Index_2 = [x for x in range(len(Int_List_2)) if Int_List_2[x]>=self.get_param('min_Int')]
                        Peak_Begin_2 = [Peak_Begin_2[x] for x in Index_2]
                        Peak_Begin_2 = [Index_forDetect[x] for x in Peak_Begin_2]
                        Peak_End_2 = [Peak_End_2[x] for x in Index_2]
                        Peak_End_2 = [Index_forDetect[x] for x in Peak_End_2]
                        Peak_Top_2 = [Peak_Top_2[x] for x in Index_2]
                        Peak_Top_2 = [Index_forDetect[x] for x in Peak_Top_2]
                        Area_List_2 = [Area_List_2[x] for x in Index_2]
                        Int_List_2 = [Int_List_2[x] for x in Index_2]
                        RT_List_2 = [Origin_RT_List_2[x] for x in Peak_Top_2]
                        RT_L_2 = [Origin_RT_List_2[x] for x in Peak_Begin_2]
                        RT_R_2 = [Origin_RT_List_2[x] for x in Peak_End_2]
                        All_Peak_Begin_2 = All_Peak_Begin_2 + Peak_Begin_2
                        All_Peak_End_2 = All_Peak_End_2 + Peak_End_2
                        All_Peak_Top_2 = All_Peak_Top_2 + Peak_Top_2
                        All_Area_List_2 = All_Area_List_2 + Area_List_2
                        All_Int_List_2 = All_Int_List_2 + Int_List_2
                        All_RT_List_2 = All_RT_List_2 + RT_List_2
                        All_RT_L_2 = All_RT_L_2 + RT_L_2
                        All_RT_R_2 = All_RT_R_2 + RT_R_2
                    Peak_Begin_2 = All_Peak_Begin_2
                    Peak_End_2 = All_Peak_End_2
                    Peak_Top_2 = All_Peak_Top_2
                    Area_List_2 = All_Area_List_2
                    Int_List_2 = All_Int_List_2
                    RT_List_2 = All_RT_List_2
                    RT_L_2 = All_RT_L_2
                    RT_R_2 = All_RT_R_2
                    Score_matrix= np.zeros([len(RT_List_1),len(RT_List_2)])
                    for i_1 in range(len(RT_List_1)):
                        for i_2 in range(len(RT_List_2)):
                            if abs(RT_List_1[i_1]-RT_List_2[i_2])<self.get_param('Target_RT_Tor'):
                                Score_matrix[i_1,i_2] = 1-abs(RT_List_1[i_1]-RT_List_2[i_2])/self.get_param('Target_RT_Tor')
                    Sm_1,Sm_2 = linear_sum_assignment(Score_matrix,True)
                    Sm_Filter = [(x,y) for x,y in zip(Sm_1,Sm_2) if Score_matrix[x,y]>0]
                    
                    Score_matrix= np.zeros([len(Sm_Filter),len(self.Calibrant_set.loc[i,'RT'])])
                    for i_1 in range(len(Sm_Filter)):
                        for i_2 in range(len(self.Calibrant_set.loc[i,'RT'])):
                            if abs(RT_List_1[Sm_Filter[i_1][0]]-self.Calibrant_set.loc[i,'RT'][i_2]*60)< self.get_param('IS_RT_Tor') and abs(RT_List_2[Sm_Filter[i_1][1]]-self.Calibrant_set.loc[i,'RT'][i_2]*60)< self.get_param('IS_RT_Tor'):
                                Score_matrix[i_1,i_2] = 0.25*(2-abs(RT_List_1[Sm_Filter[i_1][0]]-self.Calibrant_set.loc[i,'RT'][i_2]*60)/ self.get_param('IS_RT_Tor') - abs(RT_List_2[Sm_Filter[i_1][1]]-self.Calibrant_set.loc[i,'RT'][i_2]*60)/ self.get_param('IS_RT_Tor')) + 0.25*(Int_List_1[Sm_Filter[i_1][0]]/max(Int_List_1) + Int_List_2[Sm_Filter[i_1][1]]/max(Int_List_2))
                    Sm_1,Sm_2 = linear_sum_assignment(Score_matrix,True)         
                    Sm_Assign = [(x,y) for x,y in zip(Sm_1,Sm_2) if Score_matrix[x,y]>0]
                    for i_1,i_2 in Sm_Assign:
                        self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'RT_1'] = np.around(RT_List_1[Sm_Filter[i_1][0]]/60,3)
                        self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Area_1'] = int(Area_List_1[Sm_Filter[i_1][0]])
                        self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Hight_1'] = int(Int_List_1[Sm_Filter[i_1][0]])
                        self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Noise_1'] = int(Noise_1)
                        self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Edge_left_1'] = np.around(RT_L_1[Sm_Filter[i_1][0]]/60,4)
                        self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Edge_right_1'] = np.around(RT_R_1[Sm_Filter[i_1][0]]/60,4)
                        self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'RT_Calibrant'] = np.around(RT_List_1[Sm_Filter[i_1][0]]/60,2)
                        if pd.notna(self.Calibrant_set.loc[i,'2_Q1']):
                            self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'RT_2'] = np.around(RT_List_2[Sm_Filter[i_1][1]]/60,3)
                            self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Area_2'] = int(Area_List_2[Sm_Filter[i_1][1]])
                            self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Hight_2'] = int(Int_List_2[Sm_Filter[i_1][1]])
                            self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Noise_2'] = int(Noise_2)
                            self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Edge_left_2'] = np.around(RT_L_2[Sm_Filter[i_1][1]]/60,4)
                            self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Edge_right_2'] = np.around(RT_R_2[Sm_Filter[i_1][1]]/60,4)
                    if len(DefFilePath)>0:
                        for i_1,i_2 in Sm_Assign:
                            fig_1.add_trace(
                                go.Scatter(
                                    x=[x/60 for x in Origin_RT_List_1][Peak_Begin_1[Sm_Filter[i_1][0]]:Peak_End_1[Sm_Filter[i_1][0]]+1],
                                    y=Origin_Int_List_1[Peak_Begin_1[Sm_Filter[i_1][0]]:Peak_End_1[Sm_Filter[i_1][0]]+1],
                                    mode='lines',
                                    name=self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Identity keys'],
                                    line={'width':1},
                                    fill='tozeroy',
                                    ))
                        fig_1.update_layout(
                            width=600,
                            height=500,
                            titlefont={'size':10},
                            title='Q1:'+str(np.around(self.Calibrant_set.loc[i,'1_Q1'],1))+' Q3:'+str(np.around(self.Calibrant_set.loc[i,'1_Q3'],1)),
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
                        fig_1.write_html(DefFilePath+'/'+'IS-Q1-'+str(np.around(self.Calibrant_set.loc[i,'1_Q1'],1))+' Q3-'+str(np.around(self.Calibrant_set.loc[i,'1_Q3'],1))+'.html',config={'responsive': False})
                        if pd.notna(self.Calibrant.loc[i,'2_Q1']):
                            for i_1,i_2 in Sm_Assign:
                                fig_2.add_trace(
                                    go.Scatter(
                                        x=[x/60 for x in Origin_RT_List_2][Peak_Begin_2[Sm_Filter[i_1][1]]:Peak_End_2[Sm_Filter[i_1][1]]+1],
                                        y=Origin_Int_List_2[Peak_Begin_2[Sm_Filter[i_1][1]]:Peak_End_2[Sm_Filter[i_1][1]]+1],
                                        mode='lines',
                                        name=self.Calibrant.loc[self.Calibrant_set.loc[i,'TL_Index'][i_2],'Identity keys'],
                                        line={'width':1},
                                        fill='tozeroy',
                                        ))
                            fig_2.update_layout(
                                width=600,
                                height=500,
                                titlefont={'size':10},
                                title='Q1:'+str(np.around(self.Calibrant_set.loc[i,'2_Q1'],1))+' Q3:'+str(np.around(self.Calibrant_set.loc[i,'2_Q3'],1)),
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
                            fig_2.write_html(DefFilePath+'/'+'IS-Q1-'+str(np.around(self.Calibrant_set.loc[i,'2_Q1'],1))+' Q3-'+str(np.around(self.Calibrant_set.loc[i,'2_Q3'],1))+'.html',config={'responsive': False})
                        self.Calibrant.to_excel(DefFilePath+'/Calibrant-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
        if 'IS' in self.TargetList.keys() and 'Edge_left_1' not in self.TargetList.keys():
            self.Calibrant = self.Calibrant[self.Calibrant['RT_Calibrant'] != '']
            self.Calibrant.reset_index(drop=True,inplace=True)
            self.TargetList.loc[:,'RT From List'] = self.TargetList.loc[:,'RT']
            if len(self.Calibrant) > 0:
                for i in range(len(self.TargetList)):  
                        Diff_RT = self.Calibrant['RT'] - self.TargetList.loc[i,'RT']
                        Diff_Left = [Diff_RT[x] for x in range(len(Diff_RT)) if Diff_RT[x]<=0]
                        Index_Left = [x for x in range(len(Diff_RT)) if Diff_RT[x]<=0]
                        Diff_Right = [Diff_RT[x] for x in range(len(Diff_RT)) if Diff_RT[x]>=0]
                        Index_Right = [x for x in range(len(Diff_RT)) if Diff_RT[x]>=0]
                        if len(Diff_Left)>0 and len(Diff_Right) == 0:
                            Index_L = np.where(Diff_Left==max(Diff_Left))[0][0]
                            RT_Fix = self.TargetList.loc[i,'RT']*self.Calibrant.loc[Index_Left[Index_L],'RT_Calibrant']/self.Calibrant.loc[Index_Left[Index_L],'RT']
                        elif len(Diff_Left)==0 and len(Diff_Right) > 0:
                            Index_R = np.where(Diff_Right==min(Diff_Right))[0][0]
                            RT_Fix = self.TargetList.loc[i,'RT']*self.Calibrant.loc[Index_Right[Index_R],'RT_Calibrant']/self.Calibrant.loc[Index_Right[Index_R],'RT']
                        else:
                            Index_L = np.where(Diff_Left==max(Diff_Left))[0][0]
                            Index_R = np.where(Diff_Right==min(Diff_Right))[0][0]
                            RT_Fix = 0.5*(self.TargetList.loc[i,'RT']*self.Calibrant.loc[Index_Left[Index_L],'RT_Calibrant']/self.Calibrant.loc[Index_Left[Index_L],'RT']+self.TargetList.loc[i,'RT']*self.Calibrant.loc[Index_Right[Index_R],'RT_Calibrant']/self.Calibrant.loc[Index_Right[Index_R],'RT'])
                        self.TargetList.loc[i,'RT'] = RT_Fix
        MRM_List = [(x,y,z) for x,y,z in zip(self.TargetList['1_Q1'],self.TargetList['1_Q3'],self.TargetList['Polarity'])]
        MRM_List = list(set(MRM_List))
        TargetList = []
        for i in MRM_List:
            temp_TargetList = self.TargetList[(abs(self.TargetList['1_Q1']-i[0])<=self.get_param('MS1_Tor'))&(abs(self.TargetList['1_Q3']-i[1])<=self.get_param('MS1_Tor'))&(self.TargetList['Polarity']==i[2])]
            temp_TL_set = pd.DataFrame({'Identity keys':[list(temp_TargetList.loc[temp_TargetList.index,'Identity keys'])],'RT':[list(temp_TargetList.loc[temp_TargetList.index,'RT'])],'Polarity':i[2],'1_Q1':i[0],'1_Q3':i[1],'2_Q1':temp_TargetList.loc[temp_TargetList.index[0],'2_Q1'],'2_Q3':temp_TargetList.loc[temp_TargetList.index[0],'2_Q3'],'TL_Index':[list(temp_TargetList.index)]})
            TargetList.append(temp_TL_set)
        self.TargetList_set=pd.concat(TargetList,ignore_index=True)
    
    def load_ManualList(self,ManualPath,Name):
        def integral(x,y):
            x_diff = list(np.diff(x))
            y_mix = [0.5*(y[x]+y[x+1]) for x in range(len(y)-1)]
            return sum([x_diff[x]*y_mix[x] for x in range(len(x_diff))])
        def generate_and_merge_intervals(elements, x):
            counts = [1]*len(elements)
            # 1. 生成初始区间
            intervals = []
            for element, count in zip(elements, counts):
                center = element
                intervals.append((center - x, center + x, count))
            
            # 2. 按区间起始点排序
            intervals.sort(key=lambda interval: interval[0])
            
            # 3. 合并重叠区间
            merged_intervals = []
            element_counts = []
            
            for interval in intervals:
                start, end, count = interval
                
                # 如果是第一个区间，直接添加
                if not merged_intervals:
                    merged_intervals.append((start, end))
                    element_counts.append(count)
                else:
                    # 检查当前区间是否与最后一个合并区间重叠
                    last_start, last_end = merged_intervals[-1]
                    
                    # 判断是否重叠：当前区间的起始点 <= 上一个区间的结束点
                    if start <= last_end:
                        # 合并区间，取更大的结束点
                        merged_intervals[-1] = (last_start, max(last_end, end))
                        # 合并计数
                        element_counts[-1] += count
                    else:
                        # 不重叠，添加新区间
                        merged_intervals.append((start, end))
                        element_counts.append(count)
        self.TargetList = pd.read_excel(ManualPath, engine='openpyxl')
        self.TargetList['RT_1'] = ''
        self.TargetList['Area_1'] = ''
        self.TargetList['Hight_1'] = ''
        self.TargetList['Noise_1'] = ''
        self.TargetList['Edge_left_1'] = ''
        self.TargetList['Edge_right_1'] = ''
        self.TargetList['Noise_1'] = ''
        self.TargetList['Plot_RT_1'] = ''
        self.TargetList['Plot_Int_1'] = ''
        self.TargetList['Plot_RT_All_1'] = ''
        self.TargetList['Plot_Int_All_1'] = ''
        self.TargetList['RT_2'] = ''
        self.TargetList['Area_2'] = ''
        self.TargetList['Hight_2'] = ''
        self.TargetList['Noise_2'] = ''
        self.TargetList['Edge_left_2'] = ''
        self.TargetList['Edge_right_2'] = ''
        self.TargetList['Noise_2'] = ''
        self.Calibrant = self.TargetList[self.TargetList['IS'].isna()]
        self.Calibrant.reset_index(drop=True,inplace=True)
        for i in range(len(self.Calibrant)):
            if self.Calibrant.loc[i,Name+'-L'] >0:
                MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.Calibrant.loc[i,'1_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.Calibrant.loc[i,'1_Q3'])<=self.get_param('MS1_Tor'))&(self.Auto_Peak['Polarity']==self.Calibrant.loc[i,'Polarity'])]
                Origin_RT_List_1 = []
                Origin_Int_List_1 = []
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
                self.Calibrant.at[i,'Plot_RT_All_1'] = Origin_RT_List_1
                self.Calibrant.at[i,'Plot_Int_All_1'] = Origin_Int_List_1
                Noise_1 = self.Calculate_Noise(Origin_RT_List_1,Origin_Int_List_1)
                Baseline_1 = MRMProcess.calc_global_baseline_with_mask_index(Origin_RT_List_1,Origin_Int_List_1,[],[],Noise_1,q=0.20,thr_k=4.0,pad_s=1.5,pad_points=None,min_nonpeak_frac=0.05,max_iter=3,drop_zeros=True)
                L_Place = np.where(abs(np.array(Origin_RT_List_1)-self.Calibrant.loc[i,Name+'-L']*60)==min(abs(np.array(Origin_RT_List_1)-self.Calibrant.loc[i,Name+'-L']*60)))
                R_Place = np.where(abs(np.array(Origin_RT_List_1)-self.Calibrant.loc[i,Name+'-R']*60)==min(abs(np.array(Origin_RT_List_1)-self.Calibrant.loc[i,Name+'-R']*60)))
                if len(L_Place[0]) == 1:
                    L_Place = L_Place[0][0]
                else:
                    L_Place_temp = L_Place[0][0]
                    for i_Place in L_Place[0]:
                        if Origin_Int_List_1[L_Place_temp] >= Origin_Int_List_1[L_Place[0][i_Place]]:
                            L_Place_temp = L_Place[0][i_Place]
                    L_Place = L_Place_temp
                if len(R_Place[0]) == 1:
                    R_Place = R_Place[0][0]
                else:
                    R_Place_temp = R_Place[0][0]
                    for i_Place in R_Place[0]:
                        if Origin_Int_List_1[R_Place_temp] > Origin_Int_List_1[R_Place[0][i_Place]]:
                            R_Place_temp = R_Place[0][i_Place]
                    R_Place = R_Place_temp
                Int_List = Origin_Int_List_1[L_Place:R_Place+1]
                RT_List = Origin_RT_List_1[L_Place:R_Place+1]
                self.Calibrant.at[i,'Plot_RT_1'] = RT_List
                self.Calibrant.at[i,'Plot_Int_1'] = Int_List
                Area = max(0,integral(RT_List, Int_List)-0.5*(min(Baseline_1,Int_List[0])+min(Baseline_1,Int_List[-1]))*(RT_List[-1]-RT_List[0]))
                self.Calibrant.loc[i,'RT_1'] = [RT_List[x] for x in range(len(Int_List)) if Int_List[x] == max(Int_List)][0]
                self.Calibrant.loc[i,'Area_1'] = Area
                self.Calibrant.loc[i,'Hight_1'] = max(Int_List)
                self.Calibrant.loc[i,'Noise_1'] = Noise_1
                self.Calibrant.loc[i,'Edge_left_1'] = RT_List[0]
                self.Calibrant.loc[i,'Edge_right_1'] = RT_List[-1]
                self.Calibrant.loc[i,'RT_2'] = [RT_List[x] for x in range(len(Int_List)) if Int_List[x] == max(Int_List)][0]
                self.Calibrant.loc[i,'Area_2'] = Area
                self.Calibrant.loc[i,'Hight_2'] = max(Int_List)
                self.Calibrant.loc[i,'Noise_2'] = Noise_1
                self.Calibrant.loc[i,'Edge_left_2'] = RT_List[0]
                self.Calibrant.loc[i,'Edge_right_2'] = RT_List[-1]
        self.TargetList = self.TargetList[self.TargetList['IS'].notna()]
        self.TargetList.reset_index(drop=True,inplace=True)
        self.TargetList['Area_1/IS'] = ''
        self.TargetList['Area_2/IS'] = ''
        for i in range(len(self.TargetList)):
            if self.TargetList.loc[i,Name+'-L'] >0:
                MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList.loc[i,'1_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList.loc[i,'1_Q3'])<=self.get_param('MS1_Tor'))&(self.Auto_Peak['Polarity']==self.TargetList.loc[i,'Polarity'])]
                Origin_RT_List_1 = []
                Origin_Int_List_1 = []
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
                self.TargetList.at[i,'Plot_RT_All_1'] = Origin_RT_List_1
                self.TargetList.at[i,'Plot_Int_All_1'] = Origin_Int_List_1
                Noise_1 = self.Calculate_Noise(Origin_RT_List_1,Origin_Int_List_1)
                Baseline_1 = MRMProcess.calc_global_baseline_with_mask_index(Origin_RT_List_1,Origin_Int_List_1,[],[],Noise_1,q=0.20,thr_k=4.0,pad_s=1.5,pad_points=None,min_nonpeak_frac=0.05,max_iter=3,drop_zeros=True)
                L_Place = np.where(abs(np.array(Origin_RT_List_1)-self.TargetList.loc[i,Name+'-L']*60)==min(abs(np.array(Origin_RT_List_1)-self.TargetList.loc[i,Name+'-L']*60)))
                R_Place = np.where(abs(np.array(Origin_RT_List_1)-self.TargetList.loc[i,Name+'-R']*60)==min(abs(np.array(Origin_RT_List_1)-self.TargetList.loc[i,Name+'-R']*60)))
                if len(L_Place[0]) == 1:
                    L_Place = L_Place[0][0]
                else:
                    L_Place_temp = L_Place[0][0]
                    for i_Place in L_Place[0]:
                        if Origin_Int_List_1[L_Place_temp] >= Origin_Int_List_1[L_Place[0][i_Place]]:
                            L_Place_temp = L_Place[0][i_Place]
                    L_Place = L_Place_temp
                if len(R_Place[0]) == 1:
                    R_Place = R_Place[0][0]
                else:
                    R_Place_temp = R_Place[0][0]
                    for i_Place in R_Place[0]:
                        if Origin_Int_List_1[R_Place_temp] > Origin_Int_List_1[R_Place[0][i_Place]]:
                            R_Place_temp = R_Place[0][i_Place]
                    R_Place = R_Place_temp
                Int_List = Origin_Int_List_1[L_Place:R_Place+1]
                RT_List = Origin_RT_List_1[L_Place:R_Place+1]
                Area = max(0,integral(RT_List, Int_List)-0.5*(min(Baseline_1,Int_List[0])+min(Baseline_1,Int_List[-1]))*(RT_List[-1]-RT_List[0]))
                self.TargetList.loc[i,'RT_1'] = [RT_List[x] for x in range(len(Int_List)) if Int_List[x] == max(Int_List)][0]
                self.TargetList.loc[i,'Area_1'] = Area
                IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[i,'IS']]['Area_1'])[0]
                self.TargetList.loc[i,'Area_1/IS'] = Area/IS_Area
                self.TargetList.loc[i,'Hight_1'] = max(Int_List)
                self.TargetList.loc[i,'Noise_1'] = Noise_1
                self.TargetList.at[i,'Plot_RT_1'] = RT_List
                self.TargetList.at[i,'Plot_Int_1'] = Int_List
                self.TargetList.loc[i,'Edge_left_1'] = RT_List[0]
                self.TargetList.loc[i,'Edge_right_1'] = RT_List[-1]
                self.TargetList.loc[i,'RT_2'] = [RT_List[x] for x in range(len(Int_List)) if Int_List[x] == max(Int_List)][0]
                self.TargetList.loc[i,'Area_2'] = Area
                self.TargetList.loc[i,'Area_2/IS'] = Area/IS_Area
                self.TargetList.loc[i,'Hight_2'] = max(Int_List)
                self.TargetList.loc[i,'Noise_2'] = Noise_1
                self.TargetList.at[i,'Plot_RT_1'] = RT_List
                self.TargetList.at[i,'Plot_Int_1'] = Int_List
                self.TargetList.loc[i,'Edge_left_2'] = RT_List[0]
                self.TargetList.loc[i,'Edge_right_2'] = RT_List[-1]
            
    def load_TargetList_only(self,TargetPath,DefFilePath=''):
        def generate_and_merge_intervals(elements, x):
            counts = [1]*len(elements)
            # 1. 生成初始区间
            intervals = []
            for element, count in zip(elements, counts):
                center = element
                intervals.append((center - x, center + x, count))
            
            # 2. 按区间起始点排序
            intervals.sort(key=lambda interval: interval[0])
            
            # 3. 合并重叠区间
            merged_intervals = []
            element_counts = []
            
            for interval in intervals:
                start, end, count = interval
                
                # 如果是第一个区间，直接添加
                if not merged_intervals:
                    merged_intervals.append((start, end))
                    element_counts.append(count)
                else:
                    # 检查当前区间是否与最后一个合并区间重叠
                    last_start, last_end = merged_intervals[-1]
                    
                    # 判断是否重叠：当前区间的起始点 <= 上一个区间的结束点
                    if start <= last_end:
                        # 合并区间，取更大的结束点
                        merged_intervals[-1] = (last_start, max(last_end, end))
                        # 合并计数
                        element_counts[-1] += count
                    else:
                        # 不重叠，添加新区间
                        merged_intervals.append((start, end))
                        element_counts.append(count)
            
            return merged_intervals, element_counts
        self.TargetList = pd.read_excel(TargetPath, engine='openpyxl')
        MRM_List = [(x,y,z) for x,y,z in zip(self.TargetList['1_Q1'],self.TargetList['1_Q3'],self.TargetList['Polarity'])]
        MRM_List = list(set(MRM_List))
        TargetList = []
        for i in MRM_List:
            temp_TargetList = self.TargetList[(abs(self.TargetList['1_Q1']-i[0])<=self.get_param('MS1_Tor'))&(abs(self.TargetList['1_Q3']-i[1])<=self.get_param('MS1_Tor'))&(self.TargetList['Polarity']==i[2])]
            temp_TL_set = pd.DataFrame({'Identity keys':[list(temp_TargetList.loc[temp_TargetList.index,'Identity keys'])],'RT':[list(temp_TargetList.loc[temp_TargetList.index,'RT'])],'Polarity':i[2],'1_Q1':i[0],'1_Q3':i[1],'2_Q1':temp_TargetList.loc[temp_TargetList.index[0],'2_Q1'],'2_Q3':temp_TargetList.loc[temp_TargetList.index[0],'2_Q3'],'TL_Index':[list(temp_TargetList.index)]})
            TargetList.append(temp_TL_set)
        self.TargetList_set=pd.concat(TargetList,ignore_index=True)
        
    def MRMPeakDetecter(self,Auto_RT_List,Auto_Int_List,BaseLine=0):
        # 单个色谱峰识别模块，由下方另一函数调用
        def GaussSmooth(x):
            if len(x)==5:
                op = x[0]*0.07+x[1]*0.23+x[2]*0.4+x[3]*0.23+x[4]*0.07
            elif len(x)==3:
                op = x[0]*0.17 +x[1]*0.66 +x[2]*0.17
            else:
                op = sum(x)/len(x)
            return op
        def integral(x,y):
            x_diff = list(np.diff(x))
            y_mix = [0.5*(y[x]+y[x+1]) for x in range(len(y)-1)]
            return sum([x_diff[x]*y_mix[x] for x in range(len(x_diff))])
        raw_Int_List = Auto_Int_List
        Auto_Int_List = np.array(list(map(lambda x:GaussSmooth(Auto_Int_List[x-1:x+1+1]) if x in range(1,len(Auto_Int_List)-1) else Auto_Int_List[x],range(len(Auto_Int_List)))))
        Diff_List = [Auto_Int_List[0]]+list(np.diff(Auto_Int_List))
        Diff_List = np.array(list(map(lambda x:GaussSmooth(Diff_List[x-1:x+1+1]) if x in range(1,len(Diff_List)-1) else Diff_List[x],range(len(Diff_List)))))
        FD_List = np.array(list(map(lambda x: MRMProcess.TFFD(x, Auto_Int_List, Auto_RT_List), range(len(Auto_Int_List)))))
        FD_List = np.array(list(map(lambda x:GaussSmooth(FD_List[x-1:x+1+1]) if x in range(1,len(FD_List)-1) else FD_List[x],range(len(FD_List)))))
        SD_List = np.array(list(map(lambda x: MRMProcess.TFSD(x, Auto_Int_List, Auto_RT_List), range(len(Auto_Int_List)))))
        SD_List = np.array(list(map(lambda x:GaussSmooth(SD_List[x-1:x+1+1]) if x in range(1,len(SD_List)-1) else SD_List[x],range(len(SD_List)))))
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
            End_Place,Complet_EP = MRMProcess.Find_FContinuous(FD_List, FD_N, Diff_N, FDMode='N',mergeRule=self.get_param('MergeRule'),FC_Number=self.get_param('FeatureDetectPlot'))
            Peak_Place = MRMProcess.Find_SDChange(FD_Change, SD_Place)  # 峰顶位置
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
                    Begin_Filter = np.where(raw_Int_List[Peak_Begin[i_edge]-self.get_param('FeatureDetectPlot'):min(len(raw_Int_List),Peak_Begin[i_edge]+2*self.get_param('FeatureDetectPlot')+1)]==min(raw_Int_List[Peak_Begin[i_edge]-self.get_param('FeatureDetectPlot'):min(len(raw_Int_List),Peak_Begin[i_edge]+2*self.get_param('FeatureDetectPlot')+1)]))[0][-1]
                    if Peak_Begin[i_edge]-self.get_param('FeatureDetectPlot')+Begin_Filter < Peak_Top[i_edge]:
                        Peak_Begin[i_edge] = Peak_Begin[i_edge]-self.get_param('FeatureDetectPlot')+Begin_Filter
                if Peak_End[i_edge]<len(raw_Int_List)-self.get_param('FeatureDetectPlot'):
                    End_Filter = np.where(raw_Int_List[max(0,Peak_End[i_edge]-2*self.get_param('FeatureDetectPlot')):Peak_End[i_edge]+self.get_param('FeatureDetectPlot')+1]==min(raw_Int_List[max(0,Peak_End[i_edge]-2*self.get_param('FeatureDetectPlot')):Peak_End[i_edge]+self.get_param('FeatureDetectPlot')+1]))[0][0]
                    if Peak_End[i_edge]+End_Filter-2*self.get_param('FeatureDetectPlot')>Peak_Top[i_edge]:
                        Peak_End[i_edge] = Peak_End[i_edge]+End_Filter-2*self.get_param('FeatureDetectPlot')
            Area_List = []
            Int_List = []
            Final = []
            for i_edge in range(len(Peak_Begin)):
                temp_Area = integral(Auto_RT_List[Peak_Begin[i_edge]:Peak_End[i_edge]+1],raw_Int_List[Peak_Begin[i_edge]:Peak_End[i_edge]+1])
                temp_Area = temp_Area - 0.5*(min(raw_Int_List[Peak_Begin[i_edge]],BaseLine)+min(raw_Int_List[Peak_End[i_edge]],BaseLine))*(Auto_RT_List[Peak_End[i_edge]]-Auto_RT_List[Peak_Begin[i_edge]])
                Area_List.append(max(0,temp_Area))
                Int_List.append(max(np.array(raw_Int_List)[Peak_Begin[i_edge]:Peak_End[i_edge]+1]))
                if max(np.array(raw_Int_List)[Peak_Begin[i_edge]:Peak_End[i_edge]+1]) >= self.get_param('min_Int'):
                    Final.append(i_edge)
            Peak_Begin = [Peak_Begin[x] for x in Final]
            Peak_End = [Peak_End[x] for x in Final]
            Peak_Top = [Peak_Top[x] for x in Final]
            Area_List = [Area_List[x] for x in Final]
            Int_List = [Int_List[x] for x in Final]
            return Peak_Begin,Peak_End,Peak_Top,Area_List,Int_List
        else:
            return [],[],[],[],[]
                
    def quick_integral(self,Q1,Q3,Polarity,T1,T2):
        def integral(x,y):
            x_diff = list(np.diff(x))
            y_mix = [0.5*(y[x]+y[x+1]) for x in range(len(y)-1)]
            return sum([x_diff[x]*y_mix[x] for x in range(len(x_diff))])
        MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-Q1)<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-Q3)<=self.get_param('MS1_Tor'))&(self.Auto_Peak['Polarity']==Polarity)]
        Origin_RT_List_1 = []
        Origin_Int_List_1 = []
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
        L_Place = np.where(abs(np.array(Origin_RT_List_1)-T1*60)==min(abs(np.array(Origin_RT_List_1)-T1*60)))
        R_Place = np.where(abs(np.array(Origin_RT_List_1)-T2*60)==min(abs(np.array(Origin_RT_List_1)-T2*60)))
        if len(L_Place[0]) == 1:
            L_Place = L_Place[0][0]
        else:
            L_Place_temp = L_Place[0][0]
            for i_Place in L_Place[0]:
                if Origin_Int_List_1[L_Place_temp] >= Origin_Int_List_1[L_Place[0][i_Place]]:
                    L_Place_temp = L_Place[0][i_Place]
            L_Place = L_Place_temp
        if len(R_Place[0]) == 1:
            R_Place = R_Place[0][0]
        else:
            R_Place_temp = R_Place[0][0]
            for i_Place in R_Place[0]:
                if Origin_Int_List_1[R_Place_temp] > Origin_Int_List_1[R_Place[0][i_Place]]:
                    R_Place_temp = R_Place[0][i_Place]
            R_Place = R_Place_temp
        Int_List = Origin_Int_List_1[L_Place:R_Place+1]
        RT_List = Origin_RT_List_1[L_Place:R_Place+1]
        return integral(RT_List, Int_List)
    
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
        noise = np.median(abs(Int_List[Int_List>0]-Auto_Int_List[Int_List>0]))
        return noise
    
    def calc_global_baseline_with_mask_index(
        RT_List,
        Int_List,
        Peak_Begin,
        Peak_End,
        Noise,
        q=0.20,
        thr_k=4.0,
        pad_s=1.5,            
        pad_points=None,      
        min_nonpeak_frac=0.05,
        max_iter=3,
        drop_zeros=True,
    ):
        def _merge_intervals_idx(intervals):
            """合并可能重叠/相邻的索引区间 intervals=[(i0,i1),...]，输出按起点排序后的合并结果。"""
            if not intervals:
                return []
            norm = []
            for a, b in intervals:
                a = int(a); b = int(b)
                if a > b:
                    a, b = b, a
                norm.append((a, b))
            norm.sort(key=lambda t: t[0])

            merged = [norm[0]]
            for a, b in norm[1:]:
                la, lb = merged[-1]
                if a <= lb + 1:
                    merged[-1] = (la, max(lb, b))
                else:
                    merged.append((a, b))
            return merged

        def _interval_mask_idx(n, intervals, pad_pts):
            mask = np.zeros(n, dtype=bool)
            for a, b in intervals:
                a2 = max(0, a - pad_pts)
                b2 = min(n - 1, b + pad_pts)
                mask[a2:b2 + 1] = True
            return mask

        def _contiguous_true_regions(mask):
            mask = np.asarray(mask, dtype=bool)
            if mask.size == 0:
                return []
            diff = np.diff(mask.astype(int))
            starts = list(np.where(diff == 1)[0] + 1)
            ends = list(np.where(diff == -1)[0])
            if mask[0]:
                starts = [0] + starts
            if mask[-1]:
                ends = ends + [mask.size - 1]
            return list(zip(starts, ends))

        def _quantile_drop0(arr, qv):
            arr = np.asarray(arr, dtype=float)
            if arr.size == 0:
                return None
            if drop_zeros:
                nz = arr[arr != 0]
                if nz.size > 0:
                    return float(np.quantile(nz, qv))
                return float(np.quantile(arr, qv))
            else:
                return float(np.quantile(arr, qv))

        # ========== check import ==========
        x = np.asarray(RT_List, dtype=float)
        y = np.asarray(Int_List, dtype=float)
        if x.shape != y.shape:
            raise ValueError("len of RT_List and Int_List")
        n = y.size
        if n < 5:
            raise ValueError("poor points")
        if Noise <= 0:
            raise ValueError("Noise error")

        begins = np.asarray(Peak_Begin, dtype=int).reshape(-1)
        ends = np.asarray(Peak_End, dtype=int).reshape(-1)
        if begins.size != ends.size:
            raise ValueError("len of Peak_Begin and Peak_End")

        # ========== pad_pts ==========
        if pad_points is not None:
            pad_pts = int(max(0, pad_points))
        else:
            dx = np.diff(x)
            dx = dx[dx > 0]
            if dx.size == 0:
                pad_pts = 1
            else:
                dt = float(np.median(dx))  
                pad_pts = int(np.round(pad_s / dt))
                pad_pts = max(0, pad_pts)

        # ========== mark ==========
        intervals = [(int(begins[i]), int(ends[i])) for i in range(begins.size)]
        clipped = []
        for a, b in intervals:
            a = max(0, min(n - 1, a))
            b = max(0, min(n - 1, b))
            clipped.append((a, b))
        intervals = _merge_intervals_idx(clipped)

        mask_peak = _interval_mask_idx(n, intervals, pad_pts)

        logs = []
        base0 = _quantile_drop0(y, q)
        if base0 is None:
            raise ValueError("empty")
        thr = base0 + thr_k * Noise

        # ========== Iterate ==========
        for it in range(max_iter):
            nonpeak = y[~mask_peak]
            if nonpeak.size < max(10, int(min_nonpeak_frac * n)):
                base0 = _quantile_drop0(y, q)
                thr = base0 + thr_k * Noise
                logs.append(f"Ture {it+1}：poor points（{nonpeak.size}），back Q{int(q*100)}。")
                break

            base0_new = _quantile_drop0(nonpeak, q)
            if base0_new is None:
                base0_new = float(np.quantile(nonpeak, q))
            base0 = base0_new
            thr = base0 + thr_k * Noise

            cand = (~mask_peak) & (y > thr)
            if not np.any(cand):
                logs.append(f"Ture {it+1}：no new peaks")
                break

            regions = _contiguous_true_regions(cand)
            expanded = mask_peak.copy()
            for s, e in regions:
                s2 = max(0, s - pad_pts)
                e2 = min(n - 1, e + pad_pts)
                expanded[s2:e2 + 1] = True

            added = int(np.count_nonzero(expanded) - np.count_nonzero(mask_peak))
            mask_peak = expanded
            logs.append(f"Ture {it+1}：new mark points {added}")

        # ========== Final baseline ==========
        nonpeak_final = y[~mask_peak]
        if nonpeak_final.size < max(10, int(min_nonpeak_frac * n)):
            baseline = _quantile_drop0(y, q)
            logs.append("baseline")
        else:
            baseline = _quantile_drop0(nonpeak_final, q)
            logs.append("baseline")

        if baseline is None:
            baseline = float(np.quantile(y, q))
            logs.append("baseline poor points")

        info = {
            "baseline": float(baseline),
            "base0": float(base0),
            "threshold": float(thr),
            "q": q,
            "thr_k": thr_k,
            "pad_s": pad_s,
            "pad_points": int(pad_pts),
            "mask_peak_fraction": float(np.mean(mask_peak)),
            "nonpeak_points": int(nonpeak_final.size),
            "total_points": int(n),
            "drop_zeros": bool(drop_zeros),
            "logs": "；".join(logs),
        }
        return float(baseline)
    
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
        def generate_and_merge_intervals(elements, x):
            counts = [1]*len(elements)
            # 1. 生成初始区间
            intervals = []
            for element, count in zip(elements, counts):
                center = element
                intervals.append((center - x, center + x, count))
            
            # 2. 按区间起始点排序
            intervals.sort(key=lambda interval: interval[0])
            
            # 3. 合并重叠区间
            merged_intervals = []
            element_counts = []
            
            for interval in intervals:
                start, end, count = interval
                
                # 如果是第一个区间，直接添加
                if not merged_intervals:
                    merged_intervals.append((start, end))
                    element_counts.append(count)
                else:
                    # 检查当前区间是否与最后一个合并区间重叠
                    last_start, last_end = merged_intervals[-1]
                    
                    # 判断是否重叠：当前区间的起始点 <= 上一个区间的结束点
                    if start <= last_end:
                        # 合并区间，取更大的结束点
                        merged_intervals[-1] = (last_start, max(last_end, end))
                        # 合并计数
                        element_counts[-1] += count
                    else:
                        # 不重叠，添加新区间
                        merged_intervals.append((start, end))
                        element_counts.append(count)
            
            return merged_intervals, element_counts
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
            self.TargetList['Plot_RT_1'] = ''
            self.TargetList['Plot_Int_1'] = ''
            self.TargetList['Plot_RT_All_1'] = ''
            self.TargetList['Plot_Int_All_1'] = ''
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
                MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList_set.loc[i,'1_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList_set.loc[i,'1_Q3'])<=self.get_param('MS1_Tor'))&(self.Auto_Peak['Polarity']==self.TargetList_set.loc[i,'Polarity'])]
                if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                    MRM_2 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList_set.loc[i,'2_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList_set.loc[i,'2_Q3'])<=self.get_param('MS1_Tor'))&(self.Auto_Peak['Polarity']==self.TargetList_set.loc[i,'Polarity'])]
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
                    for i_plot in self.TargetList_set.loc[i,'TL_Index']:
                        self.TargetList.at[i_plot,'Plot_RT_All_1'] = Origin_RT_List_1
                        self.TargetList.at[i_plot,'Plot_Int_All_1'] = Origin_Int_List_1
                    ''' MRM-1 '''
                    Noise_1 = self.Calculate_Noise(Origin_RT_List_1,Origin_Int_List_1)
                    Baseline_1 = MRMProcess.calc_global_baseline_with_mask_index(Origin_RT_List_1,Origin_Int_List_1,[],[],Noise_1,q=0.20,thr_k=4.0,pad_s=1.5,pad_points=None,min_nonpeak_frac=0.05,max_iter=3,drop_zeros=True)
                    time_seq,Count_seq = generate_and_merge_intervals([x*60 for x in self.TargetList_set.loc[i,'RT']],self.get_param('Target_RT_Tor')+12)
                    All_Peak_Begin_1 = []
                    All_Peak_End_1 = []
                    All_Peak_Top_1 = []
                    All_Area_List_1 = []
                    All_Int_List_1 = []
                    All_RT_List_1 = []
                    All_RT_L_1 = []
                    All_RT_R_1 = []
                    All_Plot_RT_1 = []
                    All_Plot_Int_1 = []
                    for i_seq in range(len(time_seq)):
                        Index_forDetect = [x for x in range(len(Origin_RT_List_1)) if time_seq[i_seq][0]<=Origin_RT_List_1[x]<=time_seq[i_seq][1]]
                        RT_List_forDetect = [Origin_RT_List_1[x] for x in range(len(Origin_RT_List_1)) if time_seq[i_seq][0]<=Origin_RT_List_1[x]<=time_seq[i_seq][1]]
                        Int_List_forDetect = [Origin_Int_List_1[x] for x in range(len(Origin_RT_List_1)) if time_seq[i_seq][0]<=Origin_RT_List_1[x]<=time_seq[i_seq][1]]
                        Peak_Result_1 = self.MRMPeakDetecter(RT_List_forDetect,Int_List_forDetect,Baseline_1)
                        (Peak_Begin_1,Peak_End_1,Peak_Top_1,Area_List_1,Int_List_1) = Peak_Result_1
                        if len(Peak_Begin_1) > 0 and len(Peak_Begin_1) < Count_seq[i_seq]:
                            Ex_Peak_Begin_1 = [0]+Peak_End_1
                            Ex_Peak_End_1 = Peak_Begin_1+[len(Int_List_forDetect)]
                            for v in range(len(Ex_Peak_Begin_1)):
                                if Ex_Peak_Begin_1[v]-Ex_Peak_End_1[v] >= self.get_param('Points') and max(Int_List_forDetect[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]) >= self.get_param('min_Int'):
                                    temp_RT = RT_List_forDetect[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]
                                    temp_Int = Int_List_forDetect[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]
                                    Peak_Result_1 = self.MRMPeakDetecter(temp_RT,temp_Int)
                                    (Ex_Peak_Begin_1,Ex_Peak_End_1,Ex_Peak_Top_1,Ex_Area_List_1,Ex_Int_List_1) = Peak_Result_1
                                    Ex_Time_Begin_1 = [temp_RT[x] for x in Ex_Peak_Begin_1]
                                    Ex_Time_End_1 = [temp_RT[x] for x in Ex_Peak_End_1]
                                    Ex_Time_Top_1 = [temp_RT[x] for x in Ex_Peak_Top_1]
                                    Ex_Peak_Begin_1 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Begin_1]
                                    Ex_Peak_End_1 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_End_1]
                                    Ex_Peak_Top_1 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Top_1]
                                    Peak_Begin_1 += Ex_Peak_Begin_1
                                    Peak_End_1 += Ex_Peak_End_1
                                    Peak_Top_1 += Ex_Peak_Top_1
                                    Area_List_1 += Ex_Area_List_1
                                    Int_List_1 += Ex_Int_List_1
                        Index_1 = [x for x in range(len(Int_List_1)) if Int_List_1[x]>=self.get_param('min_Int')]
                        Peak_Begin_1 = [Peak_Begin_1[x] for x in Index_1]
                        Peak_Begin_1 = [Index_forDetect[x] for x in Peak_Begin_1]
                        Peak_End_1 = [Peak_End_1[x] for x in Index_1]
                        Peak_End_1 = [Index_forDetect[x] for x in Peak_End_1]
                        Peak_Top_1 = [Peak_Top_1[x] for x in Index_1]
                        Peak_Top_1 = [Index_forDetect[x] for x in Peak_Top_1]
                        Area_List_1 = [Area_List_1[x] for x in Index_1]
                        Int_List_1 = [Int_List_1[x] for x in Index_1]
                        RT_List_1 = [Origin_RT_List_1[x] for x in Peak_Top_1]
                        RT_L_1 = [Origin_RT_List_1[x] for x in Peak_Begin_1]
                        RT_R_1 = [Origin_RT_List_1[x] for x in Peak_End_1]
                        Plot_RT_1 = [Origin_RT_List_1[Peak_Begin_1[x]:Peak_End_1[x]+1] for x in range(len(Peak_Begin_1))]
                        Plot_Int_1 = [Origin_Int_List_1[Peak_Begin_1[x]:Peak_End_1[x]+1] for x in range(len(Peak_Begin_1))]
                        All_Peak_Begin_1 = All_Peak_Begin_1 + Peak_Begin_1
                        All_Peak_End_1 = All_Peak_End_1 + Peak_End_1
                        All_Peak_Top_1 = All_Peak_Top_1 + Peak_Top_1
                        All_Area_List_1 = All_Area_List_1 + Area_List_1
                        All_Int_List_1 = All_Int_List_1 + Int_List_1
                        All_RT_List_1 = All_RT_List_1 + RT_List_1
                        All_RT_L_1 = All_RT_L_1 + RT_L_1
                        All_RT_R_1 = All_RT_R_1 + RT_R_1
                        All_Plot_RT_1 = All_Plot_RT_1 + Plot_RT_1
                        All_Plot_Int_1 = All_Plot_Int_1 + Plot_Int_1
                    Peak_Begin_1 = All_Peak_Begin_1
                    Peak_End_1 = All_Peak_End_1
                    Peak_Top_1 = All_Peak_Top_1
                    Area_List_1 = All_Area_List_1
                    Int_List_1 = All_Int_List_1
                    RT_List_1 = All_RT_List_1
                    RT_L_1 = All_RT_L_1
                    RT_R_1 = All_RT_R_1
                    Plot_RT_1 = All_Plot_RT_1
                    Plot_Int_1 = All_Plot_Int_1
                    ''' MRM-2 '''
                    Noise_2 = self.Calculate_Noise(Origin_RT_List_2,Origin_Int_List_2)
                    Baseline_2 = MRMProcess.calc_global_baseline_with_mask_index(Origin_RT_List_2,Origin_Int_List_2,[],[],Noise_2,q=0.20,thr_k=4.0,pad_s=1.5,pad_points=None,min_nonpeak_frac=0.05,max_iter=3,drop_zeros=True)
                    All_Peak_Begin_2 = []
                    All_Peak_End_2 = []
                    All_Peak_Top_2 = []
                    All_Area_List_2 = []
                    All_Int_List_2 = []
                    All_RT_List_2 = []
                    All_RT_L_2 = []
                    All_RT_R_2 = []
                    for i_seq in range(len(time_seq)):
                        Index_forDetect = [x for x in range(len(Origin_RT_List_2)) if time_seq[i_seq][0]<=Origin_RT_List_2[x]<=time_seq[i_seq][1]]
                        RT_List_forDetect = [Origin_RT_List_2[x] for x in range(len(Origin_RT_List_2)) if time_seq[i_seq][0]<=Origin_RT_List_2[x]<=time_seq[i_seq][1]]
                        Int_List_forDetect = [Origin_Int_List_2[x] for x in range(len(Origin_RT_List_2)) if time_seq[i_seq][0]<=Origin_RT_List_2[x]<=time_seq[i_seq][1]]
                        Peak_Result_2 = self.MRMPeakDetecter(RT_List_forDetect,Int_List_forDetect,Baseline_2)
                        (Peak_Begin_2,Peak_End_2,Peak_Top_2,Area_List_2,Int_List_2) = Peak_Result_2
                        if len(Peak_Begin_2) > 0 and len(Peak_Begin_2) < Count_seq[i_seq]:
                            Ex_Peak_Begin_2 = [0]+Peak_End_2
                            Ex_Peak_End_2 = Peak_Begin_2+[len(Int_List_forDetect)]
                            for v in range(len(Ex_Peak_Begin_2)):
                                if Ex_Peak_Begin_2[v]-Ex_Peak_End_2[v] >= self.get_param('Points') and max(Int_List_forDetect[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]) >= self.get_param('min_Int'):
                                    temp_RT = RT_List_forDetect[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]
                                    temp_Int = Int_List_forDetect[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]
                                    Peak_Result_2 = self.MRMPeakDetecter(temp_RT,temp_Int)
                                    (Ex_Peak_Begin_2,Ex_Peak_End_2,Ex_Peak_Top_2,Ex_Area_List_2,Ex_Int_List_2) = Peak_Result_2
                                    Ex_Time_Begin_2 = [temp_RT[x] for x in Ex_Peak_Begin_2]
                                    Ex_Time_End_2 = [temp_RT[x] for x in Ex_Peak_End_2]
                                    Ex_Time_Top_2 = [temp_RT[x] for x in Ex_Peak_Top_2]
                                    Ex_Peak_Begin_2 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Begin_2]
                                    Ex_Peak_End_2 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_End_2]
                                    Ex_Peak_Top_2 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Top_2]
                                    Peak_Begin_2 += Ex_Peak_Begin_2
                                    Peak_End_2 += Ex_Peak_End_2
                                    Peak_Top_2 += Ex_Peak_Top_2
                                    Area_List_2 += Ex_Area_List_2
                                    Int_List_2 += Ex_Int_List_2
                        Index_2 = [x for x in range(len(Int_List_2)) if Int_List_2[x]>=self.get_param('min_Int')]
                        Peak_Begin_2 = [Peak_Begin_2[x] for x in Index_2]
                        Peak_Begin_2 = [Index_forDetect[x] for x in Peak_Begin_2]
                        Peak_End_2 = [Peak_End_2[x] for x in Index_2]
                        Peak_End_2 = [Index_forDetect[x] for x in Peak_End_2]
                        Peak_Top_2 = [Peak_Top_2[x] for x in Index_2]
                        Peak_Top_2 = [Index_forDetect[x] for x in Peak_Top_2]
                        Area_List_2 = [Area_List_2[x] for x in Index_2]
                        Int_List_2 = [Int_List_2[x] for x in Index_2]
                        RT_List_2 = [Origin_RT_List_2[x] for x in Peak_Top_2]
                        RT_L_2 = [Origin_RT_List_2[x] for x in Peak_Begin_2]
                        RT_R_2 = [Origin_RT_List_2[x] for x in Peak_End_2]
                        All_Peak_Begin_2 = All_Peak_Begin_2 + Peak_Begin_2
                        All_Peak_End_2 = All_Peak_End_2 + Peak_End_2
                        All_Peak_Top_2 = All_Peak_Top_2 + Peak_Top_2
                        All_Area_List_2 = All_Area_List_2 + Area_List_2
                        All_Int_List_2 = All_Int_List_2 + Int_List_2
                        All_RT_List_2 = All_RT_List_2 + RT_List_2
                        All_RT_L_2 = All_RT_L_2 + RT_L_2
                        All_RT_R_2 = All_RT_R_2 + RT_R_2
                    Peak_Begin_2 = All_Peak_Begin_2
                    Peak_End_2 = All_Peak_End_2
                    Peak_Top_2 = All_Peak_Top_2
                    Area_List_2 = All_Area_List_2
                    Int_List_2 = All_Int_List_2
                    RT_List_2 = All_RT_List_2
                    RT_L_2 = All_RT_L_2
                    RT_R_2 = All_RT_R_2
                    Score_matrix= np.zeros([len(RT_List_1),len(RT_List_2)])
                    for i_1 in range(len(RT_List_1)):
                        for i_2 in range(len(RT_List_2)):
                            if abs(RT_List_1[i_1]-RT_List_2[i_2])<self.get_param('Target_RT_Tor'):
                                Score_matrix[i_1,i_2] = 1-abs(RT_List_1[i_1]-RT_List_2[i_2])/self.get_param('Target_RT_Tor')
                    Sm_1,Sm_2 = linear_sum_assignment(Score_matrix,True)
                    Sm_Filter = [(x,y) for x,y in zip(Sm_1,Sm_2) if Score_matrix[x,y]>0]
                    
                    Score_matrix= np.zeros([len(Sm_Filter),len(self.TargetList_set.loc[i,'RT'])])
                    for i_1 in range(len(Sm_Filter)):
                        for i_2 in range(len(self.TargetList_set.loc[i,'RT'])):
                            if abs(RT_List_1[Sm_Filter[i_1][0]]-self.TargetList_set.loc[i,'RT'][i_2]*60)< self.get_param('Target_RT_Tor') and abs(RT_List_2[Sm_Filter[i_1][1]]-self.TargetList_set.loc[i,'RT'][i_2]*60)< self.get_param('Target_RT_Tor') and Peak_End_1[Sm_Filter[i_1][0]]-Peak_Begin_1[Sm_Filter[i_1][0]] >= self.get_param('Points')-1:
                                Score_matrix[i_1,i_2] = 0.3*(2-abs(RT_List_1[Sm_Filter[i_1][0]]-self.TargetList_set.loc[i,'RT'][i_2]*60)/ self.get_param('Target_RT_Tor') - abs(RT_List_2[Sm_Filter[i_1][1]]-self.TargetList_set.loc[i,'RT'][i_2]*60)/ self.get_param('Target_RT_Tor')) + 0.25*(Int_List_1[Sm_Filter[i_1][0]]/max(Int_List_1) + Int_List_2[Sm_Filter[i_1][1]]/max(Int_List_2))
                    Sm_1,Sm_2 = linear_sum_assignment(Score_matrix,True)         
                    Sm_Assign = [(x,y) for x,y in zip(Sm_1,Sm_2) if Score_matrix[x,y]>0]
                    for i_1,i_2 in Sm_Assign:
                        IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS']]['Area_1'])
                        if len(IS_Area)>0:
                            IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS']]['Area_1'])[0]
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS_Area'] = int(IS_Area)
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_1/IS'] = np.around(Area_List_1[Sm_Filter[i_1][0]]/IS_Area,5)
                        else:
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_1/IS'] = 'None'
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS_Area'] = 'None'
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'RT_1'] = np.around(RT_List_1[Sm_Filter[i_1][0]]/60,3)
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_1'] = int(Area_List_1[Sm_Filter[i_1][0]])
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Hight_1'] = int(Int_List_1[Sm_Filter[i_1][0]])
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Noise_1'] = int(Noise_1)
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_left_1'] = np.around(RT_L_1[Sm_Filter[i_1][0]]/60,4)
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_right_1'] = np.around(RT_R_1[Sm_Filter[i_1][0]]/60,4)
                        self.TargetList.at[self.TargetList_set.at[i,'TL_Index'][i_2],'Plot_RT_1'] = Plot_RT_1[Sm_Filter[i_1][0]]
                        self.TargetList.at[self.TargetList_set.at[i,'TL_Index'][i_2],'Plot_Int_1'] = Plot_Int_1[Sm_Filter[i_1][0]]
                        self.TargetList.at[self.TargetList_set.at[i,'TL_Index'][i_2],'Plot_RT_All_1'] = Origin_RT_List_1
                        self.TargetList.at[self.TargetList_set.at[i,'TL_Index'][i_2],'Plot_Int_All_1'] = Origin_Int_List_1
                        if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'RT_2'] = np.around(RT_List_2[Sm_Filter[i_1][1]]/60,3)
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_2'] = int(Area_List_2[Sm_Filter[i_1][1]])
                            IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS']]['Area_1'])
                            if len(IS_Area)>0:
                                IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS']]['Area_1'])[0]
                                self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_2/IS'] = np.around(Area_List_2[Sm_Filter[i_1][1]]/IS_Area,4)
                            else:
                                self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_2/IS'] = 'None'
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Hight_2'] = int(Int_List_2[Sm_Filter[i_1][1]])
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Noise_2'] = int(Noise_2)
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_left_2'] = np.around(RT_L_2[Sm_Filter[i_1][1]]/60,4)
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_right_2'] = np.around(RT_R_2[Sm_Filter[i_1][1]]/60,4)
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
                        fig_1.write_html(FilePath+'/'+'Q1-'+str(np.around(self.TargetList_set.loc[i,'1_Q1'],1))+' Q3-'+str(np.around(self.TargetList_set.loc[i,'1_Q3'],1))+'.html',config={'responsive': False})
                        if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                            fig_2.write_html(FilePath+'/'+'Q1-'+str(np.around(self.TargetList_set.loc[i,'2_Q1'],1))+' Q3-'+str(np.around(self.TargetList_set.loc[i,'2_Q3'],1))+'.html',config={'responsive': False})
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
                    MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList.loc[i,'1_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList.loc[i,'1_Q3'])<=self.get_param('MS1_Tor'))]
                    if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                        MRM_2 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList_set.loc[i,'2_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList_set.loc[i,'2_Q3'])<=self.get_param('MS1_Tor'))]
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
    
    def detect_Peak_MRM_only(self,FilePath=''):
        def generate_and_merge_intervals(elements, x):
            counts = [1]*len(elements)
            # 1. 生成初始区间
            intervals = []
            for element, count in zip(elements, counts):
                center = element
                intervals.append((center - x, center + x, count))
            
            # 2. 按区间起始点排序
            intervals.sort(key=lambda interval: interval[0])
            
            # 3. 合并重叠区间
            merged_intervals = []
            element_counts = []
            
            for interval in intervals:
                start, end, count = interval
                
                # 如果是第一个区间，直接添加
                if not merged_intervals:
                    merged_intervals.append((start, end))
                    element_counts.append(count)
                else:
                    # 检查当前区间是否与最后一个合并区间重叠
                    last_start, last_end = merged_intervals[-1]
                    
                    # 判断是否重叠：当前区间的起始点 <= 上一个区间的结束点
                    if start <= last_end:
                        # 合并区间，取更大的结束点
                        merged_intervals[-1] = (last_start, max(last_end, end))
                        # 合并计数
                        element_counts[-1] += count
                    else:
                        # 不重叠，添加新区间
                        merged_intervals.append((start, end))
                        element_counts.append(count)
            
            return merged_intervals, element_counts
        self.TargetList['IS_Area'] = ''
        self.TargetList['RT_1'] = ''
        self.TargetList['Area_1'] = ''
        self.TargetList['Area_1/IS'] = ''
        self.TargetList['Hight_1'] = ''
        self.TargetList['Noise_1'] = ''
        self.TargetList['Edge_left_1'] = ''
        self.TargetList['Edge_right_1'] = ''
        self.TargetList['Plot_RT_1'] = ''
        self.TargetList['Plot_Int_1'] = ''
        self.TargetList['Plot_RT_All_1'] = ''
        self.TargetList['Plot_Int_All_1'] = ''
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
            MRM_1 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList_set.loc[i,'1_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList_set.loc[i,'1_Q3'])<=self.get_param('MS1_Tor'))&(self.Auto_Peak['Polarity']==self.TargetList_set.loc[i,'Polarity'])]
            if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                MRM_2 = self.Auto_Peak[(abs(self.Auto_Peak['Pre_MZ']-self.TargetList_set.loc[i,'2_Q1'])<=self.get_param('MS1_Tor'))&(abs(self.Auto_Peak['Pro_MZ']-self.TargetList_set.loc[i,'2_Q3'])<=self.get_param('MS1_Tor'))&(self.Auto_Peak['Polarity']==self.TargetList_set.loc[i,'Polarity'])]
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
                ''' MRM-1 '''
                Noise_1 = self.Calculate_Noise(Origin_RT_List_1,Origin_Int_List_1)
                Baseline_1 = MRMProcess.calc_global_baseline_with_mask_index(Origin_RT_List_1,Origin_Int_List_1,[],[],Noise_1,q=0.20,thr_k=4.0,pad_s=1.5,pad_points=None,min_nonpeak_frac=0.05,max_iter=3,drop_zeros=True)
                time_seq,Count_seq = generate_and_merge_intervals([x*60 for x in self.TargetList_set.loc[i,'RT']],self.get_param('Target_RT_Tor')+12)
                All_Peak_Begin_1 = []
                All_Peak_End_1 = []
                All_Peak_Top_1 = []
                All_Area_List_1 = []
                All_Int_List_1 = []
                All_RT_List_1 = []
                All_RT_L_1 = []
                All_RT_R_1 = []
                All_Plot_RT_1 = []
                All_Plot_Int_1 = []
                for i_seq in range(len(time_seq)):
                    Index_forDetect = [x for x in range(len(Origin_RT_List_1)) if time_seq[i_seq][0]<=Origin_RT_List_1[x]<=time_seq[i_seq][1]]
                    RT_List_forDetect = [Origin_RT_List_1[x] for x in range(len(Origin_RT_List_1)) if time_seq[i_seq][0]<=Origin_RT_List_1[x]<=time_seq[i_seq][1]]
                    Int_List_forDetect = [Origin_Int_List_1[x] for x in range(len(Origin_RT_List_1)) if time_seq[i_seq][0]<=Origin_RT_List_1[x]<=time_seq[i_seq][1]]
                    Peak_Result_1 = self.MRMPeakDetecter(RT_List_forDetect,Int_List_forDetect,Baseline_1)
                    (Peak_Begin_1,Peak_End_1,Peak_Top_1,Area_List_1,Int_List_1) = Peak_Result_1
                    if len(Peak_Begin_1) > 0 and len(Peak_Begin_1) < Count_seq[i_seq]:
                        Ex_Peak_Begin_1 = [0]+Peak_End_1
                        Ex_Peak_End_1 = Peak_Begin_1+[len(Int_List_forDetect)]
                        for v in range(len(Ex_Peak_Begin_1)):
                            if Ex_Peak_Begin_1[v]-Ex_Peak_End_1[v] >= self.get_param('Points') and max(Int_List_forDetect[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]) >= self.get_param('min_Int'):
                                temp_RT = RT_List_forDetect[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]
                                temp_Int = Int_List_forDetect[Ex_Peak_Begin_1[v]:Ex_Peak_End_1[v]]
                                Peak_Result_1 = self.MRMPeakDetecter(temp_RT,temp_Int)
                                (Ex_Peak_Begin_1,Ex_Peak_End_1,Ex_Peak_Top_1,Ex_Area_List_1,Ex_Int_List_1) = Peak_Result_1
                                Ex_Time_Begin_1 = [temp_RT[x] for x in Ex_Peak_Begin_1]
                                Ex_Time_End_1 = [temp_RT[x] for x in Ex_Peak_End_1]
                                Ex_Time_Top_1 = [temp_RT[x] for x in Ex_Peak_Top_1]
                                Ex_Peak_Begin_1 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Begin_1]
                                Ex_Peak_End_1 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_End_1]
                                Ex_Peak_Top_1 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Top_1]
                                Peak_Begin_1 += Ex_Peak_Begin_1
                                Peak_End_1 += Ex_Peak_End_1
                                Peak_Top_1 += Ex_Peak_Top_1
                                Area_List_1 += Ex_Area_List_1
                                Int_List_1 += Ex_Int_List_1
                    Index_1 = [x for x in range(len(Int_List_1)) if Int_List_1[x]>=self.get_param('min_Int')]
                    Peak_Begin_1 = [Peak_Begin_1[x] for x in Index_1]
                    Peak_Begin_1 = [Index_forDetect[x] for x in Peak_Begin_1]
                    Peak_End_1 = [Peak_End_1[x] for x in Index_1]
                    Peak_End_1 = [Index_forDetect[x] for x in Peak_End_1]
                    Peak_Top_1 = [Peak_Top_1[x] for x in Index_1]
                    Peak_Top_1 = [Index_forDetect[x] for x in Peak_Top_1]
                    Area_List_1 = [Area_List_1[x] for x in Index_1]
                    Int_List_1 = [Int_List_1[x] for x in Index_1]
                    RT_List_1 = [Origin_RT_List_1[x] for x in Peak_Top_1]
                    RT_L_1 = [Origin_RT_List_1[x] for x in Peak_Begin_1]
                    RT_R_1 = [Origin_RT_List_1[x] for x in Peak_End_1]
                    Plot_RT_1 = [Origin_RT_List_1[Peak_Begin_1[x]:Peak_End_1[x]+1] for x in range(len(Peak_Begin_1))]
                    Plot_Int_1 = [Origin_Int_List_1[Peak_Begin_1[x]:Peak_End_1[x]+1] for x in range(len(Peak_Begin_1))]
                    All_Peak_Begin_1 = All_Peak_Begin_1 + Peak_Begin_1
                    All_Peak_End_1 = All_Peak_End_1 + Peak_End_1
                    All_Peak_Top_1 = All_Peak_Top_1 + Peak_Top_1
                    All_Area_List_1 = All_Area_List_1 + Area_List_1
                    All_Int_List_1 = All_Int_List_1 + Int_List_1
                    All_RT_List_1 = All_RT_List_1 + RT_List_1
                    All_RT_L_1 = All_RT_L_1 + RT_L_1
                    All_RT_R_1 = All_RT_R_1 + RT_R_1
                    All_Plot_RT_1 = All_Plot_RT_1 + Plot_RT_1
                    All_Plot_Int_1 = All_Plot_Int_1 + Plot_Int_1
                Peak_Begin_1 = All_Peak_Begin_1
                Peak_End_1 = All_Peak_End_1
                Peak_Top_1 = All_Peak_Top_1
                Area_List_1 = All_Area_List_1
                Int_List_1 = All_Int_List_1
                RT_List_1 = All_RT_List_1
                RT_L_1 = All_RT_L_1
                RT_R_1 = All_RT_R_1
                Plot_RT_1 = All_Plot_RT_1
                Plot_Int_1 = All_Plot_Int_1
                ''' MRM-2 '''
                Noise_2 = self.Calculate_Noise(Origin_RT_List_2,Origin_Int_List_2)
                Baseline_2 = MRMProcess.calc_global_baseline_with_mask_index(Origin_RT_List_2,Origin_Int_List_2,[],[],Noise_2,q=0.20,thr_k=4.0,pad_s=1.5,pad_points=None,min_nonpeak_frac=0.05,max_iter=3,drop_zeros=True)
                All_Peak_Begin_2 = []
                All_Peak_End_2 = []
                All_Peak_Top_2 = []
                All_Area_List_2 = []
                All_Int_List_2 = []
                All_RT_List_2 = []
                All_RT_L_2 = []
                All_RT_R_2 = []
                for i_seq in range(len(time_seq)):
                    Index_forDetect = [x for x in range(len(Origin_RT_List_2)) if time_seq[i_seq][0]<=Origin_RT_List_2[x]<=time_seq[i_seq][1]]
                    RT_List_forDetect = [Origin_RT_List_2[x] for x in range(len(Origin_RT_List_2)) if time_seq[i_seq][0]<=Origin_RT_List_2[x]<=time_seq[i_seq][1]]
                    Int_List_forDetect = [Origin_Int_List_2[x] for x in range(len(Origin_RT_List_2)) if time_seq[i_seq][0]<=Origin_RT_List_2[x]<=time_seq[i_seq][1]]
                    Peak_Result_2 = self.MRMPeakDetecter(RT_List_forDetect,Int_List_forDetect,Baseline_2)
                    (Peak_Begin_2,Peak_End_2,Peak_Top_2,Area_List_2,Int_List_2) = Peak_Result_2
                    if len(Peak_Begin_2) > 0 and len(Peak_Begin_2) < Count_seq[i_seq]:
                        Ex_Peak_Begin_2 = [0]+Peak_End_2
                        Ex_Peak_End_2 = Peak_Begin_2+[len(Int_List_forDetect)]
                        for v in range(len(Ex_Peak_Begin_2)):
                            if Ex_Peak_Begin_2[v]-Ex_Peak_End_2[v] >= self.get_param('Points') and max(Int_List_forDetect[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]) >= self.get_param('min_Int'):
                                temp_RT = RT_List_forDetect[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]
                                temp_Int = Int_List_forDetect[Ex_Peak_Begin_2[v]:Ex_Peak_End_2[v]]
                                Peak_Result_2 = self.MRMPeakDetecter(temp_RT,temp_Int)
                                (Ex_Peak_Begin_2,Ex_Peak_End_2,Ex_Peak_Top_2,Ex_Area_List_2,Ex_Int_List_2) = Peak_Result_2
                                Ex_Time_Begin_2 = [temp_RT[x] for x in Ex_Peak_Begin_2]
                                Ex_Time_End_2 = [temp_RT[x] for x in Ex_Peak_End_2]
                                Ex_Time_Top_2 = [temp_RT[x] for x in Ex_Peak_Top_2]
                                Ex_Peak_Begin_2 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Begin_2]
                                Ex_Peak_End_2 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_End_2]
                                Ex_Peak_Top_2 = [x for x in range(len(RT_List_forDetect)) if RT_List_forDetect[x] in Ex_Time_Top_2]
                                Peak_Begin_2 += Ex_Peak_Begin_2
                                Peak_End_2 += Ex_Peak_End_2
                                Peak_Top_2 += Ex_Peak_Top_2
                                Area_List_2 += Ex_Area_List_2
                                Int_List_2 += Ex_Int_List_2
                    Index_2 = [x for x in range(len(Int_List_2)) if Int_List_2[x]>=self.get_param('min_Int')]
                    Peak_Begin_2 = [Peak_Begin_2[x] for x in Index_2]
                    Peak_Begin_2 = [Index_forDetect[x] for x in Peak_Begin_2]
                    Peak_End_2 = [Peak_End_2[x] for x in Index_2]
                    Peak_End_2 = [Index_forDetect[x] for x in Peak_End_2]
                    Peak_Top_2 = [Peak_Top_2[x] for x in Index_2]
                    Peak_Top_2 = [Index_forDetect[x] for x in Peak_Top_2]
                    Area_List_2 = [Area_List_2[x] for x in Index_2]
                    Int_List_2 = [Int_List_2[x] for x in Index_2]
                    RT_List_2 = [Origin_RT_List_2[x] for x in Peak_Top_2]
                    RT_L_2 = [Origin_RT_List_2[x] for x in Peak_Begin_2]
                    RT_R_2 = [Origin_RT_List_2[x] for x in Peak_End_2]
                    All_Peak_Begin_2 = All_Peak_Begin_2 + Peak_Begin_2
                    All_Peak_End_2 = All_Peak_End_2 + Peak_End_2
                    All_Peak_Top_2 = All_Peak_Top_2 + Peak_Top_2
                    All_Area_List_2 = All_Area_List_2 + Area_List_2
                    All_Int_List_2 = All_Int_List_2 + Int_List_2
                    All_RT_List_2 = All_RT_List_2 + RT_List_2
                    All_RT_L_2 = All_RT_L_2 + RT_L_2
                    All_RT_R_2 = All_RT_R_2 + RT_R_2
                Peak_Begin_2 = All_Peak_Begin_2
                Peak_End_2 = All_Peak_End_2
                Peak_Top_2 = All_Peak_Top_2
                Area_List_2 = All_Area_List_2
                Int_List_2 = All_Int_List_2
                RT_List_2 = All_RT_List_2
                RT_L_2 = All_RT_L_2
                RT_R_2 = All_RT_R_2
                Score_matrix= np.zeros([len(RT_List_1),len(RT_List_2)])
                for i_1 in range(len(RT_List_1)):
                    for i_2 in range(len(RT_List_2)):
                        if abs(RT_List_1[i_1]-RT_List_2[i_2])<self.get_param('Target_RT_Tor'):
                            Score_matrix[i_1,i_2] = 1-abs(RT_List_1[i_1]-RT_List_2[i_2])/self.get_param('Target_RT_Tor')
                Sm_1,Sm_2 = linear_sum_assignment(Score_matrix,True)
                Sm_Filter = [(x,y) for x,y in zip(Sm_1,Sm_2) if Score_matrix[x,y]>0]
                
                Score_matrix= np.zeros([len(Sm_Filter),len(self.TargetList_set.loc[i,'RT'])])
                for i_1 in range(len(Sm_Filter)):
                    for i_2 in range(len(self.TargetList_set.loc[i,'RT'])):
                        if abs(RT_List_1[Sm_Filter[i_1][0]]-self.TargetList_set.loc[i,'RT'][i_2]*60)< self.get_param('Target_RT_Tor') and abs(RT_List_2[Sm_Filter[i_1][1]]-self.TargetList_set.loc[i,'RT'][i_2]*60)< self.get_param('Target_RT_Tor') and Peak_End_1[Sm_Filter[i_1][0]]-Peak_Begin_1[Sm_Filter[i_1][0]] >= self.get_param('Points')-1:
                            Score_matrix[i_1,i_2] = 0.3*(2-abs(RT_List_1[Sm_Filter[i_1][0]]-self.TargetList_set.loc[i,'RT'][i_2]*60)/ self.get_param('Target_RT_Tor') - abs(RT_List_2[Sm_Filter[i_1][1]]-self.TargetList_set.loc[i,'RT'][i_2]*60)/ self.get_param('Target_RT_Tor')) + 0.25*(Int_List_1[Sm_Filter[i_1][0]]/max(Int_List_1) + Int_List_2[Sm_Filter[i_1][1]]/max(Int_List_2))
                Sm_1,Sm_2 = linear_sum_assignment(Score_matrix,True)         
                Sm_Assign = [(x,y) for x,y in zip(Sm_1,Sm_2) if Score_matrix[x,y]>0]
                for i_1,i_2 in Sm_Assign:
                    self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_1/IS'] = 'None'
                    self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'IS_Area'] = 'None'
                    self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'RT_1'] = np.around(RT_List_1[Sm_Filter[i_1][0]]/60,3)
                    self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_1'] = int(Area_List_1[Sm_Filter[i_1][0]])
                    self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Hight_1'] = int(Int_List_1[Sm_Filter[i_1][0]])
                    self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Noise_1'] = int(Noise_1)
                    self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_left_1'] = np.around(RT_L_1[Sm_Filter[i_1][0]]/60,4)
                    self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_right_1'] = np.around(RT_R_1[Sm_Filter[i_1][0]]/60,4)
                    self.TargetList.at[self.TargetList_set.at[i,'TL_Index'][i_2],'Plot_RT_1'] = Plot_RT_1[Sm_Filter[i_1][0]]
                    self.TargetList.at[self.TargetList_set.at[i,'TL_Index'][i_2],'Plot_Int_1'] = Plot_Int_1[Sm_Filter[i_1][0]]
                    self.TargetList.at[self.TargetList_set.at[i,'TL_Index'][i_2],'Plot_RT_All_1'] = Origin_RT_List_1
                    self.TargetList.at[self.TargetList_set.at[i,'TL_Index'][i_2],'Plot_Int_All_1'] = Origin_Int_List_1
                    if pd.notna(self.TargetList_set.loc[i,'2_Q1']):
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'RT_2'] = np.around(RT_List_2[Sm_Filter[i_1][1]]/60,3)
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_2'] = int(Area_List_2[Sm_Filter[i_1][1]])
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Area_2/IS'] = 'None'
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Hight_2'] = int(Int_List_2[Sm_Filter[i_1][1]])
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Noise_2'] = int(Noise_2)
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_left_2'] = np.around(RT_L_2[Sm_Filter[i_1][1]]/60,4)
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Edge_right_2'] = np.around(RT_R_2[Sm_Filter[i_1][1]]/60,4)
                    x_1=list(MRM_1.loc[:,'RTList']/60)[0][Peak_Begin_1[Sm_Filter[i_1][0]]:Peak_End_1[Sm_Filter[i_1][0]]+1]
                    y_1=list(MRM_1.loc[:,'Int'])[0][Peak_Begin_1[Sm_Filter[i_1][0]]:Peak_End_1[Sm_Filter[i_1][0]]+1]
                    x_2=list(MRM_2.loc[:,'RTList']/60)[0][Peak_Begin_2[Sm_Filter[i_1][1]]:Peak_End_2[Sm_Filter[i_1][1]]+1]
                    y_2=list(MRM_2.loc[:,'Int'])[0][Peak_Begin_2[Sm_Filter[i_1][1]]:Peak_End_2[Sm_Filter[i_1][1]]+1]
                    y_1 = [y_1[x] for x in range(len(x_1)) if x_1[x] in x_2]
                    y_2 = [y_2[x] for x in range(len(x_2)) if x_2[x] in x_1]
                    self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][i_2],'Similarity'] = np.around(MRMProcess.CosineSimilarity(y_1, y_2, LOG=False),3)
        
    def detect_Blank(self,FilePath=''):
        def integral(x,y):
            x_diff = list(np.diff(x))
            y_mix = [0.5*(y[x]+y[x+1]) for x in range(len(y)-1)]
            return sum([x_diff[x]*y_mix[x] for x in range(len(x_diff))])
        def generate_and_merge_intervals(elements, x):
            counts = [1]*len(elements)
            # 1. 生成初始区间
            intervals = []
            for element, count in zip(elements, counts):
                center = element
                intervals.append((center - x, center + x, count))
            
            # 2. 按区间起始点排序
            intervals.sort(key=lambda interval: interval[0])
            
            # 3. 合并重叠区间
            merged_intervals = []
            element_counts = []
            
            for interval in intervals:
                start, end, count = interval
                
                # 如果是第一个区间，直接添加
                if not merged_intervals:
                    merged_intervals.append((start, end))
                    element_counts.append(count)
                else:
                    # 检查当前区间是否与最后一个合并区间重叠
                    last_start, last_end = merged_intervals[-1]
                    
                    # 判断是否重叠：当前区间的起始点 <= 上一个区间的结束点
                    if start <= last_end:
                        # 合并区间，取更大的结束点
                        merged_intervals[-1] = (last_start, max(last_end, end))
                        # 合并计数
                        element_counts[-1] += count
                    else:
                        # 不重叠，添加新区间
                        merged_intervals.append((start, end))
                        element_counts.append(count)
            
            return merged_intervals, element_counts
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
                Baseline_1 = MRMProcess.calc_global_baseline_with_mask_index(Origin_RT_List_1,Origin_Int_List_1,[],[],Noise_1,q=0.20,thr_k=4.0,pad_s=1.5,pad_points=None,min_nonpeak_frac=0.05,max_iter=3,drop_zeros=True)
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
                for ii in range(len(self.TargetList_set.loc[i,'RT'])):
                    Blank_RT = [Origin_RT_List_1[x] for x in range(len(Origin_RT_List_1)) if self.TargetList_set.loc[i,'RT'][ii]*60-self.get_param('Target_RT_Tor') <= Origin_RT_List_1[x] <= self.TargetList_set.loc[i,'RT'][ii]*60+self.get_param('Target_RT_Tor')]
                    Blank_Int = [Origin_Int_List_1[x] for x in range(len(Origin_RT_List_1)) if self.TargetList_set.loc[i,'RT'][ii]*60-self.get_param('Target_RT_Tor') <= Origin_RT_List_1[x] <= self.TargetList_set.loc[i,'RT'][ii]*60+self.get_param('Target_RT_Tor')]
                    if len(Blank_RT) == 0:
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'Area_1']=0
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'Area_1/IS']=0
                    else:
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'Area_1']=max(0,integral(Blank_RT,Blank_Int)-0.5*(min(Blank_Int[0],Baseline_1)+min(Blank_Int[-1],Baseline_1))*(Blank_RT[-1]-Blank_RT[0]))
                        IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'IS']]['Area_1'])
                        if len(IS_Area)>0:
                            IS_Area = list(self.Calibrant[self.Calibrant['Identity keys']==self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'IS']]['Area_1'])[0]
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'Area_1/IS'] = np.around(self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'Area_1']/IS_Area,4)
                        else:
                            self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'Area_1/IS'] = ''
                        self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'Hight_1'] = max(Blank_Int)
                    if len(FilePath)>0:
                        fig_1.add_trace(
                            go.Scatter(
                                x=[x/60 for x in Blank_RT],
                                y=Blank_Int,
                                mode='lines',
                                name=self.TargetList.loc[self.TargetList_set.loc[i,'TL_Index'][ii],'Identity keys'],
                                line={'width':1},
                                fill='tozeroy',
                                ))
                if len(FilePath)>0:
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
                    fig_1.write_html(FilePath+'/'+'Q1-'+str(np.around(self.TargetList_set.loc[i,'1_Q1'],1))+' Q3-'+str(np.around(self.TargetList_set.loc[i,'1_Q3'],1))+'.html',config={'responsive': False})
        if len(FilePath)>0:
            self.TargetList.to_excel(FilePath+'/MRM-Results-'+str(time.localtime()[1])+str(time.localtime()[2])+str(time.localtime()[3])+str(time.localtime()[4])+'.xlsx',index=False)
                    
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
        if 2<= x <= len(Auto_Int_List)-3:
            FirstDerivative = (Auto_Int_List[x+1]*8+Auto_Int_List[x-2]-Auto_Int_List[x-1] * 8-Auto_Int_List[x+2])/((Auto_RT_List[x+2]-Auto_RT_List[x-2])*3)
        elif x < 2:
            FirstDerivative = (Auto_Int_List[x]*(-3)+Auto_Int_List[x+1]*4+Auto_Int_List[x+2])/((Auto_RT_List[x+2]-Auto_RT_List[x])*3)
        elif x > len(Auto_Int_List)-3:
            FirstDerivative = (Auto_Int_List[x]*(3)-4*Auto_Int_List[x-1]-Auto_Int_List[x-2])/((Auto_RT_List[x]-Auto_RT_List[x-2])*3)
        return FirstDerivative
    def TFSD(x, Auto_Int_List, Auto_RT_List):
        if 2<= x <= len(Auto_Int_List)-3:
            SecondDerivative = (Auto_Int_List[x+1]*16+Auto_Int_List[x-1]*16-Auto_Int_List[x] * 30-Auto_Int_List[x+2]-Auto_Int_List[x-2])/((Auto_RT_List[x+2]-Auto_RT_List[x-2])*3)
        elif x < 2:
            SecondDerivative = (Auto_Int_List[x]-Auto_Int_List[x+1]+Auto_Int_List[x+2])/((Auto_RT_List[x+2]-Auto_RT_List[x])*3)
        elif x > len(Auto_Int_List)-3:
            SecondDerivative = (Auto_Int_List[x]-Auto_Int_List[x-1]+Auto_Int_List[x-2])/((Auto_RT_List[x]-Auto_RT_List[x-2])*3)
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
