# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 07:47:48 2019

@author: hcji
"""

import os
import json
import joblib
import numpy as np
import pandas as pd
import tensorflow.keras as keras
import matplotlib.pyplot as plt
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.layers import Dense, Input, Flatten, Conv1D, MaxPooling1D, concatenate, Embedding, LSTM
from tensorflow.keras import metrics, optimizers
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, mean_absolute_error, r2_score
from tqdm import tqdm
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestClassifier
from smiles_to_onehot.encoding import get_dict, one_hot_coding


def build_RI_model_descriptor(morgan, cdkdes, RI, descriptor, save_name):
    # remove nan
    keep = np.where(~ np.isnan(RI))[0]
    if descriptor == 'all':
        X = np.hstack((morgan[keep,:], cdkdes[keep,:]))
    elif descriptor == 'morgan':
        X = morgan[keep,:]
    else:
        X = cdkdes[keep,:]
    Y = RI[keep]
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.1)
    
    # scale
    scaler = MinMaxScaler()
    scaler.fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)
    joblib.dump(scaler, 'Model/RI/' + save_name + '_scaler.save') 
    
    # train model
    layer_in = Input(shape=(X.shape[1],), name="morgan_fp")
    layer_dense = layer_in
    n_nodes = X.shape[1]
    for j in range(5):
        layer_dense = Dense(int(n_nodes), activation="relu")(layer_dense)
        n_nodes = n_nodes * 0.5
    layer_output = Dense(1, activation="linear", name="output")(layer_dense)
    model = Model(layer_in, outputs = layer_output) 
    opt = optimizers.Adam(lr=0.001)
    model.compile(optimizer=opt, loss='mse', metrics=[metrics.mae])
    history = model.fit(X_train, Y_train, epochs=20, batch_size=1024, validation_split=0.11)
    '''
    # plot loss
    plt.plot(history.history['loss'])
    plt.plot(history.history['mean_absolute_error'])
    plt.plot(history.history['val_loss'])
    plt.plot(history.history['val_mean_absolute_error'])
    plt.ylabel('values')
    plt.xlabel('epoch')
    plt.legend(['loss', 'mae', 'val_loss', 'val_mae'], loc='upper left')
    # plt.show()
    plt.savefig("Result/retention_" + save_name + '_loss.png')
    '''
    # predict
    Y_predict = model.predict(X_test)
    r2 = round(r2_score(Y_predict, Y_test), 4)
    mae = round(mean_absolute_error(Y_predict, Y_test), 4)
    plt.cla()
    plt.plot(Y_test, Y_predict, '.', color = 'blue')
    plt.plot([0,4500], [0,4500], color ='red')
    plt.ylabel('Predicted RI')
    plt.xlabel('Experimental RI')        
    plt.text(0, 4000, 'R2='+str(r2), fontsize=12)
    plt.text(0, 3600, 'MAE='+str(mae), fontsize=12)
    # plt.show()
    plt.savefig("Result/retention_" + save_name + '_r2.png')
    plt.close('all')
    
    # save model
    model.save('Model/RI/' + save_name + '_model.h5')
    return {'r2': r2, 'mae': mae, 'model': model}


def build_RI_model_CNN(smiles, RI, method, save_name):
    words = get_dict(smiles, save_path='Model/RI/' + save_name + '_dict.json')
    keep = np.where(~ np.isnan(RI))[0]
    X = []
    Y = []
    for i, smi in enumerate(smiles):
        if i not in keep:
            continue
        xi = one_hot_coding(smi, words, max_len=600)
        if xi is not None:
            X.append(xi.todense())
            Y.append(RI[i])
    X = np.array(X)
    Y = np.array(Y)
    
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.1)
    layer_in = Input(shape=(X.shape[1:3]), name="smile")
    layer_conv = layer_in
    if method == 'single_channel':
        for i in range(5):
            layer_conv = Conv1D(128, kernel_size=4, activation='relu', kernel_initializer='normal')(layer_conv)
            layer_conv = MaxPooling1D(pool_size=2)(layer_conv)
        layer_dense = Flatten()(layer_conv)
    else:
        layer_conv_k1 = Conv1D(128, kernel_size=3, activation='relu', kernel_initializer='normal')(layer_in)
        layer_conv_k1 = MaxPooling1D(pool_size=2)(layer_conv_k1)
        layer_conv_k2 = Conv1D(128, kernel_size=4, activation='relu', kernel_initializer='normal')(layer_in)
        layer_conv_k2 = MaxPooling1D(pool_size=2)(layer_conv_k2)
        layer_conv_k3 = Conv1D(128, kernel_size=5, activation='relu', kernel_initializer='normal')(layer_in) 
        layer_conv_k3 = MaxPooling1D(pool_size=2)(layer_conv_k3)
        layer_conv_k1 = Conv1D(128, kernel_size=3, activation='relu', kernel_initializer='normal')(layer_conv_k1)
        layer_conv_k1 = MaxPooling1D(pool_size=2)(layer_conv_k1)
        layer_conv_k2 = Conv1D(128, kernel_size=4, activation='relu', kernel_initializer='normal')(layer_conv_k2)
        layer_conv_k2 = MaxPooling1D(pool_size=2)(layer_conv_k2)
        layer_conv_k3 = Conv1D(128, kernel_size=5, activation='relu', kernel_initializer='normal')(layer_conv_k3)
        layer_conv_k3 = MaxPooling1D(pool_size=2)(layer_conv_k3)
        layer_dense_k1 = Flatten()(layer_conv_k1)
        layer_dense_k2 = Flatten()(layer_conv_k2)
        layer_dense_k3 = Flatten()(layer_conv_k3)
        layer_dense = concatenate([layer_dense_k1, layer_dense_k2, layer_dense_k3], axis=-1)
        
    for i in range(1):
        layer_dense = Dense(32, activation="relu", kernel_initializer='normal')(layer_dense)
    layer_output = Dense(1, activation="linear", name="output")(layer_dense)
    model = Model(layer_in, outputs = layer_output) 
    opt = optimizers.Adam(lr=0.001)
    model.compile(optimizer=opt, loss='mse', metrics=[metrics.mae])
    history = model.fit(X_train, Y_train, epochs=20, validation_split=0.11)
    '''
    # plot loss
    plt.plot(history.history['loss'])
    plt.plot(history.history['mean_absolute_error'])
    plt.plot(history.history['val_loss'])
    plt.plot(history.history['val_mean_absolute_error'])
    plt.ylabel('values')
    plt.xlabel('epoch')
    plt.legend(['loss', 'mae', 'val_loss', 'val_mae'], loc='upper left')
    # plt.show()
    plt.savefig("Result/retention_" + save_name + '_loss.png')
    '''
    # predict
    Y_predict = model.predict(X_test)
    r2 = round(r2_score(Y_predict, Y_test), 4)
    mae = round(mean_absolute_error(Y_predict, Y_test), 4)
    plt.cla()
    plt.plot(Y_test, Y_predict, '.', color = 'blue')
    plt.plot([0,4500], [0,4500], color ='red')
    plt.ylabel('Predicted RI')
    plt.xlabel('Experimental RI')        
    plt.text(0, 4000, 'R2='+str(r2), fontsize=12)
    plt.text(0, 3600, 'MAE='+str(mae), fontsize=12)
    # plt.show()
    plt.savefig("Result/retention_" + save_name + '_r2.png')
    plt.close('all')
    
    # save model
    model.save('Model/RI/' + save_name + '_model.h5')
    return {'r2': r2, 'mae': mae, 'model': model}        
    

def build_RI_model_combine(morgan, cdkdes, smiles, RI, save_name):
    words = get_dict(smiles, save_path='Model/RI/' + save_name + '_dict.json')
    keep = np.where(~ np.isnan(RI))[0]
    X1 = []
    X2 = []
    Y = []
    for i, smi in enumerate(tqdm(smiles)):
        if i not in keep:
            continue
        xi = one_hot_coding(smi, words, max_len=1000)
        if xi is not None:
            xi = xi.todense()
            xj = []
            for k in xi:
                if np.sum(k) > 0:
                    xj.append(k.argmax())
                else:
                    break
            x2 = np.hstack((morgan[i], cdkdes[i]))
            X1.append(np.array(xj))
            X2.append(x2)
            Y.append(RI[i])
    


def build_RI_model_RNN(smiles, RI, save_name):
    words = get_dict(smiles, save_path='Model/RI/' + save_name + '_dict.json')
    keep = np.where(~ np.isnan(RI))[0]
    X = []
    Y = []
    for i, smi in enumerate(tqdm(smiles)):
        if i not in keep:
            continue
        xi = one_hot_coding(smi, words, max_len=1000)
        if xi is not None:
            xi = xi.todense()
            xj = []
            for k in xi:
                if np.sum(k) > 0:
                    xj.append(k.argmax())
                else:
                    break
            X.append(np.array(xj))
            Y.append(RI[i])
    # X = np.array(X)
    # Y = np.array(Y)    

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.1)
    
    X_train = pad_sequences(X_train, maxlen=1000)
    X_test = pad_sequences(X_test, maxlen=1000)
    
    model = Sequential()
    model.add(Embedding(X_train.shape[0], 300))
    model.add(LSTM(256))
    model.add(Dense(1, activation='linear'))    
    opt = optimizers.Adam(lr=0.001)
    model.compile(optimizer=opt, loss='mse', metrics=[metrics.mae])
    history = model.fit(X_train, Y_train, epochs=10, validation_split=0.11)    

    # plot loss
    plt.plot(history.history['loss'])
    plt.plot(history.history['mean_absolute_error'])
    plt.plot(history.history['val_loss'])
    plt.plot(history.history['val_mean_absolute_error'])
    plt.ylabel('values')
    plt.xlabel('epoch')
    plt.legend(['loss', 'mae', 'val_loss', 'val_mae'], loc='upper left')
    plt.show()
    
    # predict
    Y_predict = model.predict(X_test)
    r2 = round(r2_score(Y_predict, Y_test), 4)
    mae = round(mean_absolute_error(Y_predict, Y_test), 4)
    plt.cla()
    plt.plot(Y_test, Y_predict, '.', color = 'blue')
    plt.plot([0,4500], [0,4500], color ='red')
    plt.ylabel('Predicted RI')
    plt.xlabel('Experimental RI')        
    plt.text(0, 4000, 'R2='+str(r2), fontsize=15)
    plt.show()
    
    # save model
    model.save('Model/RI/' + save_name + '_model.h5')
    return {'r2': r2, 'mae': mae, 'model': model}      
