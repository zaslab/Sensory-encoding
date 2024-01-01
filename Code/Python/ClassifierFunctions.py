# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:09:38 2023

@author: ZasLab
"""

import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import sys
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
def merge_conditions(conditions, groups):
    """
    conditions: a 1-dimensional array of the condition codes (one entry per row in neuron data)
    groups: a list of tuples, each tuple in format `(group_name, (tuple of condition names))`
    """
    grouped_conditions = -1 * np.ones_like(conditions) # We define -1 as an invalid condition
    group_names = []
    for i, (group_name, group) in enumerate(groups):
        group_names.append(group_name)
        for condition in group:
            grouped_conditions[conditions == condition2num[condition]] = i
            
    return grouped_conditions, group_names

def filter_conditions(x, y):
    x = x[y >= 0]
    y = y[y >= 0]
    return x, y

def prepare_data(neuron_names, conditions):
    x = np.concatenate([wt_neurons_data[neuron_name] for neuron_name in neuron_names], axis=1)
    x, y = filter_conditions(x, conditions)
    return x, y

class SensoryData(object):
    def __init__(self, name, neuron_data, conditions):
        super(SensoryData, self).__init__()
        self.name = name
        self.neuron_data = neuron_data
        
        self.all_conditions, self.all_condition_names = merge_conditions(conditions, [
            (condition_name, (condition_name,)) for condition_name in ['IAA', 'DA', 'Gly', 'NaCl', 'Quin', 'SDS', 'rho',]
        ])
        self.valence_conditions, self.valence_condition_names = merge_conditions(conditions, [
            ('Positive', ('IAA', 'DA', 'NaCl')), # Positive valence
            ('Negative', ('Quin', 'SDS', 'Gly')), # Negative valence
        ])
        self.volatile_conditions, self.volatile_condition_names = merge_conditions(conditions, [
            ('Volatile', ('IAA', 'DA')), # Volatile valence
            ('Soluble', ('Gly', 'NaCl', 'Quin', 'SDS')), # Soluble valence
        ])
        self.conditions = {
            'all': (self.all_conditions, self.all_condition_names),
            'valence': (self.valence_conditions, self.valence_condition_names),
            'volatility': (self.volatile_conditions, self.volatile_condition_names),
        }
    
    def get_data(self, neuron_names, conditions,step):
        self._check_supported_conditions(conditions)
        
        x = np.concatenate([self.neuron_data[neuron_name] for neuron_name in neuron_names], axis=1)
        if step<3:
            x = x[:, list(range(int(step - 1),26,2))]
            
        x, y = filter_conditions(x, self.conditions[conditions][0])
        
        return x, y
    
    def get_labels(self, conditions):
        self._check_supported_conditions(conditions)
        
        return self.conditions[conditions][1]

    def _check_supported_conditions(self, conditions_name):
        if conditions_name not in self.conditions.keys():
            raise Exception(f'Supported conditions are: {list(self.conditions.keys())}')
class Test(object):
    def __init__(self, classifier_cls, classifier_params):
        super(Test, self).__init__()
        self.classifier_cls = classifier_cls
        self.classifier_params = classifier_params
        self.classifier = None
        self.y_values = []
        self.y_predictions = []
        
    @property
    def confusion_matrix(self):
        return confusion_matrix(
            self.y_values, self.y_predictions,
            sample_weight=np.ones(len(self.y_values), dtype=np.float64), # This parameter is a hack to get the confusion matrix as floats
            normalize='true'
        )
    
    def _create_classifier(self):
        self.classifier = self.classifier_cls(**self.classifier_params)
        
    def reset(self):
        # Delete the classifier object so that a new one will be created with a clean state
        del self.classifier
        self.classifier = None
    
    def train(self, x_train, y_train):
        if self.classifier is None:
            # Create a new object, this is to support the case that the classifier updates an internal state
            # instead of leaning new from scratch on the data provided
            self._create_classifier()
        self.classifier.fit(x_train, y_train)
    
    def test(self, x_test, y_test):
        if self.classifier is None:
            raise Exception('Classifier must be trained before testing!')
        
        # Keep real Y values to build confusion matrix
        self.y_values.extend(y_test)

        predictions = self.classifier.predict(x_test)
        self.y_predictions.extend(predictions)
        
    def train_and_test(self, x, y, should_reset=True):
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3, stratify=y)
        
        if should_reset:
            self.reset()
        
        self.train(x_train, y_train)
        self.test(x_test, y_test)
        
    def plot_confusion_matrix(self, condition_names, ax=None):
        fig = None
        if ax is None:
            fig, ax = plt.subplots()
            
        
        sns.heatmap(self.confusion_matrix,
                    annot=True,
                    vmin=0, vmax=1,
                    xticklabels=condition_names, yticklabels=condition_names,
                    ax=ax
                   )
        
        if fig is not None:
            fig.show()
        
