# -*- coding: utf-8 -*-

"""
Functions and classes that allow the computation of the EEG microstates based on preprocessed (no artefacts, average-referenced, filtered) EEG data epochs.

configuration : classes
    Defines the parameters used for microstate computation and visualization.
data_provider : classes
    Contains classes used to select the correct input and write the correct output for the modelmaps, sortmaps, and parameters scripts.
microstates : functions
    Compute EEG microstates for each dataset in inputhdf5 of name microstate_input.
microstates_helper : functions
    Contains helper functions for microstate computation.
modelmaps : functions
    Computes four "mean" (principal component) modelmaps for the specified input microstate maps.
sortmaps : functions
    Sorts the modelmaps based on user-defined modelmaps or the modelmaps of a particular level.		
"""

#############################################################################

__all__ = ["configuration", "data_provider", "microstates", "microstates_helper", "modelmaps", "sormaps", "parameters"]