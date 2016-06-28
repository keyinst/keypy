# -*- coding: utf-8 -*-

"""
Functions and classes that allow the computation of the EEG microstates based on preprocessed (no artefacts, average-referenced, filtered) EEG data epochs.

configuration : classes
    Defines the parameters used for microstate computation and visualization.
modelmaps : functions
    Compute EEG microstate models for each dataset in inputhdf5 of name modelmaps_input.
micorstates_helper : functions
    Contains helper functions for modelmap, meanmod, sortmaps, and parameter computation.
meanmods : functions
    Computes four "mean" (principal component) modelmaps (meanmods) for the specified input modelmaps.
meanmods_provider : classes
    Contains classes used to select the correct input and write the correct output for the meanmods script.
sortmaps : functions
    Sorts the mean models / modelmaps based on user-defined modelmaps or the mean models of a particular level.		
sortmaps_provider : classes
    Contains classes used to select the correct input and write the correct output for the sortmaps script.
parameters : functions
    Computes parameters (coverage, duration, occurrence) of spontaneous EEG based on user-defined modelmaps or the mean models of a particular level.
parameters_provider : classes
    Contains classes used to select the correct input and write the correct output for the parameters script.    	
"""

#############################################################################

__all__ = ["configuration", "modelmaps", "microstates_helper", "meanmods", "meanmods_provider", "sortmaps", "sortmaps_provider", "parameters", "parameters_provider"]