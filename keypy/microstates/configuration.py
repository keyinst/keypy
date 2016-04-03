# -*- coding: utf-8 -*-

#################
# #Create Class Configuration Object
#################

class MstConfiguration(object):
    """
    Class that defines the parameters used for microstate computation and visualization.

    Attributes
    ----------
        subtract_column_mean_at_start (bool) : bool
			subtract column mean within microstate preprocessing
        debug : bool
			set to True for more detailed printing
        use_gfp_peaks : bool
			compute microstates based on global field power peaks only
        force_avgref : bool
			recompute average reference before microstate computation
        set_gfp_all_1 : bool
			set global field power to 1 before microstate computation
        use_smoothing : bool
			smoothe global field power curve before peak extraction (default = False)
        gfp_type_smoothing : {'hamming', 'hanning'}
            `hamming` : use hamming window to smooth global field power curve before peak identification
            `hanning` : use hanning window to smooth global field power curve before peak identification		
        smoothing_window : int
        		window for smoothing in time frames to reduce high frequency noise before peak identification (depends on time frames per second), e.g. 100.
        use_fancy_peaks : bool
            Whether a particular smoothing algorithm from scipy.signal.find_peaks_cwt is applied before peak computation or not.
            Reference: Bioinformatics (2006) 22 (17): 2059-2065. doi: 10.1093/bioinformatics/btl355 http://bioinformatics.oxfordjournals.org/content/22/17/2059.long)           
        method_GFPpeak : {'GFPL1', 'GFPL2'}
            `GFPL1` : use L1-Norm to compute GFP peaks (default)
            `GFPL2` : use L2-Norm to compute GFP peaks
        original_nr_of_maps : int
            Number of maps to compute from microstate algorithm (default = 4).
        seed_number : int
            Number of seeds used for microstate algorithms (default = 3000). It should be at least 3 times larger than the average number of GFP peaks you would like to compute the microstates based on.
        max_number_of_iterations : int
            The maximal number of iterations performed to increase the variance explained by the defined N microstates (default = 500). 
        ERP : bool
            Whether microstate computation is done based on ERP (time-locked) data (in this case map polarity is considered). (default = False)
        correspondance_cutoff : double
            Pearson correlation coefficient minimum or dissimilarity maximum that is regarded as the microstate of a particular class (used for microstate class visualization across time) (default = 0 which means no cutoff).
        fixed_seed : int
            Fixate seed for testing for microstate and modelmaps algorithms (default=100). Can be None!
        dissimilarity_measure : {'correlation', 'dissimilarity'}
            Determines the measure to assess similarity between two maps (default='dissimilarity').
        modelmaps_normalization_type : {'vector_norm_1','gfp1','none'}
            Specifies how microstates and modelmaps are normalized during microstate and modelmap computation (defualt='vector_norm_1').
    """

    def __init__(self, subtract_column_mean_at_start = False, debug = True, use_gfp_peaks = True, force_avgref = True, set_gfp_all_1 = False, use_smoothing = False, gfp_type_smoothing='hamming', smoothing_window=100, use_fancy_peaks = False, method_GFPpeak = 'GFPL1', original_nr_of_maps = 4, seed_number = 3000, max_number_of_iterations = 500, ERP = False, correspondance_cutoff = 0, fixed_seed = 100, similarity_measure = 'dissimilarity', modelmaps_normalization_type ='vector_norm_1'):
        """ This constructor initializes the class members with the corresponding parameter values. 
        """
        self.subtract_column_mean_at_start = subtract_column_mean_at_start
        self.debug = debug
        self.use_gfp_peaks = use_gfp_peaks
        self.force_avgref = force_avgref
        self.set_gfp_all_1 = set_gfp_all_1
        self.use_smoothing = use_smoothing
        self.gfp_type_smoothing = gfp_type_smoothing
        self.smoothing_window = smoothing_window
        self.use_fancy_peaks = use_fancy_peaks
        self.method_GFPpeak = method_GFPpeak
        self.original_nr_of_maps = original_nr_of_maps
        self.seed_number = seed_number
        self.max_number_of_iterations = max_number_of_iterations
        self.ERP = ERP
        self.correspondance_cutoff = correspondance_cutoff
        self.fixed_seed = fixed_seed
        self.similarity_measure = similarity_measure
        self.modelmaps_normalization_type = modelmaps_normalization_type

#---------------------------------