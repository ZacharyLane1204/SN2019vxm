
##################################################
# Parameters for building the observation to fit #
##################################################

# The observed dataset is expressed as 4 array-likes containing respectively:
# - bands: a list of strings corresponding to the unique identifiers also used
#          in the 'filters' variable below 
# - fluxes: measures of fluxes (or upper limits) at the corresponding bands
#           listed in variable 'bands'.
#           The input measure should be given in units of milli-Jansky
# - errors: 1-sigma error on the measures of the fluxes (or upper limits) listed
#           in the variable 'fluxes'
# - uplims: sequence of booleans identifying whether the values listed 
#           in argument ``fluxes`` should be considered non-detection (``True``)
#           or a detection (``False``)
bands  = ['SDSS.u', 'SDSS.g', 'SDSS.r', 'SDSS.i', 
          'SDSS.z', 'WISE.W1', 'WISE.W2']
fluxes = [0.01832719496888013, 0.09174618195405225, 
          0.1423206766853215, 0.16664028511455084, 
          0.10843962323779213, 0.4988844874600117, 
          0.43251383103500923]

errors = [0.004828557725035493, 0.0017838151826790193, 
          0.0028976281183520605, 0.0043839823246877565, 
          0.015185305742367743, 0.031245292363478026, 
          0.10158178978520571]
uplims = [False, False, False, False, False, False, False]

# The photometric filters system used in the observation.
# This parameter should be an iterable containing names
# of filters already present in the database, e.g.
#
# filters = ['GOODS.b', 'GOODS.i', 'GOODS.v', 'GOODS.z']
#
# NOTE that, if the bands listed in variable ``bands`` 
# are all present in the database and have the same names
# of the filters listed with ``galapy.PhotometricSystem.print_filters()``,
# this variable can be set as ``filters = bands``
filters = bands

# Eventual custom photometric filters. 
# This parameter should be a nested dictionary with user-defined transmissions.
# Such transmissions are passed as properly formatted dictionaries
# 
# filters_custom = { 'filter_1' : { 'wavelengths' : array-like,
#                                     'photons' : array-like },
#                     'filter_2' : { 'wavelengths' : array-like,
#                                     'photons' : array-like },
#                     ...
#                     'filter_N' : { 'wavelengths' : array-like,
#                                     'photons' : array-like } }
#
# As the keywords in the lower level dictionaries suggest,
# the two arrays provided, must define the
# wavelength grid and the corresponding transmission in photon units
# (and thus they must have the same size).
# Note that the chosen keyword in the higher level dictionary
# (e.g. 'filter_1') will be used as unique identifier of
# the custom transmission.
filters_custom = None

# Method for treatment of the upper limits, when present.
# Available methods are:
# - 'simple' : a heaviside function which is 0. when the model's flux is
#              smaller than the upper-limit and +infty otherwise
# - 'chi2' : the distance between model and data is expressed as a normal chi-squared
# - 'S12' : modified version of the chi-squared integrating a gaussian, with mean = data-flux
#           and std = data-error, up to the value of the model's flux (Sawicki, 2012)
method_uplims = 'chi2'

########################################
# Parameters defining the galaxy model #
########################################

# Cosmological model to use for computing distances and ages.
# Pre-computed cosmologies present in the database:
# - WMAP7
# - WMAP9
# - Planck15
# - Planck18
# To use another cosmological model the user has to provide instead
# a dictionary with the following key-value couples:
# - 'redshift' : an iterable containing monothonically increasing
#                values of redshift
# - 'luminosity_distance' : an iterable with the luminosity distances
#                           corresponding to the redshift values in
#                           the first iterable
# - 'age' : an iterable with the age of the universe corresponding to
#           the redshift values in the first iterable
cosmo     = 'Planck18'

# The Star-Formation History model to use for building the galaxy model.
# Available models are:
# - 'insitu'
# - 'constant'
# - 'delayedexp'
# - 'lognormal'
sfh_model = 'lognormal'

# The SSP library to use for building the unattenuated stellar emission component
ssp_lib   = 'parsec22.NT'

# Whether to provide X-ray emission support (True) or not (False). 
do_Xray = False

# Whether to provide radio-emission support (True) or not (False).
do_Radio = False

# Whether to build a galaxy containing an AGN (True) or not (False).
do_AGN = False

# Sub-sampling of the wavelength grid.
# If lstep is an integer it will consider a wavelength grid entry every lstep values.
# If lstep is a sequence of integers or a mask, only the wavelength grid entries
# corresponding to the indices provided will be considered.
# If None, it will consider the whole wavelength grid (safest choice)
lstep = None

# Eventual noise model to add. The default is ``None``, i.e. no noise will be added.
# Valid choices for the noise are:
# - 'calibration_error' : accounts for eventual systematics in the calibration of
#                         the errors in the observations.
#
# If this is set to ``None`` all the other hyper-parameters in this file
# related to noise will be ignored.
noise_model = None

# Eventual keyword arguments to be passed to the noise model of choice 
# (leave empty for no keyword arguments)
noise_kwargs = {}

#############################################
# Here define the fixed and free parameters #
#############################################

# The dictionary provided here will be used to set the fixed parameters to
# the given value or to flag parameters as 'free' with some associated prior.
# All the parameters that are not used in the specified galaxy model will be ignored.
# (e.g. if the galaxy model has been built with ``do_AGN = False`` all the eventual
#  AGN-related parameters provided will be ignored)
#
# - To set the parameter 'fixed_parameter' as FIXED provide a single value:
#   parameters = {
#       ...
#       'fixed_parameter' : some_float_value,
#       ...
#   }
#
# - To set the parameter 'free_parameter' as FREE provide a tuple
#   parameters = {
#       ...
#       'free_parameter' : ( a_list, a_bool ),
#       ...
#   }
#   The list ``a_list`` contains the minimum and maximum value of the
#   UNIFORM prior from which to draw samples for the 'free_parameter'.
#   The boolean ``a_bool`` states whether the prior has to be considered
#   * logarithmic: ``a_bool = True``, therefore samples will be drawn from the interval
#                  10**min(a_list) < 'free_parameter' < 10**max(a_list)
#   * linear: ``a_bool = False``, therefore samples will be drawn from the interval
#             min(a_list) < 'free_parameter' < max(a_list)

galaxy_parameters = {
    
    ##########
    # Galaxy #
    ##########

    'age'      : ( [3., 11.], True ),
    'redshift' : 0.019,

    ##########################
    # Star Formation History #
    ##########################
    
    'sfh.tau_quench' : ( [3., 11.], True ),

    # Log-Normal model
    'sfh.psi_norm'   : ( [0., 4.], True ),
    'sfh.sigma_star' : ( [0., 5.], False ),
    'sfh.tau_star'   : ( [6., 11.], True ),
    'sfh.Mdust' : ( [6., 14.], True ),
    'sfh.Zgxy'  : ( [0., 1.], False ),
        
    ########################
    # Inter-Stellar Medium #
    ########################

    'ism.f_MC' : ( [0., 1.], False ),

    ##################
    # Molecular Clouds
    
    'ism.norm_MC' : 100.,
    'ism.N_MC'    : ( [0., 5.], True ),
    'ism.R_MC'    : ( [0., 5.], True ),
    'ism.tau_esc' : ( [4., 8.], True ),
    'ism.dMClow'  : 1.3,
    'ism.dMCupp'  : 1.6,
    
    ##############
    # Diffuse Dust
    
    'ism.norm_DD' : 1.0,
    'ism.Rdust'   : ( [0., 5.], True ),
    'ism.f_PAH'   : ( [0., 1.], False ),
    'ism.dDDlow'  : 0.7,
    'ism.dDDupp'  : 2.0,

    ###############
    # Synchrotron #
    ###############

    'syn.alpha_syn'   : 0.75,
    'syn.nu_self_syn' : 0.2,

    #####################
    # Nebular Free-Free #
    #####################

    'nff.Zgas' : 0.02,
    'nff.Zi'   : 1.,

    ###########################
    # Active Galactic Nucleus #
    ###########################

    'agn.fAGN' : ( [-3., 3.], True ),
    
    'agn.ct' : 40,
    'agn.al' : 0.,
    'agn.be' : -0.5,
    'agn.ta' : 6.,
    'agn.rm' : 60,
    'agn.ia' : 0.001,
    
}

noise_parameters = {

    #########
    # Noise #
    #########

    ###################
    # Calibration Error
    
    'f_cal' : ( [-10., 1.], True ),
}

##############################
# Parameters for the sampler #
##############################

# Choose the sampler, valid options are
# 'emcee' : Affine Invariant MCMC ensemble sampler
# 'dynesty' : Dynamic Nested Sampler
sampler = 'dynesty'

# EMCEE SAMPLER-SPECIFIC MANDATORY PARAMETERS
# - set the number of walkers (``nwalkers``)
# - set the chain length (``nsamples``)
nwalkers = 2000
nsamples = 10000

# Sampler keyword arguments.
# These are the parameters that will be passed to the constructor of
# the chosen sampler. See relative documentation for further informations:
# - emcee :
# - dynesty : e.g. to decrease the number of walkers
# sampler_kw = {'walks':25}  # default is 'walks' : 50
# ( 'walks' : >= 50 is recommended for 15 < ndim < 25 )
sampler_kw = {'walks':120}

# Sampling keyword arguments
# These are the parameters that will be passed to the routine running the
# parameter-space sampling. See relative documentation for further informations:
# - emcee :
# - dynesty :
sampling_kw = {'dlogz_init': 1e-3}

# Output directory (note that if the directory does not exist it will be created)
output_directory = ''

# An identification name, it will be pre-pended to all files stored in the output directory
# (if the string is empty the current date+time will be used)
run_id = ''

# The method used for storing results. 
# Possible choices are:
# - 'hdf5' : uses HDF5 format to save a dictionary containing all the informations
#            to build the Results class used for the analysis and visualization of the
#            sampling-run results. (This is the safest choice, also in terms of security)
# - 'pickle' : uses the standard python pickle format to directly save an instance of
#              the Results class used for the analysis and visualization of the
#              sampling-run results. (Recommended only for local usage, pickled objects
#              are not safe for distribution)
store_method = 'pickle'

# Only available if the output format chosen is HDF5 (see parameter store_method).
# If True it stores only the chains, weights and loglikelihood values obtained by the
# sampling run (along with information to re-build all the models used in running.
# If False it will save all the infos on the derived quantities as well.
# Selecting lightweight storage, it does not change the way users will ultimatelly access
# the Results object, what will change is the time spent for loading it.
# With lightweight storage the size of the output file is smaller, and the Results object
# is computed (instantiated) when loading the file (so this will take up to some minutes).
# With lightweight storage off, all the quantities are computed at the end of the sampling
# run and are ready to use but the output file can reach a size of up to some GiB. 
store_lightweight = False

# Whether to pickle the sampler raw results.
# (might be useful for analyzing run statistics)
pickle_raw = False

# Whether to pickle the sampler at the end-of-run state.
# (might be useful for extending the run)
pickle_sampler = False

#############################
# ... and that's all folks! #
#############################
