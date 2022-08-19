'''
Description:
These are general functions to be used in the river node, raster, 
and lake processors for aggregating heights and areas with their 
corresponding uncertainties


Copyright (c) 2018-, California Institute of Technology ("Caltech"). U.S.
Government sponsorship acknowledged.
All rights reserved.

Author (s): Brent Williams
'''
import sys
import numpy as np
import scipy.stats

from scipy import interpolate

from jpl.modules.constants import AGG_CLASSES
FLOAT_EPS = sys.float_info.epsilon

def simple(in_var, metric='mean', pcnt=68):
    """
    Aggregate the input variable according to desired metric/accumulator.
    
    INPUT:
    in_var = 1d list or array of a given variable for all pixels over a feature 
             (river node, raster bin, lake)
    """
    if metric == 'mean':
        out_var = np.nanmean(in_var)
    elif metric == 'median':
        out_var = np.nanmedian(in_var)
    elif metric == 'sum':
        out_var = np.nansum(in_var)
    elif metric == 'std':
        out_var = np.nanstd(in_var)
    elif metric == 'pcnt':
        out_var = np.percentile(in_var, pcnt)
    elif metric == 'count':
        out_var = np.sum(np.ones(np.shape(in_var)))
    elif metric == 'mode':
        out_var,_ = scipy.stats.mode(in_var)
    return out_var

def sig0_with_uncerts(
        sig0, good, sig0_std,
        num_rare_looks=1.0, num_med_looks=1.0,  method='rare'):
    """
    Return the aggregate sig0, sample std, and estimate of uncert
    implements methods: rare (default), medium which is the assumption
    of the kind of sig0 input from the pixel cloud.

    INPUTS:
    sig0        = 1d array of rare (or medium) sig0s over the water feature
    good        = mask for filtering out some pixels
    sig0_std    = pixel-wise sig0 error std/uncertainty from pixel cloud
    method      = type of sig0 input ('rare', or 'medium')
    num_*_looks = rare or medium number of looks (only used if method=medium)

    OUTPUTS:
    sig0_agg     = aggregated sig0
    sig0_std_out = sample std of sig0
    sig0_uncert  = estimate of sig0_agg 1-sigma uncertainty
    """
    # just do a simple average of sig0, nothing fancy
    sig0_agg = simple(sig0[good], metric='mean')

    sig0_std_out = None
    sig_uncert = None
    if method == 'rare':
        # rare method assumes sig0 is the rare sig0 and all pixels are indep
        # compute uncertainty as the std assuming
        # all sig0 measurements are independent
        num_pixels = simple(sig0[good], metric='count')
        sig_uncert = np.sqrt(simple(
            sig0_std[good]**2, metric='sum')) / num_pixels
        sig_std_out = simple(sig0[good], metric='std') 
    elif method == 'medium':
        # assumes that sig0 is the medium sig0 and computes the
        # uncert using the sample std scaled by rare and medium looks
        sig_uncert = None # may want to implement this if we want to use med sig0
        sig_std_out = height_uncert_std(
            sig0, good, num_rare_looks, num_med_looks)
    return sig0_agg, sig_std_out, sig_uncert

def height_only(height, good, height_std=1.0, method='weight'):
    """
    Return the aggregate height
    implements methods: weight (default), median, uniform 
    good is a mask used to filter out heights that are expected to be bad 
    or outliers, if desired
    
    INPUTS:
    height      = 1d array of heights over the water feature
    good        = mask for filtering out some pixels
    height_std  = pixel-wise height error std (phase_std * dh_dphase)
    method      = type of aggregator ('weight', 'uniform', or 'median')

    OUTPUTS:
    height_out  = aggregated height
    weight_norm = normalized weighting array

    Reference: implements Eq. (1) in 
    "SWOT Hydrology Height and Area Uncertainty Estimation," 
    Brent Williams, 2018, JPL Memo
    """
    if method == 'median':
        # for median, use uniform weighting for uncertainty estimation later
        height_agg = simple(height[good], metric='median')
        num_pixels = simple(height[good], metric='count')
        weight = np.ones(np.shape(height))
        return height_agg, weight/num_pixels
    elif method == 'uniform':
        weight = np.ones(np.shape(height)) 
    elif method == 'weight':
        # inverse height variance weighting
        weight = np.ones(np.shape(height))/(height_std)**2
    else:
        raise Exception("Unknown height aggregation method: {}".format(method))
    height_agg = simple(weight[good]*height[good], metric='sum')
    weight_sum = simple(weight[good], metric='sum')
    num_pixels = simple(height[good], metric='count')

    weight_sum_pixc = np.ones(np.shape(weight))
    weight_sum_pixc[good] = weight_sum

    height_out = height_agg/weight_sum
    weight_norm = weight/weight_sum_pixc
    return height_out, weight_norm

def height_uncert_std(
        height, good, num_rare_looks, num_med_looks, height_std=1.0,
        method='weight'):
    """
    Compute the sample standard deviation of the heights and scale by the
    appropriate factor instead of 1/sqrt(N), since the medium pixels are
    correlated

    INPUTS:
    height         = 1d array of heights over the water feature
    good           = mask for filtering out some pixels
    num_rare_looks = rare number of looks (either actual looks, or effective)
    num_med_looks  = medium number of looks (either actual looks, or effective)

    Reference: implements Eq. (14) in 
    "SWOT Hydrology Height and Area Uncertainty Estimation," 
    Brent Williams, 2018, JPL Memo
    """
    # need to do a weighted sample std when aggregating with weights
    # TODO: for median, should probably throw out outliers...
    weight = np.ones(np.shape(height))  # default to uniform
    if method == 'weight':
        weight = np.ones(np.shape(height)) / height_std ** 2
    height_agg = simple(weight[good]*height[good], metric='sum')
    weight_sum = simple(weight[good], metric='sum')
    height_mean = height_agg/weight_sum
    height_agg2 = simple(
        weight[good]*(height[good]-height_mean)**2.0, metric='sum')
    h_std = np.sqrt(height_agg2/weight_sum)

    num_pixels = simple(height[good], metric='count')
    # num_med_looks is rare_looks*num_pix_in_adaptive_window,
    # so need to normalize out rare to get number of independent pixels
    num_ind_pixels = simple(num_med_looks[good]/num_rare_looks[good],'mean')
    height_std_out = h_std * np.sqrt(num_ind_pixels/num_pixels)

    return height_std_out

def height_uncert_multilook(
        ifgram, power1, power2, weight_norm, good, num_rare_looks,
        looks_to_efflooks, dh_dphi, dlat_dphi, dlon_dphi):
    """
    compute height uncertainty bound by multilooking 
    everything over feature to compute average coh
    then computing phase noise std using CRB,
    then projecting through sensitivities

    INPUTS (all arrays are 1d lists for a given feature) 
    ifgram            = rare complex flattened interferogram
    power(1,2)        = the rare two channel interferogram powers
    good              = mask for filtering out some pixels if desired
    weight_norm       = the normalized weighting function
    num_rare_looks    = rare looks array
    looks_to_efflooks = scale factor to get effective looks from looks taken
    dh_dphi           = height sensitivity to phase array

    OUTPUTS:
    height_uncert_out = scalar height uncertainty for this feature using this
                        method

    Reference: implements square root of Eq. (7) in
    "SWOT Hydrology Height and Area Uncertainty Estimation,"
    Brent Williams, 2018, JPL Memo
    """
    # multilook the rare interferogram over the raster bin
    #  by averaging certain fields
    agg_real = simple(np.real(ifgram[good])*weight_norm[good])
    agg_imag = simple(np.imag(ifgram[good])*weight_norm[good])
    agg_p1 = simple(power1[good]*weight_norm[good])
    agg_p2 = simple(power2[good]*weight_norm[good])
    num_pixels = simple(power1[good],'count')

    # compute coherence
    coh = abs(agg_real + 1j *agg_imag)/np.sqrt(agg_p1*agg_p2)

    # get total num_eff_looks
    rare_looks = num_rare_looks#/looks_to_efflooks
    agg_looks = simple(rare_looks[good])

    num_looks = agg_looks * num_pixels

    # get phase noise variance using CRB
    phase_var = (0.5 / num_looks) * (1.0-coh**2)/(coh**2)
    # TODO: Probably want to set a more reasonable phase variance
    if phase_var < FLOAT_EPS: phase_var = FLOAT_EPS

    agg_dh_dphi = simple(dh_dphi[good]*weight_norm[good])
    agg_dh_dphi2 = simple(dh_dphi[good]**2*weight_norm[good])

    agg_dlat_dphi = simple(dlat_dphi[good]*weight_norm[good])
    agg_dlat_dphi2 = simple(dlat_dphi[good]**2*weight_norm[good])

    agg_dlon_dphi = simple(dlon_dphi[good]*weight_norm[good])
    agg_dlon_dphi2 = simple(dlon_dphi[good]**2*weight_norm[good])

    height_uncert_out = np.sqrt(
        phase_var) * np.abs(np.array(agg_dh_dphi2)/np.array(agg_dh_dphi))

    lat_uncert_out = np.sqrt(
        phase_var) * np.abs(np.array(agg_dlat_dphi2)/np.array(agg_dlat_dphi))

    lon_uncert_out = np.sqrt(
        phase_var) * np.abs(np.array(agg_dlon_dphi2)/np.array(agg_dlon_dphi))

    return height_uncert_out, lat_uncert_out, lon_uncert_out

def height_with_uncerts(
        height,  good, num_rare_looks, num_med_looks,
        ifgram, power1, power2, look_to_efflooks, dh_dphi,
        dlat_dphi, dlon_dphi, height_std=1.0, method='weight'):
    """
    Return the aggregate height with corresponding uncertainty
    implements methods: weight(default), median, uniform
    good is a mask used to filter out heights that are expected
    to be bad or outliers, if desired
    """
    # first aggregate the heights
    height_out, weight_norm = height_only(
        height,  good, height_std=height_std, method=method)

    # now compute uncertainties
    height_std_out = height_uncert_std(
        height, good, num_rare_looks, num_med_looks,
        height_std=height_std, method=method)

    height_uncert_out, lat_uncert_out, lon_uncert_out = height_uncert_multilook(
        ifgram, power1, power2, weight_norm, good,
        num_rare_looks, look_to_efflooks, dh_dphi, dlat_dphi, dlon_dphi)
    return (height_out, height_std_out, height_uncert_out, lat_uncert_out,
            lon_uncert_out)

def area_only(
        pixel_area, water_fraction, klass, good,
        interior_water_klasses=AGG_CLASSES['interior_water_klasses'],
        water_edge_klasses=AGG_CLASSES['water_edge_klasses'],
        land_edge_klasses=AGG_CLASSES['land_edge_klasses'],
        dark_water_klasses=AGG_CLASSES['dark_water_klasses'],
        method='composite'):
    """
    Return the aggregate area
    implements methods: simple, water_fraction, composite (default)
    good is a mask used to filter out pixels that are expected to be bad
    or outliers, if desired

    INPUTS:
    pixel_area     = 1d array of pixel_area
    water_fraction = 1d array of water fraction
    klass          = classification, with edges
    good           = mask for filtering out some pixels
    method         = type of aggregator('simple', 'water_fraction', 'composite')

    VARIABLES:
    Idw    = Binary detected-water mask representing pixels detected as water
    Idw_in = Binary detected-water mask for pixels detected as interior water
    Ide    = Binary detected-water mask for pixels detected as edge water
    I      = Binary mask for pixels detected as water or land-near-water

    OUTPUTS:
    area_out  = aggregated area

    Reference: implements Eq.s (15), (16), and (17) in 
    "SWOT Hydrology Height and Area Uncertainty Estimation," 
    Brent Williams, 2018, JPL Memo
    
    Updated to handle multiple class mappings for interior, edges, and dark
    """
    Idw_in = np.zeros(np.shape(pixel_area))
    Idw = np.zeros(np.shape(pixel_area))
    Ide = np.zeros(np.shape(pixel_area))

    for interior_water_klass in interior_water_klasses:
        # these should include all pixels aggregated as entirely water
        # including: dark water, low-coherence water etc...
        Idw_in[klass == interior_water_klass] = 1.0
        Idw[klass == interior_water_klass] = 1.0
    for water_edge_klass in water_edge_klasses:
        Idw[klass == water_edge_klass] = 1.0
        Ide[klass == water_edge_klass] = 1.0
    for land_edge_klass in land_edge_klasses:
        Ide[klass == land_edge_klass] = 1.0

    # handle current and legacy dark water classes like interior water
    for dark_water_klass in dark_water_klasses:
        Idw_in[klass == dark_water_klass] = 1.0
        Idw[klass == dark_water_klass] = 1.0

    I = np.zeros(np.shape(pixel_area))
    I[(Idw + Idw_in + Ide) > 0] = 1.0  # all pixels near water

    if method == 'simple':
        area_agg = simple(pixel_area[good] * Idw[good], metric='sum')
        num_pixels = simple(Idw[good], metric='sum')
    elif method == 'water_fraction':
        area_agg = simple(
            pixel_area[good] * water_fraction[good] * I[good], metric='sum')
        num_pixels = simple(I[good], metric='sum')
    elif method == 'composite':
        area_agg_in = simple(pixel_area[good] * Idw_in[good], metric='sum')
        area_agg_edge = simple(
            pixel_area[good] * water_fraction[good] * Ide[good], metric='sum')
        area_agg = area_agg_in + area_agg_edge
        num_pixels = simple(Idw_in[good] + Ide[good], metric='sum')
    else:
        raise Exception("Unknown area aggregation method: {}".format(method))
    return area_agg, num_pixels

def area_uncert(
        pixel_area, water_fraction, water_fraction_uncert, darea_dheight,
        klass, Pfd, Pmd, good, Pca=0.9, Pw=0.5,Ptf=0.5, ref_dem_std=10,
        interior_water_klasses=AGG_CLASSES['interior_water_klasses'],
        water_edge_klasses=AGG_CLASSES['water_edge_klasses'],
        land_edge_klasses=AGG_CLASSES['land_edge_klasses'],
        dark_water_klasses=AGG_CLASSES['dark_water_klasses'],
        method='composite'):
    '''
    Ie  = mask for edge pixels
    Pfd = Probability of false detection of water [0,1]
    Pmd = Probability of missed detection of water [0,1]
    Pca = Probability of correct assignment [0,1]
    Pw  = Probability of water (prior) [0,1]
    Ptf = Truth probability of pixel assignment 0.5 for non-informative
    ref_dem_std = 10
    method = composite (default), simple, water_fraction

    Reference: implements Eq.s (18), (26), and (30) in 
    "SWOT Hydrology Height and Area Uncertainty Estimation," 
    Brent Williams, 2018, JPL Memo
    
    TODO: add dark water uncertainty estimate
    TODO: add low coherence uncertainty estimate
    '''
    # get indicator functions
    Ide = np.zeros(np.shape(pixel_area))
    for water_edge_klass in water_edge_klasses:
        Ide[klass == water_edge_klass] = 1.0
    for land_edge_klass in land_edge_klasses:
        Ide[klass == land_edge_klass] = 1.0

    Pe = Ide  # use detected edge as probability of true edge pixels...

    I = np.zeros(np.shape(pixel_area))
    I[Ide > 0] = 1.0
    for interior_water_klass in interior_water_klasses:
        I[klass == interior_water_klass] = 1.0  # all pixels near water

    # get false and missed assignment rates from correct assignment rate 
    Pfa = 1 - Pca
    Pma = 1 - Pca

    # get the assignment rates, bias and variance 
    #Ptf = 0.5 # 0.5 equally likley to be in or out of the bin
    Pf = Pca * Ptf + Pfa *(1-Ptf)
    Bf = Pfa * (1-Ptf) - Pma * Pf
    Vf = Pfa*(1-Ptf) + Pma * Pf

    # handle pixel size uncertainty
    #ref_dem_std = 10 # m
    sigma_a = darea_dheight*ref_dem_std*pixel_area # 0.05* pixel_area

    # handle the sampling error
    sigma_s2 = 1.0/12.0

    # get detection rates
    if method == 'simple' or method == 'composite':
        Pcd = 1 - Pmd
        Pdw = Pcd*Pw + Pfd*(1-Pw)

        # also get
        Vdw = (Pfd*(1-Pw)+Pmd*Pw)
        V_dwf = (Pf*Vdw+Pw*Vf - 2*Pmd*Pw*Pfa*(1-Ptf))

        # get the aggregate sampling error
        # first term in Eq. 18  using Pe = Ie 
        var_samp_bar = simple(
            pixel_area[good]**2.0 * sigma_s2 * Pe[good] * Ptf, metric='sum')

        # handle the case where there are no edge pixels...
        num_pixels_edge = simple(Pe[good], metric='sum')
        if num_pixels_edge == 0:
            var_samp_bar = 0

        # the area uncertainty to be aggregated (2nd term in Eq. 12)
        var_pix_area_dw = sigma_a**2 * Pdw * Pf
        var_pix_area_dw_bar = simple(
            var_pix_area_dw[good]*I[good], metric='sum')

        # the detection and assignment rate uncertainty to be aggregated
        # 3rd term in Eq. 12
        var_area_dw = pixel_area**2 * V_dwf
        var_area_dw_bar = simple(var_area_dw[good] * I[good], metric='sum')

        # aggregate bias term (last term in Eq. 12)
        Bdw = Pfd * (1-Pw) - Pmd * Pw
        # use detected water prob for bias??
        Bdwf = Bdw * Ptf + Pw * Bf

        Bdwf_bar = simple(pixel_area[good] * Bdwf[good] * I[good], metric='sum')
        Bdwf2_bar = simple(
            pixel_area[good]**2 * Bdwf[good]**2 * I[good], metric='sum')
        B_term_dw = Bdwf_bar**2 - Bdwf2_bar #
        # sqrt of Eq. 12, std_dw = sigma_{A_f}
        std_dw = np.sqrt(
            var_samp_bar + var_pix_area_dw_bar + var_area_dw_bar + B_term_dw)

    if method == 'water_fraction' or method == 'composite':
        # use water_fraction (only for edge pixels if composite)
        # implements Eq. 26
        alpha2 = water_fraction_uncert
        alpha = water_fraction

        # get truth alpha and its var from noisy meausred alpha...
        sig_alpha2 = alpha2**2
        alpha_t = alpha

        Valpha=(sig_alpha2*Pf + alpha_t**2*Vf)
        var_area_alpha = pixel_area**2 * Valpha
        # 2nd term in Eq. 26
        var_area_alpha_bar = simple(var_area_alpha[good]*I[good], metric='sum')

        var_pix_area_alpha = sigma_a**2 * (sig_alpha2 + alpha**2) * Pf
        # 1st term in Eq. 26
        var_pix_area_alpha_bar = simple(
            var_pix_area_alpha[good]*I[good], metric='sum')

        B_tmp = pixel_area*alpha_t*Bf
        Balphaf_bar = simple(B_tmp[good] * I[good], metric='sum')
        Balphaf2_bar = simple(B_tmp[good]**2 * I[good], metric='sum')

        # 3rd term in Eq 26
        B_term_alpha = Balphaf_bar**2 - Balphaf2_bar
        # sqrt of Eq. 26, std_alpha = sigma_{A_f,alpha}
        std_alpha = np.sqrt(
            var_area_alpha_bar + var_pix_area_alpha_bar +
            B_term_alpha)#/abs(area_bar)

    if method == 'composite':
        # assume that Pde is constant over each feature
        N_edges = simple(Ide, metric='sum')
        if N_edges == 0:
            Pde = 0
        else:
            N_tot = simple(I, metric='sum')
            Pde = N_edges/N_tot
            if Pde>1:
                Pde =1
        Pde_x = np.zeros(np.shape(pixel_area))+Pde

        # fakely account for Pme < Pe by scaling by 10%
        var_samp_composite_bar = 0.1 * var_samp_bar
        var_area_composite_bar = (
            var_area_alpha_bar * Pde + (1-Pde) * var_area_dw_bar)
        var_pix_area_composite_bar = (
            var_pix_area_alpha_bar * Pde + (1-Pde) * var_pix_area_dw_bar)
        B_tmp = pixel_area*((1-Pde_x) * Bdwf + Pde_x* alpha_t*Bf)

        Bcomposite_bar = simple(B_tmp, metric='sum')
        Bcomposite2_bar = simple(B_tmp**2, metric='sum')
        B_term_composite = Bcomposite_bar**2 - Bcomposite2_bar
        # sqrt of Eq. 30, std_compsite = sigma_{A'_f}
        std_composite = np.sqrt(var_samp_composite_bar
                                + var_area_composite_bar
                                + var_pix_area_composite_bar
                                + B_term_composite)
        std_out = std_composite
    if method == 'simple':
        std_out = std_dw
    if method == 'water_fraction':
        std_out = std_alpha
    return std_out

def area_with_uncert(
        pixel_area, water_fraction, water_fraction_uncert, darea_dheight,
        klass, Pfd, Pmd, good, Pca=0.9, Pw=0.5, Ptf=0.5, ref_dem_std=10,
        interior_water_klasses=AGG_CLASSES['interior_water_klasses'],
        water_edge_klasses=AGG_CLASSES['water_edge_klasses'],
        land_edge_klasses=AGG_CLASSES['land_edge_klasses'],
        dark_water_klasses=AGG_CLASSES['dark_water_klasses'],
        method='composite'):

    area_agg, num_pixels = area_only(
        pixel_area, water_fraction, klass, good, method=method,
        interior_water_klasses=interior_water_klasses,
        water_edge_klasses=water_edge_klasses,
        land_edge_klasses=land_edge_klasses,
        dark_water_klasses=dark_water_klasses)

    area_unc = area_uncert(
        pixel_area, water_fraction, water_fraction_uncert, darea_dheight,
        klass, Pfd, Pmd, good, Pca=Pca, Pw=Pw, Ptf=Ptf, ref_dem_std=ref_dem_std,
        interior_water_klasses=interior_water_klasses,
        water_edge_klasses=water_edge_klasses,
        land_edge_klasses=land_edge_klasses, 
        dark_water_klasses=dark_water_klasses,
        method=method)

    # normalize to get area percent error
    area_pcnt_uncert = area_unc/abs(area_agg)*100.0
    return area_agg, area_unc, area_pcnt_uncert

def get_sensor_index(pixc):
    """ Return the sensor index for a pixel cloud from illumination time """
    f = interpolate.interp1d(pixc['tvp']['time'], range(len(pixc['tvp']['time'])))
    illumination_time = pixc['pixel_cloud']['illumination_time'].data[
        np.logical_not(pixc['pixel_cloud']['illumination_time'].mask)]
    sensor_index = (np.rint(f(illumination_time))).astype(int).flatten()
    return sensor_index

def flatten_interferogram(
        ifgram, plus_y_antenna_xyz, minus_y_antenna_xyz, target_xyz, tvp_index,
        wavelength):
    """ Return the flattened interferogram using provided geolocations """
    # Compute distance between target and sensor for each pixel
    dist_e = np.sqrt(
        (plus_y_antenna_xyz[0][tvp_index] - target_xyz[0])**2
        + (plus_y_antenna_xyz[1][tvp_index] - target_xyz[1])**2
        + (plus_y_antenna_xyz[2][tvp_index] - target_xyz[2])**2)

    dist_r = np.sqrt(
        (minus_y_antenna_xyz[0][tvp_index] - target_xyz[0])**2
        + (minus_y_antenna_xyz[1][tvp_index] - target_xyz[1])**2
        + (minus_y_antenna_xyz[2][tvp_index] - target_xyz[2])**2)

    # Compute the corresponding reference phase and flatten the interferogram
    phase_ref = -2*np.pi / wavelength*(dist_e - dist_r)
    interferogram_flatten  = ifgram*np.exp(-1.j*phase_ref)

    return interferogram_flatten

