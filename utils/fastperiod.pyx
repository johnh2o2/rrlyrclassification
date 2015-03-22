# FASTPERIOD - Cythonized versions of the following period finding algorithms
# from periodbase.py
# 1. dworetksy_period_find -> dworetsky
# 2. stetson_period_find -> stetson
# 3. lomb_scargle_periodogram -> lombscargle
# 4. spectral_window_function -> specwindow
#
# the average speedup is about 5-7x. not too bad.
#
# compile this using the following:
# $ python setup.py build_ext --inplace

import numpy as np
cimport numpy as np
from numpy import pi

def stetson_single_per(  np.ndarray[np.float64_t, ndim=1] t,
                         np.ndarray[np.float64_t, ndim=1] x, 
                         np.float64_t p, 
                         np.ndarray[np.float64_t, ndim=1] err):

    cdef:
        np.int64_t                          n_obs, j_range, j
        np.float64_t                        fold_t, strlen_numer, strlen_denom, strlen, weights
        np.ndarray[np.float64_t, ndim=1]    phase, phase_sorted, x_sorted, err_sorted
        np.ndarray[np.int64_t, ndim=1]      phase_sort_ind
    n_obs = len(x)
    j_range = n_obs-1

    fold_t = np.min(t) # fold at the first time element

    phase = (t - fold_t)/p - np.floor((t - fold_t)/p)
    phase_sort_ind = np.argsort(phase)

    phase_sorted = phase[phase_sort_ind]
    x_sorted = x[phase_sort_ind]
    
    
    err_sorted = err[phase_sort_ind]

    strlen_numer = 0.0
    strlen_denom = 0.0

    # build up the string                 
    for j in range(j_range):
        
        weights = ( ( (err_sorted[j])*(err_sorted[j]) + 
                      (err_sorted[j+1])*(err_sorted[j+1]) ) *
                    (phase_sorted[j+1] - phase_sorted[j] + 1.0/n_obs) )
        weights = 1.0/weights

        strlen_numer += weights*np.fabs(x_sorted[j] - x_sorted[j+1])
        strlen_denom += weights

    # add the last element
    weights = ( ( (err_sorted[-1])*(err_sorted[-1]) + 
                  (err_sorted[0])*(err_sorted[0]) ) *
                (phase_sorted[0] - phase_sorted[-1] + 1.0/n_obs) )
    weights = 1.0/weights

    strlen_numer += weights*np.fabs(x_sorted[-1] - x_sorted[0])
    strlen_denom += weights

    strlen = strlen_numer/strlen_denom

    return strlen

def dworetsky_single_per(np.ndarray[np.float64_t, ndim=1] t,
                         np.ndarray[np.float64_t, ndim=1] x, 
                         np.float64_t p):
    cdef:
        np.int64_t                          n_obs, j_range, j
        np.float64_t                        fold_t, strlen_numer, strlen_denom, strlen, weights
        np.ndarray[np.float64_t, ndim=1]    phase, mod_x, mod_x_sorted, phase_sorted
        np.ndarray[np.int64_t, ndim=1]      phase_sort_ind

    j_range = len(x)-1
    fold_t = np.min(t) # fold at the first time element

    mod_x = (x - np.min(x))/(2.0*(np.max(x) - np.min(x))) - 0.25

    phase = (t - fold_t)/p - np.floor((t - fold_t)/p)

    phase_sort_ind = np.argsort(phase)
    phase_sorted = phase[phase_sort_ind]
    mod_x_sorted = mod_x[phase_sort_ind]

    strlen = 0.0

    # now calculate the string length
    for j in range(j_range):
        strlen += np.sqrt( (mod_x_sorted[j+1] - mod_x_sorted[j]) * 
                           (mod_x_sorted[j+1] - mod_x_sorted[j]) +
                           (phase_sorted[j+1] - phase_sorted[j]) * 
                           (phase_sorted[j+1] - phase_sorted[j]))

    strlen += np.sqrt( (mod_x_sorted[0] - mod_x_sorted[-1]) * 
                       (mod_x_sorted[0] - mod_x_sorted[-1]) +
                       (phase_sorted[0] - phase_sorted[-1] + 1) *
                       (phase_sorted[0] - phase_sorted[-1] + 1))

    return strlen   

def dworetsky(np.ndarray[np.float64_t,ndim=1] time,
              np.ndarray[np.float64_t,ndim=1] mag,
              np.ndarray[np.float64_t,ndim=1] err,
              np.float64_t init_p,
              np.float64_t end_p,
              np.float64_t f_step):

    # declare ALL THE THINGS
    cdef np.ndarray[np.float64_t,ndim=1] mod_mag
    cdef np.ndarray[np.float64_t,ndim=1] phase
    cdef np.ndarray[long,ndim=1] phase_sort_ind
    cdef np.ndarray[np.float64_t,ndim=1] phase_sorted
    cdef np.ndarray[np.float64_t,ndim=1] mod_mag_sorted
    cdef np.ndarray[np.float64_t,ndim=1] out_periods
    cdef np.ndarray[np.float64_t,ndim=1] out_strlens

    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef long j_range

    cdef np.float64_t init_f
    cdef np.float64_t end_f
    cdef np.float64_t strlen
    cdef np.float64_t period
    cdef np.float64_t fold_time

    j_range = len(mag)-1
    fold_time = np.min(time) # fold at the first time element

    mod_mag = (mag - np.min(mag))/(2.0*(np.max(mag) - np.min(mag))) - 0.25
    
    init_f = 1.0/end_p
    end_f = 1.0/init_p

    n_freqs = np.ceil((end_f - init_f)/f_step)
    out_periods = np.zeros(n_freqs,dtype=np.float64)
    out_strlens = np.zeros(n_freqs,dtype=np.float64)

    for i in range(n_freqs):
        
        period = 1.0/init_f

        phase = (time - fold_time)/period - np.floor((time - fold_time)/period)
        
        phase_sort_ind = np.argsort(phase)
        phase_sorted = phase[phase_sort_ind]
        mod_mag_sorted = mod_mag[phase_sort_ind]

        strlen = 0.0
        
        # now calculate the string length
        for j in range(j_range):
            strlen += np.sqrt( (mod_mag_sorted[j+1] - mod_mag_sorted[j]) * 
                               (mod_mag_sorted[j+1] - mod_mag_sorted[j]) +
                               (phase_sorted[j+1] - phase_sorted[j]) * 
                               (phase_sorted[j+1] - phase_sorted[j]))

        strlen += np.sqrt( (mod_mag_sorted[0] - mod_mag_sorted[-1]) * 
                           (mod_mag_sorted[0] - mod_mag_sorted[-1]) +
                           (phase_sorted[0] - phase_sorted[-1] + 1) *
                           (phase_sorted[0] - phase_sorted[-1] + 1))

        out_periods[i] = period
        out_strlens[i] = strlen

            # print('j = %i, period = %f, strlen = %f' %
            #       (j,period,strlen))

        init_f += f_step

    return (out_periods,out_strlens)


def stetson(np.ndarray[np.float64_t,ndim=1] time,
            np.ndarray[np.float64_t,ndim=1] mag,
            np.ndarray[np.float64_t,ndim=1] err,
            np.float64_t init_p,
            np.float64_t end_p,
            np.float64_t f_step):

    cdef np.ndarray[np.float64_t,ndim=1] phase
    cdef np.ndarray[long,ndim=1] phase_sort_ind
    cdef np.ndarray[np.float64_t,ndim=1] phase_sorted
    cdef np.ndarray[np.float64_t,ndim=1] mag_sorted
    cdef np.ndarray[np.float64_t,ndim=1] out_periods
    cdef np.ndarray[np.float64_t,ndim=1] out_strlens

    cdef Py_ssize_t i
    cdef Py_ssize_t j
    cdef long n_obs
    cdef long j_range

    cdef np.float64_t init_f
    cdef np.float64_t end_f
    cdef np.float64_t strlen
    cdef np.float64_t period
    cdef np.float64_t strlen_numer
    cdef np.float64_t strlen_denom
    cdef np.float64_t weights
    cdef np.float64_t fold_time

    n_obs = len(mag)
    j_range = n_obs-1

    fold_time = np.min(time) # fold at the first time element
    
    init_f = 1.0/end_p
    end_f = 1.0/init_p

    n_freqs = np.ceil((end_f - init_f)/f_step)

    out_periods = np.zeros(n_freqs,dtype=np.float64)
    out_strlens = np.zeros(n_freqs,dtype=np.float64)

    for i in range(n_freqs):
        
        period = 1.0/init_f

        phase = (time - fold_time)/period - np.floor((time - fold_time)/period)
        phase_sort_ind = np.argsort(phase)

        phase_sorted = phase[phase_sort_ind]
        mag_sorted = mag[phase_sort_ind]
        err_sorted = err[phase_sort_ind]

        strlen_numer = 0.0
        strlen_denom = 0.0
        
        # build up the string                 
        for j in range(j_range):
            
            weights = ( ( (err_sorted[j])*(err_sorted[j]) + 
                          (err_sorted[j+1])*(err_sorted[j+1]) ) *
                        (phase_sorted[j+1] - phase_sorted[j] + 1.0/n_obs) )
            weights = 1.0/weights

            strlen_numer += weights*np.fabs(mag_sorted[j] - mag_sorted[j+1])
            strlen_denom += weights

        # add the last element
        weights = ( ( (err_sorted[-1])*(err_sorted[-1]) + 
                      (err_sorted[0])*(err_sorted[0]) ) *
                    (phase_sorted[0] - phase_sorted[-1] + 1.0/n_obs) )
        weights = 1.0/weights

        strlen_numer += weights*np.fabs(mag_sorted[-1] - mag_sorted[0])
        strlen_denom += weights

        strlen = strlen_numer/strlen_denom

        out_periods[i] = period
        out_strlens[i] = strlen
        
        init_f += f_step

    return (out_periods,out_strlens)


def lombscargle(np.ndarray[np.float64_t,ndim=1] time,
                np.ndarray[np.float64_t,ndim=1] mag,
                np.ndarray[np.float64_t,ndim=1] err,
                np.float64_t init_p,
                np.float64_t end_p,
                np.float64_t f_step):
    
    cdef np.ndarray[np.float64_t,ndim=1] out_periods
    cdef np.ndarray[np.float64_t,ndim=1] out_pgram
    cdef np.ndarray[np.float64_t,ndim=1] norm_time
    
    cdef Py_ssize_t i
    cdef long n_obs

    cdef np.float64_t init_f
    cdef np.float64_t end_f
    cdef np.float64_t period
    cdef np.float64_t stdev
    cdef np.float64_t power_top_cos
    cdef np.float64_t power_bot_cos
    cdef np.float64_t power_top_sin
    cdef np.float64_t power_bot_sin
    cdef np.float64_t power
    cdef np.float64_t normalized_periodogram
    cdef np.float64_t tau

    stdev = np.std(mag)
    n_obs = len(time)

    norm_time = time - np.min(time)
    
    init_f = 1.0/end_p
    end_f = 1.0/init_p

    n_freqs = np.ceil((end_f - init_f)/f_step)

    out_periods = np.zeros(n_freqs,dtype=np.float64)
    out_pgram = np.zeros(n_freqs,dtype=np.float64)

    for i in range(n_freqs):
    
        period = 1.0/init_f

        tau = ( (period/(4.0*pi)) *
                np.arctan( (np.sum( np.sin(4.0*pi/period) * norm_time)) /
                           (np.sum( np.cos(4.0*pi/period) * norm_time)) ) )

        power_top_cos = (np.sum(mag * np.cos((2.0*pi/period)*(norm_time-tau))) *
                         np.sum(mag * np.cos((2.0*pi/period)*(norm_time-tau))))
        power_bot_cos = np.sum( (np.cos((2.0*pi/period)*(norm_time-tau))) *
                                (np.cos((2.0*pi/period)*(norm_time-tau))) )

        power_top_sin = (np.sum(mag * np.sin((2.0*pi/period)*(norm_time-tau))) *
                         np.sum(mag * np.sin((2.0*pi/period)*(norm_time-tau))))
        power_bot_sin = np.sum( (np.sin((2.0*pi/period)*(norm_time-tau))) *
                                (np.sin((2.0*pi/period)*(norm_time-tau))) )

        power = 0.5 * ( (power_top_cos/power_bot_cos) + 
                        (power_top_sin/power_bot_sin) )

        normalized_periodogram = power/(stdev*stdev)

        out_pgram[i] = normalized_periodogram
        out_periods[i] = period

        init_f += f_step

    return (out_periods,out_pgram)


def specwindow(np.ndarray[np.float64_t,ndim=1] time,
               np.float64_t init_p,
               np.float64_t end_p,
               np.float64_t f_step):
    
    cdef np.ndarray[np.float64_t,ndim=1] out_periods
    cdef np.ndarray[np.float64_t,ndim=1] out_pgram
    cdef np.ndarray[np.float64_t,ndim=1] norm_time
    
    cdef Py_ssize_t i

    cdef np.float64_t init_f
    cdef np.float64_t end_f
    cdef np.float64_t period
    cdef np.float64_t power_top_cos
    cdef np.float64_t power_bot_cos
    cdef np.float64_t power_top_sin
    cdef np.float64_t power_bot_sin
    cdef np.float64_t power
    cdef np.float64_t tau

    norm_time = time - np.min(time)
    
    init_f = 1.0/end_p
    end_f = 1.0/init_p

    n_freqs = np.ceil((end_f - init_f)/f_step)

    out_periods = np.zeros(n_freqs,dtype=np.float64)
    out_pgram = np.zeros(n_freqs,dtype=np.float64)

    for i in range(n_freqs):
    
        period = 1.0/init_f

        tau = ( (period/(4.0*pi)) *
                np.arctan( (np.sum( np.sin(4.0*pi/period) * norm_time)) /
                           (np.sum( np.cos(4.0*pi/period) * norm_time)) ) )

        power_top_cos = (np.sum(1.0 * np.cos((2.0*pi/period)*(norm_time-tau))) *
                         np.sum(1.0 * np.cos((2.0*pi/period)*(norm_time-tau))))
        power_bot_cos = np.sum( (np.cos((2.0*pi/period)*(norm_time-tau))) *
                                (np.cos((2.0*pi/period)*(norm_time-tau))) )

        power_top_sin = (np.sum(1.0 * np.sin((2.0*pi/period)*(norm_time-tau))) *
                         np.sum(1.0 * np.sin((2.0*pi/period)*(norm_time-tau))))
        power_bot_sin = np.sum( (np.sin((2.0*pi/period)*(norm_time-tau))) *
                                (np.sin((2.0*pi/period)*(norm_time-tau))) )

        power = 0.5 * ( (power_top_cos/power_bot_cos) + 
                        (power_top_sin/power_bot_sin) )

        out_pgram[i] = power
        out_periods[i] = period
        init_f += f_step

    return (out_periods,out_pgram)
