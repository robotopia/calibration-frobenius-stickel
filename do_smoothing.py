import numpy as np
import matplotlib.pyplot as plt
import aocal
import aocal_plot
import argparse


class GeneralisedStickel:
    pass


def rotate_phases(ao1, theta_rad, chan=None):

    if np.isscalar(theta_rad) and chan is not None:
        ao1[:,:,chan,:] *= np.exp(1j*theta_rad)
        return ao1

    if np.isscalar(theta_rad) and chan is None:
        ao1 *= np.exp(1j*theta_rad)
        return ao1

    if not np.isscalar(theta_rad):
        # Attempt broadcasting
        ao1 *= np.exp(1j*theta_rad[np.newaxis, np.newaxis, :, np.newaxis])
        return ao1

def optimal_rotation(ao1, ao2, axis=None):
    mean = np.nanmean(ao1 * np.conj(ao2), axis=axis)
    return np.angle(mean)

def ao_min_diff(ao1, ao2, axis=None):
    '''
    Returns the minimum difference of two (sets of) calibration solutions,
    where "minimum" means "minimum up to rotation of the solutions by a unit
    length phasor":

        ao1 - e^(iθ) * ao2

    "axis" controls which dimensions are summed over when computing the
    optimum theta value.
    '''
    theta_min = optimal_rotation(ao1, ao2, axis=axis)
    diff = ao1 - rotate_phases(ao2, theta_min)

    return diff

def ao_min_diff2_freq(ao1, freqs=None):
    '''
    Returns the second order central finite difference over the frequency axis
    '''
    # If no freqs are given, construct one
    if freqs is None:
        freqs = np.arange(ao1.n_chan)

    # Pad the frequency axis by zeros on either side
    zero_slice = np.full((ao1.n_int, ao1.n_ant, 1, ao1.n_pol), 0.0 + 0.0j)
    padded_ao1 = np.concatenate((zero_slice, ao1, zero_slice), axis=2)

    # Do something similar with the frequencies
    df_start = freqs[1] - freqs[0]
    df_end = freqs[-1] - freqs[-2]
    padded_freqs = np.concatenate(([freqs[0] - df_start], freqs, [freqs[-1] + df_end]))

    # First order differences
    dfreq_12 = padded_freqs[np.newaxis, np.newaxis, 1:-1, np.newaxis] - padded_freqs[np.newaxis, np.newaxis, :-2, np.newaxis]  # f_n - f_{n-1}
    dfreq_23 = padded_freqs[np.newaxis, np.newaxis, 2:, np.newaxis] - padded_freqs[np.newaxis, np.newaxis, 1:-1, np.newaxis]  # f_{n+1} - f_n
    diff_12 = ao_min_diff(padded_ao1[:,:,1:-1,:], padded_ao1[:,:,:-2,:], axis=(0, 1, 3)) / dfreq_12
    diff_23 = ao_min_diff(padded_ao1[:,:,2:,:], padded_ao1[:,:,1:-1,:], axis=(0, 1, 3)) / dfreq_23

    # Second order difference
    dfreq_13 = padded_freqs[np.newaxis, np.newaxis, 2:, np.newaxis] - padded_freqs[np.newaxis, np.newaxis, :-2, np.newaxis]
    diff2_123 = (diff_23 - diff_12) / (0.5 * dfreq_13)

    return diff2_123


def test_optimal_rotation(ao):

    # Plot the original
    aocal_plot.plot(ao, plot_filename="orig", n_rows=6, ants_per_line=6)

    # Randomise the phases, and plot
    theta = np.random.random((ao.n_chan,)) * 2*np.pi
    theta[0] = 0.0
    rotate_phases(ao, theta)
    aocal_plot.plot(ao, plot_filename="rand", n_rows=6, ants_per_line=6)

    # Optimise across frequency
    # Remove nan frequencies
    delete_chans = np.isnan(np.angle(np.nanmean(ao, axis=(0,1,3))))
    chans = np.arange(ao.n_chan)

    ao = np.delete(ao, delete_chans, axis=2)
    chans = np.delete(chans, delete_chans)

    delta_theta_mins = np.array([0.0] + [-optimal_rotation(ao[:,:,i,:], ao[:,:,i-1,:]) for i in range(1, ao.shape[2])])
    theta_mins = np.cumsum(delta_theta_mins)

    # Rotate using the found optimal thetas
    ao = rotate_phases(ao, theta_mins)
    aocal_plot.plot(ao, plot_filename="test", n_rows=6, ants_per_line=6, chan_idxs=chans)


def objective(ao, ao_hat, lmbda, freqs):
    '''
    ao is the original data set
    ao_hat is the model
    lmbda is the regularisation parameter

    IMPORTANT: It is assumed that any channels which contain only nans have
    already been removed
    '''
    # Calculate the residuals
    R = ao_min_diff(ao, ao_hat, axis=(0, 1, 3)) # "axis" controls how theta_min is calculated

    # The "goodness-of-fit" cost function:
    C_fit = np.nanmean(np.abs(R)**2)  # Is actually (the square of) the Frobenius norm

    # The second derivative
    Diff2 = ao_min_diff2_freq(ao_hat, freqs=freqs)

    # The smoothness cost function
    C_smooth = np.nanmean(np.abs(Diff2)**2)

    # The final, regularised cost function
    return C_fit + lmbda*C_smooth


def test_objective(ao, lmbda):

    # Make a "model" which is just a slightly perturbed version of the original
    ao_hat = ao + 0.05*ao*np.exp(2j*np.pi*np.random.random(ao.shape))
    freqs = np.arange(ao.n_chan)

    return objective(ao, ao_hat, lmbda, freqs)


def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='Apply smoothing by regularisation to calibration solutions')

    parser.add_argument('solution_file', help='A calibration solution file in the "AOCal" format')
    parser.add_argument('--lmbda', type=float, default=1.0, help='The regularisation parameter. Small numbers make models that match the data, large numbers make smooth models. Finding the right balance is a case-by-case black art. [default = 1.0]')

    args = parser.parse_args()

    # Validate arguments
    assert args.lmbda > 0, "The 'lmbda' parameter must be > 0"

    # Open cal file and load the data
    ao = aocal.fromfile(args.solution_file)

    #test_optimal_rotation(ao)
    print(test_objective(ao, args.lmbda))


if __name__ == '__main__':
    main()
