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

def ao_min_diff(ao1, ao2, idx=None):
    '''
    Returns the minimum difference of two (sets of) calibration solutions.
    "Minimum" up to rotation of the solutions by a unit length phasor:
    ao1 - e^(iθ) * ao2
    '''
    theta_min = optimal_rotation(ao1, ao2)
    diff = ao1 - rotate_phases(ao2, theta_rad)

    if idx is not None:
        diff /= np.diff(idx)

    return diff

def ao_min_diff2(ao1, ao2, idx=None):

    diff1 = ao_min_diff(ao1, ao2, idx)
    if idx is not None:
        diff1_idx += idx[:-1] + np.diff(idx)/2

    # UP TO HERE

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

def objective(ao, ao_hat):
    '''
    ao is the original data set, where 
    '''
    pass

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

    print(ao.shape)

    test_optimal_rotation(ao)


if __name__ == '__main__':
    main()
