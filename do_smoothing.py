import numpy as np
import matplotlib.pyplot as plt
import aocal
import aocal_plot
import argparse


class GeneralisedStickel:
    pass


def rotate_phases(ao, theta_rad, chan=None):

    if np.isscalar(theta_rad) and chan is not None:
        ao[:,:,chan,:] *= np.exp(1j*theta_rad)
        return

    if np.isscalar(theta_rad) and chan is None:
        ao *= np.exp(1j*theta_rad)
        return

    if not np.isscalar(theta_rad):
        # Attempt broadcasting
        ao *= np.exp(1j*theta_rad[np.newaxis, np.newaxis, :, np.newaxis])
        return

def optimal_rotation(ao1, ao2):
    return np.angle(np.nansum(ao1 * np.conj(ao2)))

def test_optimal_rotation(ao):

    # Plot the original
    aocal_plot.plot(ao, plot_filename="orig", n_rows=6, ants_per_line=6)

    # Randomise the phases, and plot
    theta = np.random.random((ao.n_chan,)) * 2*np.pi
    theta[0] = 0.0
    rotate_phases(ao, theta)
    aocal_plot.plot(ao, plot_filename="rand", n_rows=6, ants_per_line=6)

    # Optimise across frequency
    delta_theta_min = -np.array([0.0] + [optimal_rotation(ao[:,:,i,:], ao[:,:,i-1,:]) for i in range(1, ao.n_chan)])
    theta_min = np.nancumsum(delta_theta_min)
    #print(np.abs(theta))
    #print(np.abs(theta_min))

    # Rotate using the found optimal thetas
    rotate_phases(ao, theta_min)
    aocal_plot.plot(ao, plot_filename="test", n_rows=6, ants_per_line=6)

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
