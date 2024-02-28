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

    rotate_phases(ao, np.deg2rad(90))
    aocal_plot.plot(ao, plot_filename="test", n_rows=6, ants_per_line=6)

if __name__ == '__main__':
    main()
