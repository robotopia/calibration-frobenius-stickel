import numpy as np
import matplotlib.pyplot as plt
import aocal
import argparse


class GeneralisedStickel:
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
    sol = aocal.fromfile(args.solution_file)

    # The idea

if __name__ == '__main__':
    main()
