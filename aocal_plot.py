#!/usr/bin/env python
import os, logging
from optparse import OptionParser #NB zeus does not have argparse!

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import pylab
from astropy.io import fits

import aocal

def get_tile_flavors(metafits):
    """
    return tile flavours ordered in AO order
    """
    hdus = fits.open(metafits)
    # This mirrors what cotter does (see metafitsfile.cpp MetaFitsFile::ReadTiles
    inputs = hdus[1].data
    tiles = inputs[inputs['pol'] == 'X']
    sort_ant = tiles['antenna'].argsort()
    return tiles['Flavors'][sort_ant]

def get_receiver_slot_order(metafits):
    """
    returns a dictionary of dictionaries which will give an AO ordinal index to each receiver and slot
    """
    # get tile metadata table in AO order (discard Y polarisations)
    hdus = fits.open(metafits)
    inputs = hdus[1].data
    tiles = inputs[inputs['pol'] == 'X']
    sort_ant = tiles['antenna'].argsort()
    tiles = tiles[sort_ant]

    receivers = set(tiles['Rx'])
    rec_slot_dict = {}
    for receiver in sorted(receivers):
        rec_slot_dict[receiver] = {}
        rec_slots = set(tiles[tiles['Rx'] == receiver]['Slot'])
        for slot in sorted(rec_slots):
            # tiles with this rec, slot
            ant = tiles[(tiles['Rx'] == receiver) & (tiles['Slot'] == slot)]['antenna']
            if not len(ant) == 1:
                raise RuntimeError("Rx %d Slot %d does not map to a single antenna, but to %s" % (receiver, slot, str(ant)))
            rec_slot_dict[receiver][slot] = ant[0]
    return rec_slot_dict

def iter_rec_slot(rec_slot_dict):
    """
    iterate over all receivers and slots in order
    """
    for receiver in sorted(rec_slot_dict.keys()):
        for slot in sorted(rec_slot_dict[receiver].keys()):
            yield rec_slot_dict[receiver][slot]

POLS = {0: "XX", 1: "XY", 2: "YX", 3: "YY"}
POL_COLOR = {"XX": "#0000FF", "XY": "#AAAAFF", "YX": "#FFAAAA", "YY": "#FF0000"}
POL_ZORDER = {"XX": 3, "XY": 1, "YX": 2, "YY": 4}

def nanaverage(a, axis=None, weights=None):
    """
    weighted average treating NaN as zero
    """
    if weights is None:
        weights = np.ones(a.shape)
    return np.nansum(a*weights, axis=axis)/np.nansum(weights, axis=axis)

def plot(ao, plot_filename=None, refant=None, n_rows=8, plot_title="", amp_max=None, format="png", outdir=None, ants_per_line=8, marker=',', markersize=2, verbose=0, metafits=None, chan_idxs=None):
    """
    plot aocal
    """

    if verbose == 1:
        logging.basicConfig(level=logging.INFO)
    elif verbose > 1:
        logging.basicConfig(level=logging.DEBUG)

    if chan_idxs is None:
        chan_idxs = np.arange(ao.n_chan)

    n_cols = ao.n_ant//n_rows
    gs = gridspec.GridSpec(n_rows, n_cols)
    gs.update(hspace=0.0, wspace=0.0)
    ao_amp = np.abs(ao)
    if refant is not None:
        if isinstance(refant, int):
            if refant < 0:
                logging.info("using average as reference antenna")
                ant_avg = nanaverage(ao, axis=1, weights=ao_amp**-2)#correct phase
                ant_avg /= np.abs(ant_avg) # normalised
                print(nanaverage(ao_amp, axis=1, weights=ao_amp**-2))
                ant_avg *= nanaverage(ao_amp, axis=1, weights=ao_amp**-2) #standard scalar avg for amp
                ao = ao / ant_avg[:, np.newaxis, :, :]
                plot_title += " refant=average"
            else:
                logging.info("using antenna %d as reference antenna", refant)
                ao = ao / ao[:, refant, :, :][:, np.newaxis, :, :]
                plot_title += " refant=%d" % refant
        else:
                logging.info("using %s average as reference antenna", refant)
                flavors = get_tile_flavors(metafits).reshape(1, -1)
                if not refant in flavors:
                    raise RuntimeError("refant flavor not found")
                # ao[flavors == refant] has shape (n_refant, n_freq, n_pol)
                ant_avg = nanaverage(ao[flavors == refant], axis=0, weights=ao_amp[flavors == refant]**-2) # correct phase
                print(ant_avg.shape)
                ao = ao / ant_avg[np.newaxis, np.newaxis, :, :]
                plot_title += " refant=%s" % refant
    else:
        logging.info("no reference antenna")

    # Scale amplitude plots to the same, maximum gain 
    if amp_max is None:
        amp_max = 2*np.median(np.nan_to_num(np.abs(ao[..., [0,-1]])))
        logging.info("amp_max=%.1f" % amp_max)

    # Order by receiver/slot if metafits is requested
    if metafits is not None:
        rec_slot_dict = get_receiver_slot_order(metafits)
        ant_iter = iter_rec_slot(rec_slot_dict)
    else:
        ant_iter = range(ao.n_ant)

    for timestep in range(ao.n_int):

        phsfig = plt.figure(figsize=(24.0, 13.5))
        ampfig = plt.figure(figsize=(24.0, 13.5))
        logging.debug("vertical_index, horizontal_index")
        for a, antenna in enumerate(ant_iter):
            vertical_index = a // n_cols
            horizontal_index = a % n_cols
            logging.debug("%02d,%02d", vertical_index, horizontal_index)

            # Phase plot
            ax = phsfig.add_subplot(gs[vertical_index, horizontal_index])

            # Amplitude plot
            ax1 = ampfig.add_subplot(gs[vertical_index, horizontal_index])

            ax.text(0.05, 0.05, f'Ant #{antenna}',
                    horizontalalignment='left',
                    verticalalignment='bottom',
                    transform=ax.transAxes)
            ax1.text(0.05, 0.05, f'Ant #{antenna}',
                     horizontalalignment='left',
                     verticalalignment='bottom',
                     transform=ax1.transAxes)

            for pol in range(ao.n_pol):
                polstr = POLS[pol]
                amps = np.abs(ao[timestep, antenna, :, pol])
                angles= np.angle(ao[timestep, antenna, :, pol], deg=True)

                # Phase plot
                ax.plot(chan_idxs, angles, color=POL_COLOR[polstr], zorder=POL_ZORDER[polstr], linestyle='None', marker=marker, markersize=markersize, label=polstr)
                # Amplitude plot
                ax1.plot(chan_idxs, amps, color=POL_COLOR[polstr], zorder=POL_ZORDER[polstr], linestyle='None', marker=marker, markersize=markersize, label=polstr)

                #if np.all(np.isnan(amps)):
                    # all flagged
                    #rect = ax.patch  # a Rectangle instance
                    #rect.set_facecolor('black')
                    #rect.set_alpha('0.4')
                    #rect = ax1.patch  # a Rectangle instance
                    #rect.set_facecolor('black')
                    #rect.set_alpha('0.4')
            ax.set_autoscale_on(False)
            ax.set_xlim([-1, chan_idxs[-1] + 1])
            ax.set_ylim([-180,180])
            if horizontal_index == 0:
                ax.set_yticks([-90,0,90])
                ax.set_ylabel("Phase (deg)")
                ax1.set_ylabel("Amplitude")
            else:
                ax.set_yticks([])
                ax1.set_yticks([])
            if vertical_index == n_rows - 1:
                ax.set_xlabel("Frequency")
                ax1.set_xlabel("Frequency")
            else:
                ax.set_xticks([])
                ax1.set_xticks([])

            ax1.set_xlim([-1, chan_idxs[-1] + 1])
            ax1.set_ylim([0, amp_max])
            #ax1.set_autoscale_on(False)

        phsfig.suptitle(plot_title,fontsize=16)
        ampfig.suptitle(plot_title, fontsize=16)

        if plot_filename is None:
            ampfig.show()
            phsfig.show()
        else:
            if outdir is not None:
                plot_filename = os.path.join(outdir, os.path.basename(plot_filename))
            if ao.n_int > 1:
                int_str = "_t%04d" % timestep
            else:
                int_str = ""
            ampfig.tight_layout()
            phsfig.tight_layout()
            ampfig.savefig("%s%s_amp.%s" % (plot_filename, int_str, format))
            phsfig.savefig("%s%s_phase.%s" % (plot_filename, int_str, format))

if __name__ == '__main__':
    parser = OptionParser(usage = "usage: %prog binfile" +
    """
    Plot calibration solutions
    """)
    parser.add_option("--refant", default=None, dest="refant", type="int", help="divide solutions through by reference antenna. Negative means divide through by mean antenna")
    parser.add_option("--ref_flavor", default=None, dest="ref_flavor", type="str", help="divide solutions through by all antennas of given flavour. Negative means divide through by mean antenna")
    parser.add_option("-m", "--metafits", default=None, dest="metafits", help="metafits file (for ordering by receiver or ref_flavour)")
    parser.add_option("-v", "--verbose", action="count", default=0, dest="verbose", help="-v info, -vv debug")
    parser.add_option("--outdir", default=None, dest="outdir", help="output directory [default: same as binfile]")
    parser.add_option("--title", default="", dest="plot_title", help="plot title")
    parser.add_option("--format", default="png", dest="format", help="plot format [default: %default]")
    parser.add_option("--suffix", default="", dest="suffix", help="suffix to add to plot names")
    parser.add_option("--amp_max", default=None, dest="amp_max", type="float", help="Maximum of y axis of amplitude plots")
    parser.add_option("--marker", default=',', dest="marker", type="string", help="matplotlib marker [default: %default]")
    parser.add_option("--markersize", default=2, dest="markersize", type="int", help="matplotlib markersize [default: %default]")
    opts, args = parser.parse_args()

    if len(args) != 1:
        parser.error("incorrect number of arguments")
    if opts.ref_flavor is not None and opts.metafits is None:
        parser.error("ref_flavor requires metafits")
    if opts.refant is not None and opts.ref_flavor is not None:
        logging.warn("refant and ref_flavor given, using refant")
    ref = opts.refant if opts.refant is not None else opts.ref_flavor


    ao = aocal.fromfile(args[0])
    plot(ao, os.path.splitext(args[0])[0]+opts.suffix, ref, plot_title = opts.plot_title, outdir=opts.outdir, format=opts.format, amp_max=opts.amp_max, marker=opts.marker, markersize=opts.markersize, verbose=opts.verbose, metafits=opts.metafits)
