#!/usr/bin/env python3

"""
Display NPZ data from stationMaster in an interactive GUI sort of way.
"""

import os
import sys
import numpy as np
import argparse
import tempfile

from lsl.common import stations
from lsl.misc.mathutils import to_dB

import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import matplotlib
matplotlib.use('TkAgg')
matplotlib.interactive(True)

from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk, FigureCanvasTkAgg
from matplotlib.figure import Figure


class TBW_GUI(object):
    """
    Object responsible for drawing and interacting with the two matplotlib
    canvas in the main GUI.

    Arguments:
    * frame - parent window used for displaying

    Keyword:
    * antennas - list of Antenna objects that compose the stations
    * freq - array of frequencies (in Hz) for the spectral data
    * spec - 2-D (antennas by channels) array of spectra
    * specTemplate - master spectrum for comparisons
    """

    def __init__(self, frame, antennas=None, freq=None, spec=None, specTemplate=None, resFreq=None):
        self.frame = frame
        self.press = None
        self.bestX = -1
        self.bestY = 0

        self.filename = ''
        self.color = 0
        self.antennas = antennas
        self.freq = freq
        self.spec = spec
        self.specTemplate = specTemplate
        self.resFreq = resFreq
        self.avgPower = None
        self.dataRange = None

        self.ax1 = None
        self.ax2 = None
        self.oldMark = None

    def loadData(self, filename):
        """
        Given the filename of an NPZ file created by stationMaster.py, load
        in the NPZ file and set all of the various data attributes used for
        plotting.
        """

        dataDict = np.load(filename)
        self.freq = dataDict['freq']
        masterSpectra = dataDict['masterSpectra']

        self.spec = masterSpectra.mean(axis=0)
        self.specTemplate = np.median(self.spec, axis=0)
        self.resFreq = dataDict['resFreq']
        try:
            self.avgPower = dataDict['avgPower']
        except KeyError:
            self.avgPower = None
        try:
            self.dataRange = dataDict['dataRange']
        except KeyError:
            self.dataRange = None

        try:
            self.adcHistogram = dataDict['adcHistogram']
        except KeyError:
            self.adcHistogram = None
        # Set the station
        try:
            ssmifContents = dataDict['ssmifContents']
            if ssmifContents.shape == ():
                station = stations.lwa1
                self.station = station.name
                self.antennas = station.antennas
            else:
                fh, tempSSMIF = tempfile.mkstemp(suffix='.txt', prefix='ssmif-')
                fh = open(tempSSMIF, 'w')
                for line in ssmifContents:
                    fh.write('%s\n' % line)
                fh.close()

                station = stations.parse_ssmif(tempSSMIF)
                self.station = station.name
                self.antennas = station.antennas
                os.unlink(tempSSMIF)

        except KeyError:
            station = stations.lwa1
            self.station = station.name
            self.antennas = station.antennas

        # Set default colorbars
        self.limits = []
        self.limits.append([0, 2])      # 0 - Median Comparison
        self.limits.append([1, 2])      # 1 - RFI-46 Index
        self.limits.append([1, 2])      # 2 - RFI-64 Index
        self.limits.append([1, 160])    # 3 - RFI-76 Index
        self.limits.append([0, 3])      # 4 - Wiggle Index
        self.limits.append([1, 3])      # 5 - Antenna Status
        self.limits.append([31, 50])    # 6 - Resonance Point (MHz)

        # Save the filename and data
        path, basename = os.path.split(filename)
        self.filename = basename
        self.date = dataDict['date']
        self.date = self.date.tostring().decode()

        try:
            self.disconnect()
        except:
            pass

        # Clear the mark and plots
        self.oldMark = None
        self.frame.figure1.clf()
        self.frame.figure2.clf()

        self.connect()

    def draw(self):
        """
        Draw the station stand field using the selected colorization scheme.
        The default scheme is based on the mean ratio between the spectrum
        and the master template (self.specTemplate) between 32 and 50 MHz.
        """

        self.frame.config(cursor='watch')
        self.frame.update()

        ## Get Stand positions from the Antenna objects.  Select only one
        ## polarization since that is all we need
        standPos = np.array([[ant.stand.x, ant.stand.y, ant.stand.z] for ant in self.antennas if ant.pol == 0])

        ## Get the stand quality figure-of-merit
        compLow = 32e6
        compHigh = 50e6

        if self.color == 0:
            # Color by mean ratio between the spectrum and the master template
            # (self.specTemplate) between 32 and 50 MHz (default)
            specDiff = np.zeros(self.spec.shape[0])
            toCompare = np.where( (self.freq>compLow) & (self.freq<compHigh) )[0]
            for i in range(self.spec.shape[0]):
                specDiff[i] = (self.spec[i,toCompare] / self.specTemplate[toCompare]).mean()

            cbTitle = '%.0f to %.0f MHz Mean Deviation' % (compLow/1e6, compHigh/1e6)
        elif self.color == 1:
            # Color by the value of the RFI-46 index.  This index is the maximum
            # ratio of the spectrum and the master template between 45 and 47 MHz.
            # The value of RFI-46 is also corrected for any systematic offset
            # between the spectrum and the template by looking at the 75 to 77 MHz
            # region.
            specDiff = np.zeros(self.spec.shape[0])
            rfi1 = np.where( (self.freq>45e6) & (self.freq<47e6) )[0]
            corr = np.where( (self.freq>75e6) & (self.freq<77e6) )[0]

            for i in range(self.spec.shape[0]):
                specDiff[i] = (self.spec[i,rfi1] / self.specTemplate[rfi1]).max()
                specDiff[i] /= (self.spec[i,corr] / self.specTemplate[corr]).mean()

            cbTitle = 'RFI-46 Index'
        elif self.color == 2:
            # Color by the value of the RFI-64 index.  This index is the maximum
            # ratio of the spectrum and the master template between 63 and 65 MHz.
            # The value of RFI-64 is also corrected for any systematic offset
            # between the spectrum and the template by looking at the 75 to 77 MHz
            # region.
            specDiff = np.zeros(self.spec.shape[0])
            rfi2 = np.where( (self.freq>63e6) & (self.freq<65e6) )[0]
            corr = np.where( (self.freq>75e6) & (self.freq<77e6) )[0]

            for i in range(self.spec.shape[0]):
                specDiff[i] = (self.spec[i,rfi2] / self.specTemplate[rfi2]).max()
                specDiff[i] /= (self.spec[i,corr] / self.specTemplate[corr]).mean()

            cbTitle = 'RFI-64 Index'
        elif self.color == 3:
            # Color by the value of the RFI-76 index.  This index is the maximum
            # ratio of the spectrum and the master template between 79 and 81 MHz.
            # The value of RFI-64 is also corrected for any systematic offset
            # between the spectrum and the template by looking at the 63 to 65 MHz
            # region.
            specDiff = np.zeros(self.spec.shape[0])
            rfi2 = np.where( (self.freq>75e6) & (self.freq<77e6) )[0]
            corr = np.where( (self.freq>63e6) & (self.freq<65e6) )[0]

            for i in range(self.spec.shape[0]):
                specDiff[i] = (self.spec[i,rfi2] / self.specTemplate[rfi2]).max()
                specDiff[i] /= np.median(self.spec[i,rfi2] / self.specTemplate[rfi2]).max()

            cbTitle = 'RFI-76 Index'
        elif self.color == 4:
            # Color by the wiggle index.  This is determined by fitting a line to
            # the log(spectra) between 65 and 70 MHz and looking at the RMS.
            specDiff = np.zeros(self.spec.shape[0])
            wgl = np.where( (self.freq>65e6) & (self.freq<70e6) )[0]
            for i in range(self.spec.shape[0]):
                junk = np.log10( self.spec[i,wgl] / self.specTemplate[wgl] )
                specDiff[i] = junk.std()
                specDiff[i] = 17 - int(self.antennas[i].arx.asp_channel % 16) + 1
                specDiff[i] *= self.antennas[i].cable.length / 10.0
                print(i, specDiff[i])

            cbTitle = 'Wiggle Index'
        elif self.color == 5:
            # Color by antenna status code.
            specDiff = np.zeros(self.spec.shape[0])
            for i in range(self.spec.shape[0]):
                specDiff[i] = self.antennas[i].status

            cbTitle = 'Antenna Status'
        else:
            # Color by the estimated resonance point frequency.  This is done
            # by finding the best-fit polynomial in between orders 3 and 12
            # for the 31 to 70 MHz spectral region.  The best-fit polynomial is
            # then evaluated to find its maximum value and that is used as the
            # resonance point.
            if self.resFreq is None:
                specDiff = np.zeros(self.spec.shape[0])
                toCompare = np.where( (self.freq>31e6) & (self.freq<70e6) )[0]
                for i in range(self.spec.shape[0]):
                    bestOrder = 0
                    bestRMS = 1e34
                    for j in range(3, 12):
                        coeff = np.polyfit(self.freq[toCompare]/1e6, to_dB(self.spec[i,toCompare]), j)
                        fit = np.polyval(coeff, self.freq[toCompare]/1e6)
                        rms = ((fit - to_dB(self.spec[i,toCompare]))**2).sum()
                        if rms < bestRMS:
                            bestOrder = j
                            bestRMS = rms

                    coeff = np.polyfit(self.freq[toCompare]/1e6, to_dB(self.spec[i,toCompare]), bestOrder)
                    fit = np.polyval(coeff, self.freq[toCompare]/1e6)
                    specDiff[i] = self.freq[toCompare[np.where( fit == fit.max() )[0]]] / 1e6
                self.resFreq = specDiff
            else:
                specDiff = self.resFreq

            cbTitle = 'Est. Resonance Point (MHz)'

        # Clip range
        specDiff = np.where( specDiff < self.limits[self.color][1], specDiff, self.limits[self.color][1])
        specDiff = np.where( specDiff > self.limits[self.color][0], specDiff, self.limits[self.color][0])

        self.frame.figure1.clf()
        self.ax1 = self.frame.figure1.gca()
        self.ax1.set_aspect('equal')
        # Stands
        m = self.ax1.scatter(standPos[:,0], standPos[:,1]+0.8, c=specDiff[0::2], s=45.0, alpha=0.80, marker='^')
        self.ax1.scatter(standPos[:,0], standPos[:,1]-0.8, c=specDiff[1::2], s=45.0, alpha=0.80, marker='v')

        if self.station == 'LWA1':
            ## Add the fence as a dashed line
            self.ax1.plot([-59.827, 59.771, 60.148, -59.700, -59.827],
                          [59.752, 59.864, -59.618, -59.948, 59.752], linestyle='--', color='k')

            ## Add the shelter
            self.ax1.plot([55.863, 58.144, 58.062, 55.791, 55.863],
                          [45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')

            ## Set the limits to just zoom in on the main station and the plot title
            self.ax1.set_xlim([-65, 65])
            self.ax1.set_ylim([-65, 65])

        elif self.station == 'LWASV':
            ## Add the awning
            self.ax1.plot([-30.625, -20.911, -13.498, -23.212, -30.625],
                          [107.211, 78.239, 81.013, 109.368, 107.211], linestyle='-', color='k')

            ## Add the shelter
            self.ax1.plot([-29.347, -28.836, -22.956, -23.467, -29.347],
                          [86.869, 84.712, 86.253, 88.41, 86.869], linestyle='-', color='k')

            ## Set the limits to just zoom in on the main station and the plot title
            self.ax1.set_xlim([-75, 75])
            self.ax1.set_ylim([-50,100])

        if self.date is None:
            self.ax1.set_title("Filename: '%s'" % self.filename)
        else:
            self.ax1.set_title('Date: UT %s' % self.date)

        ## Set the color bar, its title, and the axis labels
        cm = self.frame.figure1.colorbar(m, ax=self.ax1)
        cm.ax.set_ylabel(cbTitle)
        if cbTitle == 'Antenna Status':
            cm.set_ticks([1, 2, 3])
            cm.set_ticklabels(['Bad', 'Suspect', 'Good'])
        self.ax1.set_xlabel(r'$\Delta$ X [m]')
        self.ax1.set_ylabel(r'$\Delta$ Y [m]')

        if self.oldMark is not None and self.bestX != -1:
            # Redraw the mark at the previously selected stand position
            xy = [self.antennas[self.bestX-1].stand.x, self.antennas[self.bestX-1].stand.y]
            self.oldMark = self.ax1.plot([xy[0], xy[0]], [xy[1], xy[1]], linestyle=' ', marker='o', ms=15.0, mfc='None', color='k')

        ## Draw it
        self.frame.canvas1.draw()
        self.frame.canvas1.flush_events()

        self.frame.config(cursor='')

    def drawSpectrum(self, clickX, clickY, preferStand=None):
        """
        Get the spectra (both polarizations) for the antennas connected to
        the selected stand.
        """

        if preferStand is None:
            ## Figure out who is who and which antennas are closest to the
            ## clicked point.  This can be a little slow so the results are
            ## saved to the bestX and bestX attributes (depending on pol.)
            dist = 1e9
            for ant in self.antennas:
                cDist = np.sqrt( (ant.stand.x - clickX)**2 + (ant.stand.y - clickY)**2 )
                if cDist <= dist:
                    dist = cDist
                    if ant.pol == 0:
                        self.bestX = ant.digitizer
                    else:
                        self.bestY = ant.digitizer

        else:
            ## Right now 259 and 260 are at 0,0,0 and sit on top of each other.  Using
            ## the preferStand keyword, we can break this at least for searches
            for ant in self.antennas:
                if preferStand == ant.stand.id:
                    if ant.pol == 0:
                        self.bestX = ant.digitizer
                    else:
                        self.bestY = ant.digitizer

        ## Plot the spectra.  This plot includes the median composite
        ## (self.specTemplate) in green, the X polarization in blue, and
        ## the Y polarization in red.
        self.frame.figure2.clf()
        self.ax2 = self.frame.figure2.gca()
        self.ax2.plot(self.freq/1e6, to_dB(self.specTemplate), alpha=0.6, color='g', label='Composite')
        self.ax2.plot(self.freq/1e6, to_dB(self.spec[self.bestX-1,:]), color='b', label='X')
        self.ax2.plot(self.freq/1e6, to_dB(self.spec[self.bestY-1,:]), color='r', label='Y')

        ## Set the title, axis labels and add a legend
        self.ax2.set_title('Stand #%i' % self.antennas[self.bestX-1].stand.id)
        self.ax2.set_xlabel('Frequency [MHz]')
        self.ax2.set_ylabel('PSD [dB/RBW]')
        self.ax2.legend(loc=0)

        ## Draw and save the click (Why?)
        self.frame.canvas2.draw()
        self.xClick = clickX
        self.yClick = clickY

    def makeMark(self, clickX, clickY):
        """
        Mark the closest stand to the clicked point.  This needs to be called
        after drawSpectrum() since that function figures out which stand is
        closest.
        """

        if self.oldMark is not None:
            try:
                for line in self.oldMark:
                    line.remove()
            except:
                pass

        ## Figure out who is who
        xy = [self.antennas[self.bestX-1].stand.x, self.antennas[self.bestX-1].stand.y]

        self.oldMark = self.ax1.plot([xy[0], xy[0]], [xy[1], xy[1]], linestyle=' ', marker='o', ms=15.0, mfc='None', color='k')

        ## Set the limits to just zoom in on the main stations
        if self.station == 'LWA1':
            self.ax1.set_xlim([-65, 65])
            self.ax1.set_ylim([-65, 65])
        elif self.station == 'LWASV':
            self.ax1.set_xlim([-75, 75])
            self.ax1.set_ylim([-50,100])

        self.frame.canvas1.draw()

    def connect(self):
        """
        Connect to all the events we need to interact with the plots.
        """

        self.cidpress   = self.frame.figure1.canvas.mpl_connect('button_press_event',  self.on_press)
        self.cidmotion  = self.frame.figure1.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cidmotion2 = self.frame.figure2.canvas.mpl_connect('motion_notify_event', self.on_motion2)

    def on_press(self, event):
        """
        On a button press we will update the spectrum and mark the closest
        stand.
        """

        if event.inaxes and self.frame.toolbar.mode == '':
            clickX = event.xdata
            clickY = event.ydata

            self.drawSpectrum(clickX, clickY)
            self.makeMark(clickX, clickY)

    def on_motion(self, event):
        """
        Deal with motion events in the stand field window.  This involves
        setting the status bar with the current x and y coordinates as well
        as the stand number of the selected stand (if any).
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            self.frame.statusbar.config(text="x=%.3f m, y=%.3f m" % (clickX, clickY))
        else:
            self.frame.statusbar.config(text="")

    def on_motion2(self, event):
        """
        Deal with motion events in the spectrum window.  This involves
        setting the status bar with the current frequency and power values
        as well as the stand number of the selected.  If not stand has been
        selected, nothing is shown.
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            if self.bestX == -1:
                self.frame.statusbar.config(text="")
            else:
                self.frame.statusbar.config(text="freq=%.3f MHz, PSD=%.3f dB/RBW" % (clickX, clickY))
        else:
            self.frame.statusbar.config(text="")


    def disconnect(self):
        """
        Disconnect all the stored connection ids.
        """

        self.frame.figure1.canvas.mpl_disconnect(self.cidpress)
        self.frame.figure1.canvas.mpl_disconnect(self.cidmotion)
        self.frame.figure2.canvas.mpl_disconnect(self.cidmotion2)


class MainWindow(tk.Tk):
    def __init__(self):
        super().__init__()

        self.dirname = ''
        self.data = None
        self.cAdjust = None

        self.title("Station Master GUI")
        self.geometry("1200x600")

        self.initUI()
        self.initEvents()

    def initUI(self):
        """
        Start the user interface.
        """

        # Status bar at bottom
        self.statusbar = ttk.Label(self, text="", relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(fill=tk.X, side=tk.BOTTOM)

        # Menu bar
        menubar = tk.Menu(self)

        fileMenu = tk.Menu(menubar, tearoff=0)
        colorMenu = tk.Menu(menubar, tearoff=0)
        detailMenu = tk.Menu(menubar, tearoff=0)
        powerMenu = tk.Menu(menubar, tearoff=0)
        selectMenu = tk.Menu(menubar, tearoff=0)

        ## File menu
        fileMenu.add_command(label='Open', command=self.onOpen, accelerator='Ctrl+O')
        fileMenu.add_command(label='Show SSMIF Status', command=self.onSSMIF)
        fileMenu.add_separator()
        fileMenu.add_command(label='Quit', command=self.onQuit, accelerator='Ctrl+Q')

        ## Color menu with radio buttons
        self.color_var = tk.IntVar(value=0)
        colorMenu.add_radiobutton(label='Median Comparison', variable=self.color_var,
                                   value=0, command=self.onColorChange)
        colorMenu.add_radiobutton(label='Resonance Point', variable=self.color_var,
                                   value=6, command=self.onColorChange)
        colorMenu.add_radiobutton(label='RFI-46 Index', variable=self.color_var,
                                   value=1, command=self.onColorChange)
        colorMenu.add_radiobutton(label='RFI-64 Index', variable=self.color_var,
                                   value=2, command=self.onColorChange)
        colorMenu.add_radiobutton(label='RFI-76 Index', variable=self.color_var,
                                   value=3, command=self.onColorChange)
        colorMenu.add_radiobutton(label='Wiggle Index', variable=self.color_var,
                                   value=4, command=self.onColorChange)
        colorMenu.add_radiobutton(label='Antenna Status', variable=self.color_var,
                                   value=5, command=self.onColorChange)
        colorMenu.add_separator()
        colorMenu.add_command(label='Adjust Contrast', command=self.onAdjust)

        ## Detail menu
        detailMenu.add_command(label='Antenna', command=self.onAntenna)
        detailMenu.add_command(label='Stand', command=self.onStand)
        detailMenu.add_command(label='FEE', command=self.onFEE)
        detailMenu.add_command(label='Cable', command=self.onCable)
        detailMenu.add_separator()
        detailMenu.add_command(label='Shelter RFI Index', command=self.onRFI)
        detailMenu.add_separator()
        detailMenu.add_command(label='Change Antenna/FEE Status', command=self.onStatus)

        ## Power menu
        powerMenu.add_command(label='Plot ADC Histogram', command=self.onHistogram)
        powerMenu.add_command(label='Plot Power', command=self.onavgPower)
        powerMenu.add_command(label='Plot Data Range', command=self.onDataRange)
        powerMenu.add_separator()
        powerMenu.add_command(label='Summary', command=self.onavgPowerSummary)

        ## Select menu
        selectMenu.add_command(label='Antenna ID', command=self.onSelectAntenna)
        selectMenu.add_command(label='Stand ID', command=self.onSelectStand)
        selectMenu.add_separator()
        selectMenu.add_command(label='Digitizer Number', command=self.onSelectDigitizer)

        # Add menus to menubar
        menubar.add_cascade(label='File', menu=fileMenu)
        menubar.add_cascade(label='Color Coding', menu=colorMenu)
        menubar.add_cascade(label='Details', menu=detailMenu)
        menubar.add_cascade(label='Average Power', menu=powerMenu)
        menubar.add_cascade(label='Find', menu=selectMenu)
        self.config(menu=menubar)

        # Main content area with two plot panels
        main_frame = ttk.Frame(self)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Left panel - station field plot
        panel1 = ttk.Frame(main_frame)
        panel1.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.figure1 = Figure(figsize=(6,6))
        self.canvas1 = FigureCanvasTkAgg(self.figure1, master=panel1)
        self.canvas1.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Right panel - spectrum plot with toolbar
        panel2 = ttk.Frame(main_frame)
        panel2.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.figure2 = Figure(figsize=(6,4), tight_layout=True)
        self.canvas2 = FigureCanvasTkAgg(self.figure2, master=panel2)
        self.canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Toolbar for the spectrum plot
        toolbar_frame = ttk.Frame(panel2)
        toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas2, toolbar_frame)
        self.toolbar.update()

    def initEvents(self):
        """
        Set all of the various events in the main window.
        """

        # Keyboard shortcuts
        self.bind('<Control-o>', lambda e: self.onOpen())
        self.bind('<Control-q>', lambda e: self.onQuit())

        # Keyboard events for navigation
        self.bind('<KeyRelease>', self.onKeyPress)

        # Window resize
        self._resize_job = None
        self.bind('<Configure>', self._on_configure)

        # Window close
        self.protocol("WM_DELETE_WINDOW", self.onQuit)

    def _on_configure(self, event):
        """Handle window resize with debouncing."""
        # Only handle resize events for the main window, not child widgets
        if event.widget == self and self.data is not None:
            if self._resize_job:
                self.after_cancel(self._resize_job)
            self._resize_job = self.after(150, self._do_resize)

    def _do_resize(self):
        """Perform the actual resize - just redraw canvases."""
        try:
            self.canvas1.draw_idle()
            self.canvas2.draw_idle()
        except:
            pass

    def onOpen(self):
        """
        Open and load in a new NPZ file created by stationMaster.py.
        """

        filepath = filedialog.askopenfilename(
            initialdir=self.dirname,
            filetypes=[('NPZ Files', '*.npz'), ('All Files', '*.*')]
        )
        if filepath:
            self.dirname = os.path.dirname(filepath)
            self.filename = os.path.basename(filepath)
            self.data = TBW_GUI(self)
            self.data.loadData(filepath)
            self.data.draw()

            if self.cAdjust is not None:
                try:
                    self.cAdjust.destroy()
                except:
                    pass
                self.cAdjust = None

    def onSSMIF(self):
        """
        Display the SSMIF antenna and FEE status codes.
        """

        if self.data is not None:
            DisplaySSMIF(self)

    def onQuit(self):
        """
        Quit station master GUI.
        """

        self.destroy()

    def onColorChange(self):
        """
        Handle color coding menu selection changes.
        """

        new_color = self.color_var.get()
        if self.data is not None and self.data.color != new_color:
            self.data.color = new_color
            self.data.draw()
            if self.cAdjust is not None:
                try:
                    self.cAdjust.destroy()
                except:
                    pass
                self.cAdjust = None

    def onAdjust(self):
        """
        Bring up the colorbar adjustment dialog window.
        """

        # Don't allow contrast adjustment for Antenna Status (color=5)
        # since it uses discrete categorical values (Bad/Suspect/Good)
        if self.data is not None and self.data.color != 5:
            ContrastAdjust(self)

    def onAntenna(self):
        """
        Display meta-data for the selected antenna pair.  This includes:
        * ID numbers
        * polarizations
        * resonance point estimated from polynomial fits
        * DP1 board number
        * DP1 digitizer number
        * status code
        """

        if self.data is None:
            return
        if self.data.bestX == -1:
            messagebox.showinfo("Antenna Details", "Please select a stand first.")
            return

        ant1 = self.data.antennas[self.data.bestX-1]
        ant2 = self.data.antennas[self.data.bestY-1]

        toCompare = np.where( (self.data.freq>31e6) & (self.data.freq<70e6) )[0]

        i = self.data.bestX-1
        bestOrder = 0
        bestRMS = 1e34
        for j in range(3, 12):
            coeff = np.polyfit(self.data.freq[toCompare]/1e6, to_dB(self.data.spec[i,toCompare]), j)
            fit = np.polyval(coeff, self.data.freq[toCompare]/1e6)
            rms = ((fit - to_dB(self.data.spec[i,toCompare]))**2).sum()
            if rms < bestRMS:
                bestOrder = j
                bestRMS = rms

        coeff = np.polyfit(self.data.freq[toCompare]/1e6, to_dB(self.data.spec[i,toCompare]), bestOrder)
        fit = np.polyval(coeff, self.data.freq[toCompare]/1e6)
        res1 = self.data.freq[toCompare[np.where( fit == fit.max() )[0]]] / 1e6
        if len(res1) < 1:
            res1 = 0.0

        i = self.data.bestY-1
        bestOrder = 0
        bestRMS = 1e34
        for j in range(3, 12):
            coeff = np.polyfit(self.data.freq[toCompare]/1e6, to_dB(self.data.spec[i,toCompare]), j)
            fit = np.polyval(coeff, self.data.freq[toCompare]/1e6)
            rms = ((fit - to_dB(self.data.spec[i,toCompare]))**2).sum()
            if rms < bestRMS:
                bestOrder = j
                bestRMS = rms

        coeff = np.polyfit(self.data.freq[toCompare]/1e6, to_dB(self.data.spec[i,toCompare]), bestOrder)
        fit = np.polyval(coeff, self.data.freq[toCompare]/1e6)
        res2 = self.data.freq[toCompare[np.where( fit == fit.max() )[0]]] / 1e6
        if len(res2) < 1:
            res2 = 0.0

        outString = """Antenna: %i
Polarization: %i

Est. Resonance: %.3f MHz

DP1 Board: %i
Digitizer: %i

ARX Board: %s
Channel:   %i

Status: %i

---

Antenna: %i
Polarization: %i

Est. Resonance: %.3f MHz

DP1 Board: %i
Digitizer: %i

ARX Board: %s
Channel:   %i

Status: %i
""" % (ant1.id, ant1.pol, res1, ant1.board, ant1.digitizer, ant1.arx.id, ant1.arx.channel, ant1.status,
       ant2.id, ant2.pol, res2, ant2.board, ant2.digitizer, ant2.arx.id, ant2.arx.channel, ant2.status)

        messagebox.showinfo("Antenna Details", outString)

    def onStand(self):
        """
        Display meta-data about the selected stand.  This includes:
        * ID number
        * position
        * distance from the center of the shelter
        * distance from the fence
        """

        if self.data is None:
            return
        if self.data.bestX == -1:
            messagebox.showinfo("Stand Detail", "Please select a stand first.")
            return

        std = self.data.antennas[self.data.bestX-1].stand
        if self.data.station == 'LWA1':
            shlDist = np.sqrt( (std.x - 56.965)**2 + (std.y - 86.623)**2 )
        elif self.data.station == 'LWASV':
            shlDist = np.sqrt( (std.x + 26.790)**2 + (std.y - 48.908)**2 )
        else:
            shlDist = np.nan

        if self.data.station == 'LWA1':
            fenDistA = np.zeros(4)
            k = 0
            for p1,p2 in zip([(-59.827,59.752), (59.771,59.864), (60.148,-59.618), (-59.700,-59.948)], [(59.771,59.864), (60.148,-59.618), (-59.700,-59.948), (-59.827,59.752)]):
                x1 = p1[0]
                y1 = p1[1]
                x2 = p2[0]
                y2 = p2[1]

                a = (y2-y1)/(x2-x1)
                b = (y2*x1-y1*x2)/(x2-x1)

                x3 = std.x
                y3 = std.y

                x4 = (x3/a + y3 - b)*a / (a**2+1)
                y4 = a*x4 + b

                fenDistA[k] = np.sqrt( (x3-x4)**2 + (y3-y4)**2 )
                k += 1

            # Catch things outside the fence
            if abs(std.x) > 60 or abs(std.y) > 60:
                k = 0
                for p1 in [(-59.827,59.752), (59.771,59.864), (60.148,-59.618), (-59.700,-59.948)]:
                    x1 = p1[0]
                    y1 = p1[1]

                    x3 = std.x
                    y3 = std.y

                    fenDistA[k] = np.sqrt( (x3-x1)**2 + (y3-y1)**2 )
                    k += 1

            fenDist = fenDistA.min()
        elif self.data.station == 'LWASV':
            fenDist = np.nan
        else:
            fenDist = np.nan

        outString = """Stand: %i

Coordinates:
x = %.3f m
y = %.3f m
z = %.3f m

Shelter Distance: %.3f m
Fence Distance: %.3f m
""" % (std.id, std.x, std.y, std.z, shlDist, fenDist)

        messagebox.showinfo("Stand Detail", outString)

    def onCable(self):
        """
        Display meta-data about the RPD cables connecting the the selected
        stand/antennas back to the shelter.  This includes:
        * ID names
        * lengths
        * delay at 10 MHz
        * gain at 10 MHz
        """

        if self.data is None:
            return
        if self.data.bestX == -1:
            messagebox.showinfo("Cable Details", "Please select a stand first.")
            return

        ant1 = self.data.antennas[self.data.bestX-1]
        ant2 = self.data.antennas[self.data.bestY-1]

        outString = """Antenna: %i
Cable: %s

Length: %.1f m

Delay @ 10 MHz: %.1f ns
Gain @ 10 MHz: %.1f dB

---

Antenna: %i
Cable: %s

Length: %.1f m

Delay @ 10 MHz: %.1f ns
Gain @ 10 MHz: %.1f dB
""" % (ant1.id, ant1.cable.id, ant1.cable.length, ant1.cable.delay(10e6, ns=True), to_dB(ant1.cable.gain(10e6)),
       ant2.id, ant2.cable.id, ant2.cable.length, ant2.cable.delay(10e6, ns=True), to_dB(ant2.cable.gain(10e6)))

        messagebox.showinfo("Cable Details", outString)

    def onFEE(self):
        """
        Display meta-data about the FEE installed in the selected stand.
        This includes:
        * ID name
        * which antennas are connected to which ports
        * gain settings
        * status
        """

        if self.data is None:
            return
        if self.data.bestX == -1:
            messagebox.showinfo("FEE Detail", "Please select a stand first.")
            return

        fee = self.data.antennas[self.data.bestX-1].fee
        portXA = self.data.antennas[self.data.bestX-1].id
        portXP = self.data.antennas[self.data.bestX-1].fee_port
        portYA = self.data.antennas[self.data.bestY-1].id
        portYP = self.data.antennas[self.data.bestY-1].fee_port

        outString = """FEE: %s

Ports:
%i = antenna %i (X pol.)
%i = antenna %i (Y pol.)

Gains:
1 = %.3f dB
2 = %.3f dB

Status: %i
""" % (fee.id, portXP, portXA, portYP, portYA, fee.gain1, fee.gain2, fee.status)

        messagebox.showinfo("FEE Detail", outString)

    def onRFI(self):
        if self.data is None:
            return
        if self.data.bestX == -1:
            messagebox.showinfo("Shelter RFI Details", "Please select a stand first.")
            return

        ant1 = self.data.antennas[self.data.bestX-1]
        ant2 = self.data.antennas[self.data.bestY-1]

        rfi1 = np.where( (self.data.freq>45e6) & (self.data.freq<47e6) )[0]
        rfi2 = np.where( (self.data.freq>63e6) & (self.data.freq<65e6) )[0]
        corr = np.where( (self.data.freq>75e6) & (self.data.freq<77e6) )[0]

        a1r1 = (self.data.spec[self.data.bestX-1,rfi1] / self.data.specTemplate[rfi1]).max()
        a1r2 = (self.data.spec[self.data.bestX-1,rfi2] / self.data.specTemplate[rfi2]).max()
        c1 = (self.data.spec[self.data.bestX-1,corr] / self.data.specTemplate[corr]).mean()
        a2r1 = (self.data.spec[self.data.bestY-1,rfi1] / self.data.specTemplate[rfi1]).max()
        a2r2 = (self.data.spec[self.data.bestY-1,rfi2] / self.data.specTemplate[rfi2]).max()
        c2 = (self.data.spec[self.data.bestY-1,corr] / self.data.specTemplate[corr]).mean()

        outString = """Antenna: %i
Polarization: %i

RFI-46 Index:
raw = %.3f
corrected = %.3f

RFI-64 Index:
raw = %.3f
corrected = %.3f

---

Antenna: %i
Polarization: %i

RFI-46 Index:
raw = %.3f
corrected = %.3f

RFI-64 Index:
raw = %.3f
corrected = %.3f
""" % (ant1.id, ant1.pol, a1r1, a1r1/c1, a1r2, a1r2/c1,
       ant2.id, ant2.pol, a2r1, a2r1/c2, a2r2, a2r2/c2)

        messagebox.showinfo("Shelter RFI Details", outString)

    def onStatus(self):
        """
        Display the change antenna/FEE status dialog.
        """

        if self.data is not None and self.data.bestX > 0:
            StatusChangeDialog(self)


    def onHistogram(self):
        """"
        Display the ADC histogram plots.
        """

        if self.data is None:
            return
        if self.data.adcHistogram is None:
            messagebox.showinfo("ADC Histogram", "ADC histogram data is not available in this file.")
            return
        if self.data.bestX < 1:
            messagebox.showinfo("ADC Histogram", "Please select a stand first.")
            return
        ADCHistogramDisplay(self)

    def onavgPower(self):
        """
        Display the average power plots.
        """

        if self.data is None:
            return
        if self.data.avgPower is None:
            messagebox.showinfo("Average Power", "Average power data is not available in this file.")
            return
        if self.data.bestX < 1:
            messagebox.showinfo("Average Power", "Please select a stand first.")
            return
        AvgPowerDisplay(self)

    def onDataRange(self):
        """
        Display the data range plots.
        """

        if self.data is None:
            return
        if self.data.dataRange is None:
            messagebox.showinfo("Data Range", "Data range information is not available in this file.")
            return
        if self.data.bestX < 1:
            messagebox.showinfo("Data Range", "Please select a stand first.")
            return
        DataRangeDisplay(self)

    def onavgPowerSummary(self):
        """
        Display a message box with the average power summary.
        """

        if self.data is None:
            return
        if self.data.avgPower is None:
            messagebox.showinfo("Average Power Summary", "Average power data is not available in this file.")
            return
        if self.data.bestX == -1:
            messagebox.showinfo("Average Power Summary", "Please select a stand first.")
            return

        ant1 = self.data.antennas[self.data.bestX-1]
        dat1 = self.data.avgPower[self.data.bestX-1,:]
        ant2 = self.data.antennas[self.data.bestY-1]
        dat2 = self.data.avgPower[self.data.bestY-1,:]

        outString = """Antenna: %i,  Polarization: %i
Global Mean: %.2f +/- %.2f
Global Range: %.2f to %.2f

Antenna: %i,  Polarization: %i
Global Mean: %.2f +/- %.2f
Global Range: %.2f to %.2f""" % (ant1.id, ant1.pol, dat1.mean(), dat1.std(), dat1.min(), dat1.max(),
                                 ant2.id, ant2.pol, dat2.mean(), dat2.std(), dat2.min(), dat2.max())

        messagebox.showinfo("Average Power Details", outString)

    def onSelectAntenna(self):
        """
        Bring up a dialog box to find an antenna based on its ID number.
        """

        if self.data is None:
            return

        try:
            currAnt = self.data.antennas[self.data.bestX-1].id
        except:
            currAnt = 1

        box = SelectBox(self, mode='antenna', current=currAnt)
        self.wait_window(box)
        if box.result is not None:
            try:
                antID = int(box.result)
                if antID < 1 or antID > 520:
                    return
                if self.data.antennas is None:
                    return
                for ant in self.data.antennas:
                    if ant.id == antID:
                        self.data.drawSpectrum(ant.stand.x, ant.stand.y)
                        self.data.makeMark(ant.stand.x, ant.stand.y)
                        break
            except ValueError:
                pass

    def onSelectStand(self):
        """
        Bring up a dialog box to find a stand based on its ID number.
        """

        if self.data is None:
            return

        try:
            currStand = self.data.antennas[self.data.bestX-1].stand.id
        except:
            currStand = 1

        box = SelectBox(self, mode='stand', current=currStand)
        self.wait_window(box)
        if box.result is not None:
            try:
                stdID = int(box.result)
                if stdID < 1 or stdID > 260:
                    return
                if self.data.antennas is None:
                    return
                for ant in self.data.antennas:
                    if ant.stand.id == stdID:
                        self.data.drawSpectrum(ant.stand.x, ant.stand.y, preferStand=stdID)
                        self.data.makeMark(ant.stand.x, ant.stand.y)
                        break
            except ValueError:
                pass

    def onSelectDigitizer(self):
        """
        Bring up a dialog box to find a antenna associated with a particular
        digitizer number.
        """

        if self.data is None:
            return

        try:
            currDig = self.data.antennas[self.data.bestX-1].digitizer
        except:
            currDig = 1

        box = SelectBox(self, mode='digitizer', current=currDig)
        self.wait_window(box)
        if box.result is not None:
            try:
                digID = int(box.result)
                if digID < 1 or digID > 520:
                    return
                if self.data.antennas is None:
                    return
                for ant in self.data.antennas:
                    if ant.digitizer == digID:
                        self.data.drawSpectrum(ant.stand.x, ant.stand.y)
                        self.data.makeMark(ant.stand.x, ant.stand.y)
                        break
            except ValueError:
                pass

    def onKeyPress(self, event):
        """
        Use the arrow keys to move around the array.
        """

        if self.data is None:
            return

        keysym = event.keysym
        try:
            currStand = self.data.antennas[self.data.bestX-1].stand.id
        except:
            currStand = 1

        if keysym == 'Left' and currStand > 1:
            ## Left decreases the stand number by 1
            currStand -= 1
        elif keysym == 'Right' and currStand < len(self.data.antennas):
            ## Right increases the stand number by 1
            currStand += 1
        else:
            return

        for ant in self.data.antennas:
            if ant.stand.id == currStand:
                self.data.drawSpectrum(ant.stand.x, ant.stand.y, preferStand=currStand)
                self.data.makeMark(ant.stand.x, ant.stand.y)
                break

    def resizePlots(self):
        """Handle window resize events.

        Note: With tk's pack geometry manager (fill=BOTH, expand=True),
        the canvas widgets resize automatically. We just need to redraw
        to update the display, not manually set figure sizes.
        """
        if self.data is None:
            return

        try:
            self.canvas1.draw_idle()
            self.canvas2.draw_idle()
        except:
            pass

    def GetToolBar(self):
        # You will need to override GetToolBar if you are using an
        # unmanaged toolbar in your frame
        return self.toolbar


class ContrastAdjust(tk.Toplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title('Contrast Adjustment')
        self.transient(parent)
        self.resizable(False, False)

        self.parent = parent

        self.initUI()
        self.initEvents()

        self.parent.cAdjust = self

    def initUI(self):
        panel = ttk.Frame(self, padding="15")
        panel.pack(fill=tk.BOTH, expand=True)

        mode_names = {
            0: 'Median Comparison',
            1: 'RFI-46 Index',
            2: 'RFI-64 Index',
            3: 'RFI-76 Index',
            4: 'Wiggle Index',
            5: 'Antenna Status',
            6: 'Resonance Point'
        }
        mode = mode_names.get(self.parent.data.color, 'Unknown')

        ttk.Label(panel, text='Color Coding Mode: %s' % mode).grid(
            row=0, column=0, columnspan=4, sticky='w', padx=5, pady=5)

        # Upper limit
        ttk.Label(panel, text='Upper Limit:').grid(row=1, column=0, sticky='w', padx=5)
        self.upr_var = tk.StringVar()
        self.uprText = ttk.Entry(panel, textvariable=self.upr_var, state='readonly', width=10)
        self.uprText.grid(row=1, column=1, padx=5)
        ttk.Button(panel, text='-', width=3, command=self.onUpperDecrease).grid(row=1, column=2, padx=2)
        ttk.Button(panel, text='+', width=3, command=self.onUpperIncrease).grid(row=1, column=3, padx=2)

        # Lower limit
        ttk.Label(panel, text='Lower Limit:').grid(row=2, column=0, sticky='w', padx=5)
        self.lwr_var = tk.StringVar()
        self.lwrText = ttk.Entry(panel, textvariable=self.lwr_var, state='readonly', width=10)
        self.lwrText.grid(row=2, column=1, padx=5)
        ttk.Button(panel, text='-', width=3, command=self.onLowerDecrease).grid(row=2, column=2, padx=2)
        ttk.Button(panel, text='+', width=3, command=self.onLowerIncrease).grid(row=2, column=3, padx=2)

        # Range
        ttk.Label(panel, text='Range:').grid(row=3, column=0, sticky='w', padx=5)
        self.rng_var = tk.StringVar()
        self.rngText = ttk.Entry(panel, textvariable=self.rng_var, state='readonly', width=10)
        self.rngText.grid(row=3, column=1, padx=5)

        # Separator
        ttk.Separator(panel, orient=tk.HORIZONTAL).grid(
            row=4, column=0, columnspan=4, sticky='ew', pady=10)

        # OK button
        ttk.Button(panel, text='Ok', command=self.onOk).grid(
            row=5, column=3, sticky='e', padx=5, pady=5)

        # Set current values
        color = self.parent.data.color
        self.upr_var.set('%.1f' % self.parent.data.limits[color][1])
        self.lwr_var.set('%.1f' % self.parent.data.limits[color][0])
        self.rng_var.set('%.1f' % self._getRange(color))

    def initEvents(self):
        self.protocol("WM_DELETE_WINDOW", self.onOk)

    def _redraw(self):
        """Redraw the main plot and force update."""
        self.parent.data.draw()

    def onUpperDecrease(self):
        color = self.parent.data.color
        self.parent.data.limits[color][1] -= self._getIncrement(color)
        self.upr_var.set('%.1f' % self.parent.data.limits[color][1])
        self.rng_var.set('%.1f' % self._getRange(color))
        self._redraw()

    def onUpperIncrease(self):
        color = self.parent.data.color
        self.parent.data.limits[color][1] += self._getIncrement(color)
        self.upr_var.set('%.1f' % self.parent.data.limits[color][1])
        self.rng_var.set('%.1f' % self._getRange(color))
        self._redraw()

    def onLowerDecrease(self):
        color = self.parent.data.color
        self.parent.data.limits[color][0] -= self._getIncrement(color)
        self.lwr_var.set('%.1f' % self.parent.data.limits[color][0])
        self.rng_var.set('%.1f' % self._getRange(color))
        self._redraw()

    def onLowerIncrease(self):
        color = self.parent.data.color
        self.parent.data.limits[color][0] += self._getIncrement(color)
        self.lwr_var.set('%.1f' % self.parent.data.limits[color][0])
        self.rng_var.set('%.1f' % self._getRange(color))
        self._redraw()

    def onOk(self):
        self.parent.cAdjust = None
        self.destroy()

    def _getRange(self, color):
        return (self.parent.data.limits[color][1] - self.parent.data.limits[color][0])

    def _getIncrement(self, color):
        return 0.1*self._getRange(color)


class AvgPowerDisplay(tk.Toplevel):
    """
    Window for displaying the average power with time for the selected stand.
    """

    def __init__(self, parent):
        super().__init__(parent)
        self.title('Time-Averaged Power')
        self.geometry('500x400')
        self.transient(parent)

        self.parent = parent

        self.initUI()
        self.initPlot()

    def _nextTen(self, value):
        """
        Round a positive value to the next highest multiple of ten.
        """

        return 10*np.ceil(value/10.0)

    def initUI(self):
        """
        Start the user interface.
        """

        # Status bar
        self.statusbar = ttk.Label(self, text="", relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(fill=tk.X, side=tk.BOTTOM)

        # Toolbar frame
        toolbar_frame = ttk.Frame(self)
        toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)

        # Plot frame
        plot_frame = ttk.Frame(self)
        plot_frame.pack(fill=tk.BOTH, expand=True)

        self.figure = Figure(tight_layout=True)
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()

    def initPlot(self):
        """
        Populate the figure/canvas areas with a plot.  We only need to do this
        once for this type of window.
        """

        avgPower = self.parent.data.avgPower
        bestX = self.parent.data.bestX
        bestY = self.parent.data.bestY

        if avgPower is None:
            return False
        if bestX < 1:
            return False

        self.figure.clf()
        self.ax1 = self.figure.gca()

        ant1 = self.parent.data.antennas[bestX-1]
        ant2 = self.parent.data.antennas[bestY-1]

        # Average power plot
        tScale = float(round(avgPower.shape[1] / 61.2244898))
        t = np.arange(0,avgPower.shape[1])/tScale + 0.5/tScale
        self.ax1.errorbar(t, avgPower[bestX-1,:], xerr=0.5/tScale, linestyle=' ', marker='+', label='Pol. %i' % ant1.pol, capsize=0)
        self.ax1.errorbar(t, avgPower[bestY-1,:], xerr=0.5/tScale, linestyle=' ', marker='+', label='Pol. %i' % ant2.pol, capsize=0)

        # Set ranges
        self.ax1.set_xlim([0, 61])
        self.ax1.set_ylim([0, self._nextTen(avgPower.max())])

        # Labels
        self.ax1.set_title('Stand #%i' % ant1.stand.id)
        self.ax1.set_xlabel('Time [ms]')
        self.ax1.set_ylabel('Average Power [counts]')

        # Legend
        self.ax1.legend(loc=0)
        self.tScale = tScale

        ## Draw and save the click (Why?)
        self.canvas.draw()
        self.connect()

    def connect(self):
        """
        Connect to all the events we need to interact with the plots.
        """

        self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_motion(self, event):
        """
        Deal with motion events in the stand field window.  This involves
        setting the status bar with the current x and y coordinates as well
        as the stand number of the selected stand (if any).
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            try:
                tScale = self.tScale
                ap1 = self.parent.data.avgPower[self.parent.data.bestX-1,int(clickX*tScale)]
                ap2 = self.parent.data.avgPower[self.parent.data.bestY-1,int(clickX*tScale)]
                self.statusbar.config(text="t=%.2f ms, X pol. Power=%.2f counts, Y pol. Power=%.2f" % (clickX, ap1, ap2))
            except IndexError:
                self.statusbar.config(text="")
        else:
            self.statusbar.config(text="")

    def disconnect(self):
        """
        Disconnect all the stored connection ids.
        """

        self.figure.canvas.mpl_disconnect(self.cidmotion)

    def GetToolBar(self):
        return self.toolbar


class DataRangeDisplay(tk.Toplevel):
    """
    Window for displaying the time series mean, min, and maximum raw data
    values.
    """

    def __init__(self, parent):
        super().__init__(parent)
        self.title('Range of Raw Data')
        self.geometry('500x400')
        self.transient(parent)

        self.parent = parent

        self.initUI()
        self.initPlot()

    def _nextTen(self, value):
        """
        Round a positive value to the next highest multiple of ten.
        """

        return 10*np.ceil(value/10.0)

    def initUI(self):
        """
        Start the user interface.
        """

        # Status bar
        self.statusbar = ttk.Label(self, text="", relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(fill=tk.X, side=tk.BOTTOM)

        # Toolbar frame
        toolbar_frame = ttk.Frame(self)
        toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)

        # Plot frame
        plot_frame = ttk.Frame(self)
        plot_frame.pack(fill=tk.BOTH, expand=True)

        self.figure = Figure(tight_layout=True)
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()

    def initPlot(self):
        """
        Populate the figure/canvas areas with a plot.  We only need to do this
        once for this type of window.
        """

        dataRange = self.parent.data.dataRange
        bestX = self.parent.data.bestX
        bestY = self.parent.data.bestY

        if dataRange is None:
            return False
        if bestX < 1:
            return False

        self.figure.clf()
        self.ax1 = self.figure.gca()

        ant1 = self.parent.data.antennas[bestX-1]
        ant2 = self.parent.data.antennas[bestY-1]

        # Data Range
        tScale = float(round(dataRange.shape[1] / 61.2244898))
        t = np.arange(0,dataRange.shape[1])/tScale + 0.5/tScale
        eb1 = np.zeros((2,t.size))
        eb1[0,:] = dataRange[bestX-1,:,1] - dataRange[bestX-1,:,0]
        eb1[1,:] = dataRange[bestX-1,:,2] - dataRange[bestX-1,:,1]
        self.ax1.errorbar(t, dataRange[bestX-1,:,1], xerr=0.5/tScale, yerr=eb1, linestyle=' ', marker='+', label='Pol. %i' % ant1.pol)
        eb2 = np.zeros((2,t.size))
        eb2[0,:] = dataRange[bestY-1,:,1] - dataRange[bestY-1,:,0]
        eb2[1,:] = dataRange[bestY-1,:,2] - dataRange[bestY-1,:,1]
        self.ax1.errorbar(t, dataRange[bestY-1,:,1], xerr=0.5/tScale, yerr=eb2, linestyle=' ', marker='+', label='Pol. %i' % ant2.pol)

        # Set ranges
        self.ax1.set_xlim([0, 61])
        self.ax1.set_ylim([-2048, 2047])

        # Labels
        self.ax1.set_title('Stand #%i' % ant1.stand.id)
        self.ax1.set_xlabel('Time [ms]')
        self.ax1.set_ylabel('Data Range [counts]')

        # Legend
        self.ax1.legend(loc=0)
        self.tScale = tScale

        ## Draw and save the click (Why?)
        self.canvas.draw()
        self.connect()

    def connect(self):
        """
        Connect to all the events we need to interact with the plots.
        """

        self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_motion(self, event):
        """
        Deal with motion events in the stand field window.  This involves
        setting the status bar with the current x and y coordinates as well
        as the stand number of the selected stand (if any).
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            try:
                tScale = self.tScale
                dr1 = self.parent.data.dataRange[self.parent.data.bestX-1,int(clickX*tScale),:]
                dr2 = self.parent.data.dataRange[self.parent.data.bestY-1,int(clickX*tScale),:]
                self.statusbar.config(text="t=%.2f ms, X pol. Range: %+i to %+i, Y Pol. Range: %+i to %+i counts" % (clickX, dr1[0], dr1[2], dr2[0], dr2[2]))
            except IndexError:
                self.statusbar.config(text="")
        else:
            self.statusbar.config(text="")

    def disconnect(self):
        """
        Disconnect all the stored connection ids.
        """

        self.figure.canvas.mpl_disconnect(self.cidmotion)

    def GetToolBar(self):
        return self.toolbar


class ADCHistogramDisplay(tk.Toplevel):
    """
    Window for displaying the ADC histogram for the selected stand.
    """

    def __init__(self, parent):
        super().__init__(parent)
        self.title('ADC Histogram')
        self.geometry('500x400')
        self.transient(parent)

        self.parent = parent

        self.initUI()
        self.initPlot()

    def _nextThousand(self, value):
        """
        Round a positive value to the next highest multiple of a thousand.
        """

        return 1000*np.ceil(value/1000.0)

    def initUI(self):
        """
        Start the user interface.
        """

        # Status bar
        self.statusbar = ttk.Label(self, text="", relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(fill=tk.X, side=tk.BOTTOM)

        # Toolbar frame
        toolbar_frame = ttk.Frame(self)
        toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)

        # Plot frame
        plot_frame = ttk.Frame(self)
        plot_frame.pack(fill=tk.BOTH, expand=True)

        self.figure = Figure(tight_layout=True)
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()

    def initPlot(self):
        """
        Populate the figure/canvas areas with a plot.  We only need to do this
        once for this type of window.
        """

        adcHistogram = self.parent.data.adcHistogram
        bestX = self.parent.data.bestX
        bestY = self.parent.data.bestY

        if adcHistogram is None:
            return False
        if bestX < 1:
            return False

        self.figure.clf()
        self.ax1 = self.figure.gca()

        ant1 = self.parent.data.antennas[bestX-1]
        ant2 = self.parent.data.antennas[bestY-1]

        # Histogram plot
        histBins = list(range(-2048, 2049))
        left, right = histBins[:-1], histBins[1:]
        v = np.array([left,right]).T.flatten()
        hX = np.array([adcHistogram[bestX-1,:], adcHistogram[bestX-1,:]]).T.flatten()
        hY = np.array([adcHistogram[bestY-1,:], adcHistogram[bestY-1,:]]).T.flatten()

        self.ax1.plot(v, hX, label='Pol. %i' % ant1.pol)
        self.ax1.plot(v, hY, label='Pol. %i' % ant2.pol)

        # Calculate and display the RMS
        rmsX = np.sqrt( (np.array(histBins[:-1])**2 * adcHistogram[bestX-1,:]).sum() / adcHistogram[bestX-1,:].sum() )
        rmsY = np.sqrt( (np.array(histBins[:-1])**2 * adcHistogram[bestY-1,:]).sum() / adcHistogram[bestY-1,:].sum() )
        self.ax1.text(0.08, 0.90, 'RMS$_%i$=%.1f' % (ant1.pol, rmsX), transform=self.ax1.transAxes)
        self.ax1.text(0.08, 0.85, 'RMS$_%i$=%.1f' % (ant2.pol, rmsY), transform=self.ax1.transAxes)

        # Set ranges
        self.ax1.set_xlim([-2048, 2047])
        hMax = max([self._nextThousand(adcHistogram[bestX-1,:].max()), self._nextThousand(adcHistogram[bestY-1,:].max())])
        hMax += 5000
        self.ax1.set_ylim([0, hMax])

        # Labels
        self.ax1.set_title('Stand #%i' % ant1.stand.id)
        self.ax1.set_xlabel('ADC Value')
        self.ax1.set_ylabel('Number')

        # Legend
        self.ax1.legend(loc=0)
        self.histBins = histBins

        ## Draw and save the click (Why?)
        self.canvas.draw()
        self.connect()

    def connect(self):
        """
        Connect to all the events we need to interact with the plots.
        """

        self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_motion(self, event):
        """
        Deal with motion events in the stand field window.  This involves
        setting the status bar with the current x and y coordinates as well
        as the stand number of the selected stand (if any).
        """

        if event.inaxes:
            clickX = int(event.xdata)
            clickY = event.ydata

            try:
                idx = self.histBins.index(clickX)

                ap1 = self.parent.data.adcHistogram[self.parent.data.bestX-1,idx]
                ap2 = self.parent.data.adcHistogram[self.parent.data.bestY-1,idx]
                self.statusbar.config(text="v=%i, X pol. Count=%i counts, Y pol. Count=%i" % (clickX, ap1, ap2))
            except (IndexError, ValueError):
                self.statusbar.config(text="")
        else:
            self.statusbar.config(text="")

    def disconnect(self):
        """
        Disconnect all the stored connection ids.
        """

        self.figure.canvas.mpl_disconnect(self.cidmotion)

    def GetToolBar(self):
        return self.toolbar


class SelectBox(tk.Toplevel):
    """
    Window for displaying the a simple dialog to find an antenna/stand/digitizer.
    """

    def __init__(self, parent, mode='antenna', current=None):
        super().__init__(parent)
        self.title('Find %s by ID' % mode.capitalize())
        self.geometry('200x125')
        self.transient(parent)
        self.grab_set()

        self.parent = parent
        self.mode = mode
        self.current = current
        self.result = None

        self.initUI()
        self.initEvents()

    def initUI(self):
        """
        Start the user interface.
        """

        # Main frame with label frame
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        lf = ttk.LabelFrame(main_frame, text='%s ID' % self.mode.capitalize())
        lf.pack(fill=tk.X, pady=5)

        self.input_var = tk.StringVar()
        self.input = ttk.Entry(lf, textvariable=self.input_var, width=20)
        self.input.pack(padx=10, pady=5)

        if self.mode == 'stand':
            ttk.Label(lf, text='Limits: 1 - 260, inclusive').pack(padx=10, pady=2)
        else:
            ttk.Label(lf, text='Limits: 1 - 520, inclusive').pack(padx=10, pady=2)

        # Button frame
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill=tk.X, pady=10)

        ttk.Button(btn_frame, text='Ok', command=self.onOK).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text='Close', command=self.onCancel).pack(side=tk.LEFT, padx=5)

        if self.current is not None:
            self.input_var.set('%i' % self.current)

        self.input.focus_set()

    def initEvents(self):
        """
        Bind the events needed to make this run.
        """

        self.bind('<Return>', lambda e: self.onOK())
        self.bind('<Escape>', lambda e: self.onCancel())
        self.bind('<Up>', self.onKeyUp)
        self.bind('<Down>', self.onKeyDown)
        self.protocol("WM_DELETE_WINDOW", self.onCancel)

    def onKeyUp(self, event=None):
        try:
            value = int(self.input_var.get(), 10)
            maxValue = 260 if self.mode == 'stand' else 520
            value = min([maxValue, value+1])
            self.input_var.set('%i' % value)
        except ValueError:
            pass

    def onKeyDown(self, event=None):
        try:
            value = int(self.input_var.get(), 10)
            value = max([1, value-1])
            self.input_var.set('%i' % value)
        except ValueError:
            pass

    def onOK(self):
        self.result = self.input_var.get()
        self.destroy()

    def onCancel(self):
        self.result = None
        self.destroy()


class StatusChangeDialog(tk.Toplevel):
    """
    Window for changing the status of an antenna or FEE.
    """

    def __init__(self, parent):
        super().__init__(parent)
        self.title('Change Antenna/FEE Status')
        self.geometry('450x225')
        self.transient(parent)

        self.parent = parent

        self.initUI()
        self.initEvents()

    def initUI(self):
        bestX = self.parent.data.bestX
        bestY = self.parent.data.bestY

        ant1 = self.parent.data.antennas[bestX-1]
        ant2 = self.parent.data.antennas[bestY-1]

        panel = ttk.Frame(self, padding="10")
        panel.pack(fill=tk.BOTH, expand=True)

        row = 0

        # Stand label
        ttk.Label(panel, text='Stand #%i' % ant1.stand.id, font=('TkDefaultFont', 10, 'bold')).grid(
            row=row, column=0, columnspan=5, sticky='w', pady=5)
        row += 1

        # X Pol section
        ttk.Label(panel, text='X Pol.').grid(row=row, column=0, columnspan=5, sticky='w')
        row += 1

        ttk.Label(panel, text='Antenna Status').grid(row=row, column=0, sticky='w')
        self.ant1_status = tk.IntVar(value=ant1.status)
        ttk.Radiobutton(panel, text='Not Installed', variable=self.ant1_status, value=0).grid(row=row, column=1)
        ttk.Radiobutton(panel, text='Bad', variable=self.ant1_status, value=1).grid(row=row, column=2)
        ttk.Radiobutton(panel, text='Suspect', variable=self.ant1_status, value=2).grid(row=row, column=3)
        ttk.Radiobutton(panel, text='Good', variable=self.ant1_status, value=3).grid(row=row, column=4)
        row += 1

        # Y Pol section
        ttk.Label(panel, text='Y Pol.').grid(row=row, column=0, columnspan=5, sticky='w')
        row += 1

        ttk.Label(panel, text='Antenna Status').grid(row=row, column=0, sticky='w')
        self.ant2_status = tk.IntVar(value=ant2.status)
        ttk.Radiobutton(panel, text='Not Installed', variable=self.ant2_status, value=0).grid(row=row, column=1)
        ttk.Radiobutton(panel, text='Bad', variable=self.ant2_status, value=1).grid(row=row, column=2)
        ttk.Radiobutton(panel, text='Suspect', variable=self.ant2_status, value=2).grid(row=row, column=3)
        ttk.Radiobutton(panel, text='Good', variable=self.ant2_status, value=3).grid(row=row, column=4)
        row += 1

        # Separator
        ttk.Separator(panel, orient=tk.HORIZONTAL).grid(row=row, column=0, columnspan=5, sticky='ew', pady=10)
        row += 1

        # FEE Status
        ttk.Label(panel, text='FEE Status').grid(row=row, column=0, sticky='w')
        self.fee_status = tk.IntVar(value=ant1.fee.status)
        ttk.Radiobutton(panel, text='Not Installed', variable=self.fee_status, value=0).grid(row=row, column=1)
        ttk.Radiobutton(panel, text='Bad', variable=self.fee_status, value=1).grid(row=row, column=2)
        ttk.Radiobutton(panel, text='Suspect', variable=self.fee_status, value=2).grid(row=row, column=3)
        ttk.Radiobutton(panel, text='Good', variable=self.fee_status, value=3).grid(row=row, column=4)
        row += 1

        # Separator
        ttk.Separator(panel, orient=tk.HORIZONTAL).grid(row=row, column=0, columnspan=5, sticky='ew', pady=10)
        row += 1

        # Buttons
        btn_frame = ttk.Frame(panel)
        btn_frame.grid(row=row, column=0, columnspan=5, sticky='e')
        ttk.Button(btn_frame, text='Ok', command=self.onOk).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text='Cancel', command=self.onCancel).pack(side=tk.LEFT, padx=5)

    def initEvents(self):
        self.protocol("WM_DELETE_WINDOW", self.onCancel)

    def onOk(self):
        # Antenna 1
        self.parent.data.antennas[self.parent.data.bestX-1].status = self.ant1_status.get()

        # Antenna 2
        self.parent.data.antennas[self.parent.data.bestY-1].status = self.ant2_status.get()

        # FEE
        fee_stat = self.fee_status.get()
        self.parent.data.antennas[self.parent.data.bestX-1].fee.status = fee_stat
        self.parent.data.antennas[self.parent.data.bestY-1].fee.status = fee_stat

        # Refresh if we are in the antenna status color coding
        if self.parent.data.color == 5:
            self.parent.data.draw()

        self.destroy()

    def onCancel(self):
        self.destroy()


class DisplaySSMIF(tk.Toplevel):
    """
    Text display window for printing out the new SSMIF entries for antenna status.
    FEE status is currently not supported because of a limitation in the LSL SSMIF
    parser.
    """

    def __init__(self, parent):
        super().__init__(parent)
        self.title('SSMIF Status Codes')
        self.geometry('600x600')
        self.transient(parent)

        self.parent = parent

        self.initUI()
        self.initEvents()
        self.generateText()

    def initUI(self):
        main_frame = ttk.Frame(self, padding="5")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Text widget with scrollbar
        text_frame = ttk.Frame(main_frame)
        text_frame.pack(fill=tk.BOTH, expand=True)

        scrollbar = ttk.Scrollbar(text_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.textCtrl = tk.Text(text_frame, wrap=tk.WORD, yscrollcommand=scrollbar.set)
        self.textCtrl.pack(fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.textCtrl.yview)

        # OK button
        ttk.Button(main_frame, text='Ok', command=self.onOk).pack(side=tk.RIGHT, pady=5)

    def initEvents(self):
        self.protocol("WM_DELETE_WINDOW", self.onOk)

    def generateText(self):
        # Sort antennas by antenna number
        ants = sorted(self.parent.data.antennas, key=lambda x: x.id)

        #
        # Antenna status codes
        #

        self.textCtrl.insert(tk.END, '# -----------------------------\n# --- Antenna Status ---\n# -----------------------------\n# Status codes 0-3 summarized defined at end of this document (and in MCS0031)\n# This refers to the *antenna*, not the FEE or some combination of the two.\n# This will be set to 3 ("OK") for any antenna n <= 2*N_STD not identified.\n# *** ANT_STAT[antenna_id] goes here:\n')
        for ant in ants:
            if ant.status == 3:
                continue

            if ant.id < 10:
                self.textCtrl.insert(tk.END, 'ANT_STAT[%i]    %i\n' % (ant.id, ant.status))
            elif ant.id < 100:
                self.textCtrl.insert(tk.END, 'ANT_STAT[%i]   %i\n' % (ant.id, ant.status))
            else:
                self.textCtrl.insert(tk.END, 'ANT_STAT[%i]  %i\n' % (ant.id, ant.status))
        self.textCtrl.insert(tk.END, '\n\n')

        # Make read-only and scroll to top
        self.textCtrl.config(state='disabled')
        self.textCtrl.see('1.0')

    def onOk(self):
        self.destroy()


def main(args):
    app = MainWindow()
    if args.filename is not None:
        app.data = TBW_GUI(app)
        app.data.loadData(args.filename)
        app.data.draw()
    app.mainloop()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='display NPZ data from stationMaster in an interactive GUI sort of way',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='?', default=None,
                        help='filename to display')
    args = parser.parse_args()
    main(args)
