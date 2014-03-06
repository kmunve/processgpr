'''
Created on 20.10.2010

@author: Karsten Mueller
'''
# Built-in
import os
# Third-party
import numpy as np
from numpy.fft import fft, fftfreq, ifft
from enthought.traits.api import Bool, Int, Str, Button, HasTraits, Instance, \
    Float, Enum
from enthought.traits.ui.api import View, Item, VGroup, VGrid, HSplit
from enthought.traits.ui.editors import TextEditor
from enthought.enable.api import ComponentEditor
from enthought.chaco.api import LinearMapper, BasePlotContainer, \
    OverlayPlotContainer, VPlotContainer, \
    create_line_plot, add_default_axes
from enthought.chaco.tools.api import PanTool, SimpleZoom, RangeSelection, \
    RangeSelectionOverlay
# Local imports                                    
from processgpr.viewer.gui_util import wxOpenFile


class TraceViewer(HasTraits):
    """
    Allows visualization of a single trace and spectral analysis.
    """
    # Load Panel Items
    gpr_type = Enum('HUBRA', 'CBAND', 'FIMBUL', label='GPR type')
    file = Str()
    dir = Str()
    browse = Button(label='...')
    sample_min = Int(0, label="Sample min/max")
    sample_max = Int(-1)
    trace_no = Int(0)
    trace_max = Int(-1, label="Max:")
    show_button = Button(label="Show")
    show = Bool(False)
    # Processing Panel Items
    removegate_button = Button(label="Remove Gating")
    removegate_keep = Bool(False)
    compensate_button = Button(label="Compensate Spreading")
    compensate_keep = Bool(False)
    compensate_power = Bool(False, label='Power')
    compensate_h0 = Float(label="h0")
    place_holder = Str("")
    # Spectrum Panel Items
    Nfft_button = Button(label="Get #samples")
    Nfft_edit = Int(32768, label="#FFT")
    spectrum_button = Button(label="Spectrum")
    cspectrum_button = Button(label="Complex Spectrum")
    clac_spectrum = Bool(False)
    zphi_button = Button(label="Phase-center",
                         enabled_when="calc_spectrum is True")
    freq_min_edit = Float(500.0, label="Freq. range")
    freq_max_edit = Float(3000.0)
    # Status Panel
    status_txt = Str("TraceViewer started ...")
    # Plot Panel Items
    container = Instance(BasePlotContainer)
    rangeselect = Instance(RangeSelection)
    #
    #    load_panel = VGrid(
    #                       Item(name='file'), Item(name='browse', show_label=False),
    #                       Item(name='sample_min'), Item(name='sample_max', show_label=False),
    #                       Item(name='trace_no'), Item(name='trace_max', style='readonly'),
    #                       Item(name='show_button', show_label=False),
    #                       show_border=True, label="Choose data file and range")

    load_panel = VGroup(Item(name='gpr_type'),
                        VGrid(
                            Item(name='file'),
                            Item(name='browse', show_label=False),
                            Item(name='sample_min'),
                            Item(name='sample_max', show_label=False),
                            Item(name='trace_no'),
                            Item(name='trace_max', style='readonly'),
                            Item(name='show_button', show_label=False)),
                        show_border=True, label="Choose data file and range")

    proc_panel = VGrid(
        Item(name='removegate_keep', show_label=False),
        Item(name='removegate_button', show_label=False),
        Item(name='place_holder', show_label=False,
             style='readonly'),
        Item(name='compensate_keep', show_label=False),
        Item(name='compensate_button', show_label=False),
        Item(name='compensate_h0', show_label=True),
        Item(name='compensate_power', show_label=True),
        show_border=True, columns=3, label="Processing steps"
    )

    spec_panel = VGrid(
        Item(name='Nfft_button', show_label=False),
        Item(name='Nfft_edit', show_label=True),
        Item(name='place_holder', show_label=False,
             style='readonly'),
        Item(name='spectrum_button', show_label=False),
        Item(name='place_holder', show_label=False,
             style='readonly'),
        Item(name='place_holder', show_label=False,
             style='readonly'),
        Item(name='cspectrum_button', show_label=False),
        Item(name='place_holder', show_label=False,
             style='readonly'),
        Item(name='place_holder', show_label=False,
             style='readonly'),
        Item(name='zphi_button', show_label=False),
        Item(name='freq_min_edit', show_label=True),
        Item(name='freq_max_edit', show_label=False),
        show_border=True, columns=3, label="Spectrum"
    )

    status_panel = VGroup(
        Item(name='status_txt', show_label=False, editor=TextEditor(),
             style='readonly',
             full_size=False, height=200),
        show_border=True, label="Status"
    )

    view = View(HSplit(
        VGroup(
            load_panel, proc_panel, spec_panel, status_panel
        ),
        #                      '_',
        Item(name="container", show_label=False, editor=ComponentEditor())
        #                          visible_when="show==True")),
    ),
                icon=r"../src/recources/icons/edit-32.png",
                title="Trace Viewer",
                width=1000, height=0.9,
                resizable=True)

    def __init__(self, file="", dir=os.getcwd()):
        self.dir = dir
        self.file = file

    def _create_traceplot(self, y=None):
        x = self.x
        if y == None:
            y = self.y
        traceplot = create_line_plot((x, y), index_bounds=None, value_bounds=None,
                                     orientation='v', color='blue', width=1.0, dash='solid',
                                     value_mapper_class=LinearMapper,
                                     bgcolor='transparent', border_visible=True,
                                     add_grid=False, add_axis=False, index_sort='ascending')
        add_default_axes(traceplot, orientation='flipped', vtitle='Time [ns]', htitle='Signal')
        traceplot.origin = "top left"
        traceplot.tools.append(PanTool(traceplot, drag_button="right"))
        zoom = SimpleZoom(component=traceplot, tool_mode="box", drag_button="left", always_on=True)
        traceplot.overlays.append(zoom)

        container = OverlayPlotContainer(padding=40, padding_left=60)
        container.add(traceplot)
        self.container = container

    def _create_spectrumplot(self, x, y):
        spectrumplot = create_line_plot((x, y), index_bounds=None, value_bounds=None,
                                        orientation='v', color='green', width=1.0, dash='solid',
                                        value_mapper_class=LinearMapper,
                                        bgcolor='transparent', border_visible=True,
                                        add_grid=False, add_axis=False, index_sort='ascending')
        add_default_axes(spectrumplot, orientation='flipped', vtitle='Frequency [MHz]', htitle='Amplitude')
        spectrumplot.origin = "top left"
        spectrumplot.tools.append(PanTool(spectrumplot, drag_button="right"))
        zoom = SimpleZoom(component=spectrumplot, tool_mode="box", drag_button="left", always_on=True)
        spectrumplot.overlays.append(zoom)

        container = OverlayPlotContainer(padding=40, padding_left=60)
        container.add(spectrumplot)
        self.container = container

    def _create_complexspectrumplot(self, x, y1, y2):
        amplitudeplot = create_line_plot((x, y1), index_bounds=None, value_bounds=None,
                                         orientation='h', color='green', width=1.0, dash='solid',
                                         value_mapper_class=LinearMapper,
                                         bgcolor='white', border_visible=True,
                                         add_grid=False, add_axis=False, index_sort='ascending')
        add_default_axes(amplitudeplot, orientation='normal', htitle='Frequency [MHz]', vtitle='Amplitude')

        amplitudeplot.tools.append(PanTool(amplitudeplot, drag_button="right"))
        zoom = SimpleZoom(component=amplitudeplot, tool_mode="box", drag_button="left", always_on=True)
        amplitudeplot.overlays.append(zoom)

        phaseplot = create_line_plot((x, y2), index_bounds=None,
                                     value_bounds=None,
                                     orientation='h', color='red',
                                     width=1.0, dash='solid',
                                     value_mapper_class=LinearMapper,
                                     bgcolor='white', border_visible=True,
                                     add_grid=False, add_axis=False,
                                     index_sort='ascending')
        add_default_axes(phaseplot, orientation='normal',
                         htitle='Frequency [MHz]', vtitle='Unwrapped phase')

        #        phaseplot.tools.append(PanTool(phaseplot, drag_button="right"))
        zoom = SimpleZoom(component=phaseplot, tool_mode="box",
                          drag_button="left", always_on=True,
                          #                          enter_zoom_key=KeySpec('z')
        )
        phaseplot.overlays.append(zoom)
        self.rangeselect = RangeSelection(phaseplot, left_button_selects=False)
        self.rangeselect.on_trait_change(self.phase_center, "selection_completed")
        phaseplot.active_tool = self.rangeselect
        phaseplot.overlays.append(RangeSelectionOverlay(component=phaseplot))

        container = VPlotContainer(padding=40, padding_left=60, spacing=40)
        container.add(phaseplot)
        container.add(amplitudeplot)
        self.container = container

    def _file_changed(self):
        try:
            from tables import openFile, NoSuchNodeError

            h5file = openFile(os.path.join(self.dir, self.file), mode='r')
            try:
                sample_max, trace_max = h5file.root.data.traces.shape
                self.sample_max = int(sample_max)
                self.trace_max = int(trace_max)
            except NoSuchNodeError:
                self.sample_max, self.trace_max = h5file.root.profile.traces.shape
            h5file.close()
        except IOError:
            print "The chosen file has an unknown format."

    def _browse_fired(self):
        file, dir = wxOpenFile(self.dir, multi=False)
        self.dir = dir
        self.file = file

    def _show_button_fired(self):
        gpr = eval(self.gpr_type + '()')
        gpr.readH5(os.path.join(self.dir, self.file),
                   self.trace_no, self.trace_no + 1,
                   self.sample_min, self.sample_max)
        gpr.gatefreq = 1523123.0
        self.x = gpr.t
        self.y = gpr.data[:, 0]
        self.r = gpr.r
        self.deltaT = gpr.deltaT
        self.gpr = gpr
        self._create_traceplot()
        self.status_txt = "%s\nTrace %i from file %s displayed" % (self.status_txt, self.trace_no, self.file)
        self.clac_spectrum = False

    def _removegate_button_fired(self):
        self.gpr.removeGating()
        if self.removegate_keep:
            self.y = self.gpr.data[:, 0]
            self._create_traceplot()
            self.status_txt = "%s\nGating gain removed permanently" % (self.status_txt)
        else:
            y = self.gpr.data[:, 0]
            self._create_traceplot(y)
            self.status_txt = "%s\nGating gain removed in plot" % (self.status_txt)

    def _compensate_button_fired(self):
        print self.gpr.data.shape
        print self.gpr.r.shape
        self.gpr.compensateVolSpread(h0=self.compensate_h0, ba=50.0,
                                     power=self.compensate_power)
        if self.compensate_keep:
            self.y = self.gpr.data[:, 0]
            self._create_traceplot()
            self.status_txt = "%s\nSpreading corrected permanently" % (self.status_txt)
        else:
            y = self.gpr.data[:, 0]
            self._create_traceplot(y)
            self.status_txt = "%s\nSpreading corrected in plot" % (self.status_txt)

    def _Nfft_button_fired(self):
        self.Nfft_edit = self.y.size

    def _spectrum_button_fired(self):
        Nfft = self.y.size
        freq = fftfreq(Nfft, self.deltaT)
        Y = abs(fft(self.y)) ** 2
        x = freq[0:Nfft / 2] * 1e-6
        y = Y[0:Nfft / 2]
        self._create_spectrumplot(x, y)

    def _cspectrum_button_fired(self):
        #        Nfft = self.y.size
        Nfft = self.Nfft_edit
        freq = fftfreq(Nfft, self.deltaT)
        #        c0 = 299792458.0
        #        k = 2*pi*freq[0:Nfft/2]/(c0/sqrt(1.8))
        Y = ifft(self.y, Nfft)
        f = freq[0:Nfft / 2] * 1e-6  # [MHz]
        y1 = abs(Y)[0:Nfft / 2]
        mag = y1 / y1.max()  # normalized amplitude spectrum
        pha = np.unwrap(np.angle(Y)[0:Nfft / 2])  # unwrapped phase spectrum

        self.freq = f
        self.pha = pha
        self._create_complexspectrumplot(f, mag, pha)
        self.clac_spectrum = True

    def _zphi_button_fired(self):
        c0 = 299792458.0
        smaller = self.freq > self.freq_min_edit
        larger = self.freq < self.freq_max_edit
        mask = smaller * larger
        p1 = np.polyfit(self.freq[mask] * 1e6, self.pha[mask], 1)

        self.status_txt = "%s\nPhase center = %5.2f m" % \
                          (self.status_txt, -p1[0] * c0 / (4 * np.pi))

    def _get_index(self):
        smaller = self.freq > self.rangeselect.selection[0]
        larger = self.freq < self.rangeselect.selection[1]
        ndx = smaller * larger
        return ndx

    def phase_center(self):
        c0 = 299792458.0
        if self.rangeselect.selection is not None:
            #        mask = mag>0.04
            mask = self._get_index()
            p1 = np.polyfit(self.freq[mask] * 1e6, self.pha[mask], 1)
            self.status_txt = "%s\nPhase center = %5.2f m" % \
                              (self.status_txt, -p1[0] * c0 / (4 * np.pi))


if __name__ == "__main__":
    file = 'L_KV08_stake78_pha_stk10.h5'
    dir = r'C:\hubradata\KV08'
    tp = TraceViewer(file, dir)
    tp.configure_traits()
    