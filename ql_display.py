# This is open-source software licensed under a BSD license.
# Please see the file LICENSE.txt for details.
import pdb

import glob
import os
import re
import time
from pathlib import Path
from subprocess import Popen
from typing import Any

from ginga.misc import Bunch
from ginga import GingaPlugin
from ginga.util import paths, iohelper, grc
from ginga.gw import Widgets
from ginga.AstroImage import AstroImage
from ginga.canvas.types.layer import DrawingCanvas, ConstructedCanvas
from ginga.canvas.types.basic import Polygon
import numpy as np
from pypeit import edgetrace, io, slittrace

from astropy.io import fits as pyfits


__all__ = ['QLDisplay']
_patt = re.compile(r'"([^ "]+)"')


class QLDisplay(GingaPlugin.GlobalPlugin):

    def __init__(self, fv):
        print("Initializing QLDisplay")
        super(QLDisplay, self).__init__(fv)

        keywords = [('Object', 'OBJECT'),
                    ('Date', 'DATE-OBS'),
                    ('Time UT', 'UT'),
                    ('OFNAME', 'OFNAME'),
                    ('Target', 'TARGNAME')
                    ]
        columns = [('Type', 'icon'),
                   ('Name', 'name'),
                   ('OFNAME', 'OFNAME'),
                   ('Target', 'Target'),
                   ('Object', 'Object'),]

        self.jumpinfo = []

        # setup plugin preferences
        self.settings = {
            "home_path" : paths.home,
            "scan_fits_headers" : True,
            "scan_limit" : 100,
            "keywords" : keywords,
            "columns" : columns,
            "color_alternate_rows" : True,
            "max_rows_for_col_resize" : 5000
        }

        homedir = self.settings.get('home_path', None)
        if homedir is None or not os.path.isdir(homedir):
            homedir = paths.home
        # self.do_scanfits = self.settings.get('scan_fits_headers', False)
        # self.scan_limit = self.settings.get('scan_limit', 100)
        # self.keywords = self.settings.get('keywords', keywords)
        # self.columns = self.settings.get('columns', columns)
        # self.moving_cursor = False        

        # Make icons
        icondir = self.fv.iconpath
        self.folderpb = self.fv.get_icon(icondir, 'folder.png')
        self.filepb = self.fv.get_icon(icondir, 'file.png')
        self.fitspb = self.fv.get_icon(icondir, 'fits.png')

        icons = {
            'folderpb': self.folderpb,
            'filepb': self.filepb,
            'fitspb': self.fitspb,
        }

        self.data_interface = LocalInterface(self.settings, icons, logger=self.logger)
        

    def get_channel_info(self, fitsimage):
        chname = self.fv.get_channelName(fitsimage)
        chinfo = self.fv.get_channelInfo(chname)
        return chinfo

    def build_gui(self, container):

        vbox = Widgets.VBox()
        vbox.set_margins(2, 2, 2, 2)

        # configure the table:
        color_alternate = self.settings.get('color_alternate_rows', True)
        self.treeview = Widgets.TreeView(sortable=True, selection='multiple',
                                 use_alt_row_color=color_alternate,
                                 dragable=True)

        # set header
        col = 0
        self._name_idx = 0
        for hdr, attrname in self.settings.get('columns'):
            if attrname == 'name':
                self._name_idx = col
            col += 1
        self.treeview.setup_table(self.settings.get('columns'), 1, 'name')

        vbox.add_widget(self.treeview, stretch=1)

        # Create the text entry
        self.entry = Widgets.TextEntry()
        vbox.add_widget(self.entry, stretch=0)

        btns = Widgets.HBox()
        btns.set_spacing(3)

        self.btn_close = Widgets.Button("Close")
        self.btn_close.set_tooltip("Close this plugin")
        btns.add_widget(self.btn_close, stretch=0)
        self.btn_help = Widgets.Button("Help")
        self.btn_help.set_tooltip("Show documentation for this plugin")
        btns.add_widget(self.btn_help, stretch=0)
        self.btn_refresh = Widgets.Button("Refresh")
        self.btn_refresh.set_tooltip("Refresh the file list from the directory")
        btns.add_widget(self.btn_refresh, stretch=0)
        self.btn_load = Widgets.Button("Load")
        self.btn_load.set_tooltip("Load files selected in file pane")
        btns.add_widget(self.btn_load, stretch=0)

        
        vbox.add_widget(btns, stretch=0)

        hbox = Widgets.HBox()
        # hbox.set_spacing(3)

        self.slit_list = Widgets.ComboBox()

        hbox.add_widget(self.slit_list)

        self.btn_reduce = Widgets.Button("Reduce Slit")
        self.btn_reduce.set_tooltip("Reduce the selected slit")
        hbox.add_widget(self.btn_reduce, stretch=0)
        container.add_widget(vbox, stretch=1)

        # self.canvas.ui_set_active(True)
        # self.canvas.add_callback('cursor-down', self.tess)

        container.add_widget(hbox, stretch=0)


        self.treeview.add_callback('activated', self.table_dblclick_cb)
        self.treeview.add_callback('selected', self.table_selected_cb)
        self.entry.add_callback('activated', self.text_entry_cb)
        self.btn_close.add_callback('activated', lambda w: self.close())
        self.btn_help.add_callback('activated', lambda w: self.help())
        self.btn_refresh.add_callback('activated', lambda w: self.refresh_cb())
        self.btn_load.add_callback('activated', lambda w: self.load_cb())
        self.btn_reduce.add_callback('activated', lambda w: self.reduce_slit())

    ###
    # CALLBACKS
    ###
        
    # Table:
    def table_dblclick_cb(self, widget, res_dict):
        self.logger.info("Table double clicked")
        # First, call the interface callback
        # Then, pull the hdul in and figure out how to display it
        self.data_interface.item_dblclicked_cb(widget, res_dict)

        self.draw_graphics()

    def table_selected_cb(self, widget, res_dict):
        self.logger.info("Table selected")
        # First, call the table selected callback
        # Then, pull the path of the selected file and display it in the text entry
        self.data_interface.item_selected_cb(widget, res_dict)
        file_path = self.data_interface.get_selected_file()
        self.entry.set_text(file_path)

    # Text Entry:
    def text_entry_cb(self, widget):
        self.logger.info("Text entry activated")
        self.data_interface.browse_cb(widget)
        self.draw_graphics()

    def help_cb(self):
        pass
    
    def refresh_cb(self):
        self.logger.info("Refresh button activated")
        self.data_interface.refresh()

    def load_cb(self):
        self.logger.info("Load button activated")

        path = str(self.entry.get_text()).strip()
        self.data_interface.load_cb(path, self.treeview.get_selected())

        self.entry.set_text(self.data_interface.curpath)
        
        self.draw_graphics()

    def reduce_slit_cb(self):
        pass

    ###
    # END CALLBACKS
    ###

    def draw_graphics(self):
        if self.data_interface.browse_cb_is_dir:
            self.update_tree_and_text()
        else:
            if Path(self.data_interface.curpath).name.startswith('Slits'):
                self.draw_all_slits()
            else:
                hdul = self.data_interface.get_raw_hdul()
                self.show_deimos_raw(hdul)

    def update_tree_and_text(self):
        self.entry.set_text(self.data_interface.curpath)

        self.treeview.set_tree(self.data_interface.tree_dict)

        if self.data_interface.resize_table_columns:
            self.treeview.set_optimal_column_widths()
            self.logger.debug("Resized columns for {0} row(s)".format(len(self.data_interface.tree_dict)))

    def interface_is_local(self):
        return self.data_interface.type == 'local'

    def close(self):
        if self.fitsimage is None:
            self.fv.stop_global_plugin(str(self))
        else:
            self.fv.stop_local_plugin(self.chname, str(self))
        return True

    def _add_info(self, channel, filelist):
        for path in filelist:
            name = iohelper.name_image_from_path(path)
            info = Bunch.Bunch(name=name, path=path)
            self.fv.gui_call(channel.add_image_info, info)

    def start(self):
        self.win = None
        self.logger.info("Starting QL plugin with browsing to home path: ")
        self.logger.info(self.data_interface.settings.get('home_path'))
        self.data_interface.browse(self.data_interface.settings.get('home_path'))
        self.update_tree_and_text()

        self.fv.get_channel_on_demand('pypeit_ql').fitsimage.get_canvas().add_callback('cursor-down', self.tess) 
        self.fv.change_channel('pypeit_ql')


    def stop(self):
        pass

    def redo(self, *args):
        return True

    def __str__(self):
        return 'qldisplay'

    def clear_canvas(self, cname):
        """
        Clear the ginga canvas

        Args:
            cname (str):  Channel name

        """
        # viewer = connect_to_ginga()
        ch = self.fv.channel(cname)
        canvas = self.fv.canvas(ch._chname)
        canvas.clear()


    def clear_all(self, allow_new=False):
        """
        Clear all of the ginga canvasses.

        Args:
            allow_new (:obj:`bool`, optional):
                Allow a subprocess to be called to execute a new ginga viewer if one
                is not already running.  See :func:`connect_to_ginga`.
        """
        # viewer = self.fv
        shell = self.fv
        chnames = shell.get_channel_names()
        for ch in chnames:
            shell.delete_channel(ch)

    def reduce_slit(self):
        
        slitspatid = self.slit_list.get_text()[1:]
        self.logger.info(f"Reducing slit {slitspatid}")
        det = ""
        if self.msc_idx == 0:
            det = "1,5"
        elif self.msc_idx == 1:
            det = "2,6"
        elif self.msc_idx == 2: 
            det = "3,7"
        elif self.msc_idx == 3:
            det = "4,8"

        reducer = LocalInterface(logger=self.logger, raw_path=self.data_dir, raw_files=self.raw_file, sci_files=self.raw_file, setup_calib_dir=self.cal_dir, redux_path="/Users/mbrodheim/drp/deimos_test/ql_test", det=det, slitspatnum=f"MSC0{self.msc_idx + 1}:{slitspatid}")
        reducer.reduce_slit()
        # self.logger.info(f"Command is: pypeit_ql keck_deimos --raw_path {self.data_dir} --raw_files {self.raw_file} --sci_files {self.raw_file} --setup_calib_dir {self.cal_dir} --redux_path /Users/mbrodheim/drp/deimos_test/ql_test --det {det} --slitspatnum MSC0{self.msc_idx + 1}:{slitspatid}")
        # Popen(["pypeit_ql", "keck_deimos", "--raw_path", self.data_dir, "--raw_files", self.raw_file, "--sci_files", self.raw_file, "--setup_calib_dir", self.cal_dir, "--redux_path", "/Users/mbrodheim/drp/deimos_test/ql_test", "--det", det, "--slitspatnum", f"MSC0{self.msc_idx + 1}:{slitspatid}"])


    def tess(self, canvas, pnt, x, y):
        self.logger.info("x: {0}, y: {1}".format(x, y))
        
        for msc_idx in self.data_interface.slit_data.keys():
            slits = self.data_interface.slit_data[msc_idx]
            if slits is None:
                continue
            offset = int(msc_idx) * slits.nspat
            left_bound_at_y = slits.left_init[np.round(y).astype(int)] + offset
            right_bound_at_y = slits.right_init[np.round(y).astype(int)] + offset

            for i in range(slits.nslits):
                if left_bound_at_y[i] < x < right_bound_at_y[i]:
                    slit_id = slits.spat_id[i]
                    self.logger.info("Found slit {0}".format(slit_id))
                    self.slit_list.show_text(f"S{slit_id}")
                    break


    def draw_all_slits(self):
        for msc_idx in self.data_interface.slit_data.keys():
                    
            if self.data_interface.slit_data[msc_idx] is None:
                continue
            slittrace = self.data_interface.slit_data[msc_idx]
            spatial_ids = slittrace.spat_id
            left_init = slittrace.left_init.T
            y_values_left = np.arange(slittrace.nspec)[::10]
            print(len(y_values_left))
            right_init = slittrace.right_init.T
            y_values_right = np.arange(slittrace.nspec)[::-10]
            print(len(y_values_right))

            for idx, spat_id in enumerate(spatial_ids):
                self.slit_list.append_text("S" + str(spat_id))

                x_vals = np.concatenate((left_init[idx][::10], right_init[idx][::-10]), axis=0) + (int(msc_idx) * slittrace.nspat)
                y_vals = np.concatenate((y_values_left, y_values_right), axis=0)
                slit_boundard_coords = (x_vals, y_vals)
                print(f"Slit {spat_id}")
                print(f"X len: {len(slit_boundard_coords[0])}")
                print(f"Y len: {len(slit_boundard_coords[1])}")
                print(f"nspec: {slittrace.nspec}")
                poly = Polygon(list(zip(slit_boundard_coords[0], slit_boundard_coords[1])), color='green', linewidth=1, fill=True, fillalpha=.1)
                
                canvas = DrawingCanvas()
                canvas.add(poly)
                self.fv.get_channel(None).fitsimage.set_canvas(canvas)

    def show_deimos_raw(self, hdul):
        # hdul = pyfits.open(filename)
    
        hdr0 = hdul[0].header
        
        ext = np.arange(1,9)
        
        binning = hdr0['BINNING'].split(',')
        
        precol =   int(hdr0['PRECOL'])   // int(binning[0])
        postpix =  int(hdr0['POSTPIX'])  // int(binning[0])
        preline =  int(hdr0['PRELINE'])  // int(binning[1])
        postline = int(hdr0['POSTLINE']) // int(binning[1])
        
        alldata = []
        for i in ext:
            data = hdul[i].data
            
            height, width = hdul[i].shape
            

            # get bias
            x1 = 0
            x2 = height
            y1 = width - postpix
            y2 = width
            
            bias = np.median(data[x1:x2, y1:y2], axis = 1)
            bias = np.array(bias, dtype = np.int64)            
            # bias subtraction
            data = data - bias[:, None]       
            
            # remove overscan
            data = data[:, precol: width - postpix]
            
            # append all 8 biased arrays into a list
            alldata.append(data)
        

        
        # creating CCD mosaic rows
        r0 = np.concatenate(alldata[:4], axis=1)
        r0 = np.flipud(r0)
        r1 = []
        for arr in alldata[4:]:
            arr = np.fliplr(arr)
            r1.append(arr)
        r1 = np.concatenate(r1, axis=1)

        fulldata = np.concatenate((r1, r0), axis=0)
        # fulldata = np.rot90(fulldata)
        
        img = AstroImage(logger=self.logger)
        img.load_data(fulldata)
        img['name'] = 'deimos_raw'
        # img.update_metadata({'name' : 'deimos_raw'})
        # self.fv.set_image(img)

        channel = self.fv.get_channel_on_demand("pypeit_ql")
        self.fv.change_channel(channel.name)
        self.logger.info("Image name: {0}".format(img.name))
        channel.add_image(img)

        # Display the image
        # channel = self.fv.gui_call(self.fv.get_channel_on_demand, "pypeit_ql")

        # # Note: this little hack needed to let window resize in time for
        # # file to auto-size properly
        # self.fv.gui_do(self.fv.change_channel, channel.name)
        # self.fv.gui_do(self.fv.add_image, "deimos_raw", img, chname=channel.name)
        return 0
    
class DRPInterface():
    """Abstract class for DRP interface.

    The interface connects to the GUI in a few places:
    - The treeview (file explorer)
    - The text entry (shows the currently selected file)
    - The reduce button, which initiates the reduction process
    """

    def __init__(self, settings, icons, logger):
        self.logger = logger
        self.settings = settings

        # Icons. Not sure if this should actually be here!
        self.folderpb = icons['folderpb']
        self.filepb = icons['filepb']
        self.fitspb = icons['fitspb']

        # Raw HDUL for the base image display
        self.raw_hdul = None

        # Current path
        self.curpath = os.path.join(self.settings.get('home_path'), '*')

        # Slit information
        self.slit_data = {
            '0' : None,
            '1' : None,
            '2' : None,
            '3' : None
        }

        # After selecting something, was it a file or a directory?
        self.browse_cb_is_dir = False # Used to determine if the text entry is a directory or file

    def set_slit_info(self, raw_path, raw_files, sci_files, setup_calib_dir, redux_path, det, slitspatnum) -> None:
        self.raw_path = raw_path
        self.raw_files = raw_files
        self.sci_files = sci_files
        self.setup_calib_dir = setup_calib_dir
        self.redux_path = redux_path
        self.det = det
        self.slitspatnum = slitspatnum

    def get_raw_hdul(self):
        """Serves as the output for loading a raw file.
        """
        return self.raw_hdul
    
    def get_selected_file(self):
        """Returns the path of the selected file, for displaying to the user
        """
        return self.selected_file
    
    def get_reduction_info(self):
        """Returns info needed to reduce a slit.
        """
        return {
            "raw_path" : self.raw_path,
            "raw_files" : self.raw_files,
            "sci_files" : self.sci_files,
            "setup_calib_dir" : self.setup_calib_dir,
            "redux_path" : self.redux_path,
            "det" : self.det,
            "slitspatnum" : self.slitspatnum
        }
    
    def reduce_slit(self):
        raise NotImplementedError
    
    def item_dblclicked_cb(self, widget, res_dict):
        raise NotImplementedError
    
    def item_selected_cb(widget, res_dict):
        raise NotImplementedError
    
    def browse_cb(self, widget):
        raise NotImplementedError
    
    def load_cb(self, text_entry_text, treeview_selected_dict):
        raise NotImplementedError

    def refresh():
        raise NotImplementedError

class LocalInterface(DRPInterface):
    """Provides an interface for local file exploration and reduction.

    Contains both interfaces for the ginga viewer to populate the file explorer,
    and commands to initiate data reduction.
    """

    def __init__(self,  settings, icons, logger):
        super().__init__(settings, icons, logger)


    def get_hdul_from_path(self, path=None):
        """Takes a path to a fits file and returns the hdul

        Parameters
        ----------
        path : path-like
            Path to the fits file. If none, uses self.raw_file
        """

        if path is None:
            path = self.raw_file
        
        hdul = pyfits.open(path)
        return hdul

    def load_and_display_files(self, paths):
        if self.fitsimage is not None:
            self.fv.gui_do(self.fitsimage.make_callback, 'drag-drop', paths)
        else:
            channel = self.fv.get_channel_info()
            if channel is None:
                chname = None
            else:
                chname = channel.name
            self.fv.gui_do(self.fv.open_uris, paths, chname=channel.name)

    def load_cb(self, text_entry_text, treeview_selected_dict):
        # Load from text box

        path = text_entry_text

        if os.path.isdir(path):
            self.browse(path)
            self.browse_cb_is_dir = True
        elif os.path.isfile(path):
            self.open_file(path)
            self.browse_cb_is_dir = False
        else:
            # Load from tree view
            #curdir, curglob = os.path.split(self.curpath)
            select_dict = treeview_selected_dict
            path = [info.path for key, info in select_dict.items()][0]
            self.logger.debug('Loading {0}'.format(path))

            # Open directory
            if os.path.isdir(path):
                # path = os.path.join(paths[0], '*')
                # self.entry.set_text(path)
                self.curpath = path
                self.browse(path)
                return

            # Exclude directories
            paths = [path for path in paths if os.path.isfile(path)][0]

            # Load files
            self.load_and_display_files(paths)

    def makelisting(self):

        tree_dict = {}
        for bnch in self.jumpinfo:
            icon = self.file_icon(bnch)
            bnch.setvals(icon=icon)
            entry_key = bnch.name

            if entry_key is None:
                raise Exception("No key for tuple")

            tree_dict[entry_key] = bnch

            import pprint
            print(entry_key)
            pprint.pprint(bnch)

        self.tree_dict = tree_dict

        # Do we need to resize column widths?
        n_rows = len(tree_dict)
        if n_rows < self.settings.get('max_rows_for_col_resize', 5000):
            self.resize_table_columns = True

    def get_path_from_item(self, res_dict):
        paths = [info.path for key, info in res_dict.items()]
        path = paths[0]
        return path

    def item_dblclicked_cb(self, widget, res_dict):
        self.logger.info("Item double clicked")
        path = self.get_path_from_item(res_dict)
        self.open_file(path)

    # def item_drag_cb(self, widget, drag_pkg, res_dict):
    #     urls = [Path(info.path).as_uri() for info in res_dict.values()]
    #     self.logger.info("urls: %s" % (urls))
    #     # destination can collect selection in two ways
    #     drag_pkg.set_urls(urls)
    #     drag_pkg.set_text('\n'.join(urls))

    def browse_cb(self, widget):
        path = str(widget.get_text()).strip()

        if os.path.isdir(path):
            self.browse(path)
            self.browse_cb_is_dir = True
        elif os.path.isfile(path):
            self.open_file(path)
            self.browse_cb_is_dir = False

        # # Load file(s) -- image*.fits, image*.fits[ext]
        # retcode = self.open_files(path)

        # # Open directory
        # if not retcode:
        #     self.browse(path)

    def item_selected_cb(self, widget, res_dict):
        paths = [info.path for info in res_dict.values()]
        n_paths = len(paths)
        if n_paths <= 0:
            return
        elif n_paths == 1:
            txt = paths[0]
        else:
            txt = ' '.join(['"{}"'.format(s) for s in paths])
        self.selected_file = txt
        self.curpath = paths[0]

    def file_icon(self, bnch):
        if bnch.type == 'dir':
            pb = self.folderpb
        elif bnch.type == 'fits':
            pb = self.fitspb
        else:
            pb = self.filepb
        return pb
    

    def open_file(self, path):
        self.logger.info("open_file path: %s" % (path))
        self.curpath = path

        if path == '..':
            curdir, curglob = os.path.split(self.curpath)
            path = os.path.join(curdir, path, curglob)

        if os.path.isdir(path):
            path = os.path.join(path, '*')
            p = Path(path)
            # if p.name.startswith('keck_deimos'):
            # #     pass
            #     self.logger.info("Found a setup directory")
            #     p = p / "Calibrations"
            #     # edge_files = list(p.glob("Edges*"))
            #     slits_files = list(p.glob("Slits*"))
                # self.logger.info(f"Found {len(edge_files)} edges files in {p.absolute()}")
                # self.logger.info(f"Found {len(slits_files)} slits files in {p.absolute()}")
                # if len(edge_files) != len(slits_files):
                #     self.logger.info("Number of edges files does not match number of slits files")
                #     return
                
                # self.mosaics = dict()
                # for i in range(len(edge_files)): # load each slits and edges file
                #     msc_idx = edge_files[i].name.split('MSC')[1]
                #     self.msc_idx = int(msc_idx.split('.')[0]) - 1
                #     self.logger.info(f"Opening mosaic number {self.msc_idx}...")

                #     self.logger.info(f"Opening edges {edge_files[i]}")
                #     # self.edges = edgetrace.EdgeTraceSet.from_file(edge_files[i])
                #     self.edges = None
                #     self.logger.info(f"Opening slits {slits_files[i]}")
                #     self.slits = slittrace.SlitTraceSet.from_file(slits_files[i])
                #     # pdb.set_trace()
                #     for slit in self.slits.spat_id:
                #         self.slit_list.append_text(f"S{slit}")
                    
                #     self.mosaics[self.msc_idx] = {"slits" : self.slits, "edges" : self.edges, "idx_name" : self.msc_idx + 1}
                # self.draw_slits()
            # else:
            #     self.browse(path)
            self.logger.info("Browsing to directory: %s" % (path))
            self.browse_cb_is_dir = True
            self.browse(path)

        elif os.path.exists(path):
            p = Path(path)
            self.browse_cb_is_dir = False
            
            if p.name.startswith('Slits'):
                # open the slits file
                msc_idx = int(p.name.split('MSC')[1].split('.')[0]) - 1
                self.slit_data[f'{msc_idx}'] = slittrace.SlitTraceSet.from_file(p)
                    
                # extract the left and right edges
                # create polygons and plot em
            elif p.suffix == '.fits':
                self.raw_file = p.name
                self.data_dir = p.parent

                # Set raw file here:
                self.raw_file = path
                self.logger.info("Getting fits image from path: {0}".format(path))
                self.raw_hdul = self.get_hdul_from_path(path)

        else:
            self.browse(path)

    def get_info(self, path):
        dirname, filename = os.path.split(path)
        name, ext = os.path.splitext(filename)
        ftype = 'file'
        if os.path.isdir(path):
            ftype = 'dir'
        elif os.path.islink(path):
            ftype = 'link'
        elif ext.lower() == '.fits':
            ftype = 'fits'
        
        na_dict = {attrname: 'N/A' for colname, attrname in self.settings.get('columns')}
        bnch = Bunch.Bunch(na_dict)
        try:
            filestat = os.stat(path)
            bnch.update(dict(path=path,
                            name=filename,
                            type=ftype
                            ))
        except OSError as e:
            # TODO: identify some kind of error with this path
            bnch.update(dict(path=path, name=filename, type=ftype,
                            st_mode=0, st_size=0,
                            st_mtime=0))

        return bnch

    def browse(self, path):
        self.logger.info("browse path: %s" % (path))
        self.curpath = path
        if os.path.isdir(path):
            dirname = path
            globname = None
        else:
            dirname, globname = os.path.split(path)
        dirname = os.path.abspath(dirname)

        # check validity of leading path name
        if not os.path.isdir(dirname):
            self.fv.show_error("Not a valid path: %s" % (dirname))
            return

        if not globname:
            globname = '*'
        path = os.path.join(dirname, globname)

        # Make a directory listing
        self.logger.debug("globbing path: %s" % (path))
        filelist = list(glob.glob(path))
        filelist.sort(key=lambda s: s.lower())
        filelist.insert(0, os.path.join(dirname, '..'))

        self.jumpinfo = list(map(self.get_info, filelist))
        import pprint
        pprint.pprint(self.jumpinfo)
        self.curpath = path

        if self.settings.get('scan_fits_headers', False):
            print("Scanning fits headers!")
            num_files = len(self.jumpinfo)
            if num_files <= self.settings.get('scan_limit'):
                self.scan_fits()
            else:
                self.logger.warning(
                    "Number of files (%d) is greater than scan limit (%d)"
                    "--skipping header scan" % (num_files, self.scan_limit))

        self.makelisting()

    def scan_fits(self):
        # Scan each FITS file and add header items
        self.logger.info("scanning files for header keywords...")
        start_time = time.time()
        for bnch in self.jumpinfo:
            if (not bnch.type == 'fits'):
                continue
            try:
                with pyfits.open(bnch.path, 'readonly') as in_f:
                    kwds = {attrname: in_f[0].header.get(kwd, 'N/A')
                            for attrname, kwd in self.settings.get('keywords')}
                bnch.update(kwds)
            except Exception as e:
                self.logger.warning(
                    "Error reading FITS keywords from "
                    "'%s': %s" % (bnch.path, str(e)))
                continue
        elapsed = time.time() - start_time
        self.logger.info("done scanning--scan time: %.2f sec" % (elapsed))

    def refresh(self):
        self.browse(self.curpath)

    def scan_headers(self):
        self.browse(self.curpath)

    def make_thumbs(self):
        path = self.curpath
        self.logger.info("Generating thumbnails for '%s'..." % (
            path))
        filelist = glob.glob(path)
        filelist.sort(key=lambda s: s.lower())

        if self.fitsimage is not None:
            # we were invoked as a local plugin
            channel = self.channel
        else:
            chviewer = self.fv.getfocus_viewer()
            # find out our channel
            chname = self.fv.get_channel_name(chviewer)
            channel = self.fv.get_channel(chname)

        self.fv.nongui_do(self._add_info, channel, filelist)

    def reduce_slit(self):
        self.logger.info(f"Command is: pypeit_ql keck_deimos --raw_path {self.data_dir} --raw_files {self.raw_file} --sci_files {self.raw_file} --setup_calib_dir {self.cal_dir} --redux_path /Users/mbrodheim/drp/deimos_test/ql_test --det {self.det} --slitspatnum MSC0{self.msc_idx + 1}:{self.slitspatid}")
        Popen(["pypeit_ql", "keck_deimos", "--raw_path", self.data_dir, "--raw_files", self.raw_file, "--sci_files", self.raw_file, "--setup_calib_dir", self.cal_dir, "--redux_path", "/Users/mbrodheim/drp/deimos_test/ql_test", "--det", self.det, "--slitspatnum", f"MSC0{self.msc_idx + 1}:{self.slitspatid}"])
        

class RemoteInterface(DRPInterface):
    """Provides an interface for remote file exploration and reduction.
    """

    def __init__(self, ) -> None:
        pass
        
    

# Append module docstring with config doc for auto insert by Sphinx.
from ginga.util.toolbox import generate_cfg_example  # noqa
if __doc__ is not None:
    __doc__ += generate_cfg_example('plugin_QLDisplay', package='ginga')

# END