"""
This module contains non I{Traits} native user interfaces.
GUIs are created using wxPython or PyQt directly.

@author: Karsten Mueller
"""

import os

#def QtOpenFile(self, dir=os.getcwd(), multi=True):
#    import sys
#    from PyQt4 import QtGui
#    
#    extfilter = '''Raw C-band Files (*.gYY);;Processed C-band Files (*.bYY);;Ramac Files (*.rd3);;All Files (*.*)'''
#    
#    app = QtGui.QApplication(sys.argv)
#    dlg = QtGui.QFileDialog
#    if multi:
#        pathlist = dlg.getOpenFileNames(None, "Open Files", dir, extfilter)
#        pathlist.sort()
#        filelist = []
#        for path in pathlist:
#            filelist.append(os.path.basename(str(path)))
#        dir = os.path.dirname(str(pathlist[0]))
#        return filelist
#    else:
#        pathlist = dlg.getOpenFileName(None, "Open File", dir, extfilter)
#        filelist = os.path.basename(str(pathlist))
#        dir = os.path.dirname(str(pathlist))
#        return filelist
#    os.chdir(dir)        


def wxOpenFile(dir=os.getcwd(), multi=True, filter=0):
    import wx

    if filter == 0:
        extfilter = "HDF 5 Files (*.h5)|*.h5|All Files (*.*)|*.*|Raw C-band Files (*.gYY)|*.gYY|Raw C-band Files (*.gYX)|*.gYX|Raw C-band Files (*.gXX)|*.gXX|Raw C-band Files (*.gXY)|*.gXY|Processed C-band Files (*.bYY)|*.bYY|Ramac Files (*.rd3)|*.rd3"
    elif filter == 1:
        extfilter = "Pick files (*.pck)|*.pck|All Files (*.*)|*.*"
    else:
        extfilter = "All Files (*.*)|*.*"

    #    app = wx.PySimpleApp()

    if multi:
        dlg = wx.FileDialog(None, message='Select files ...', defaultDir=dir, wildcard=extfilter, style=wx.FD_MULTIPLE)
        if dlg.ShowModal() == wx.ID_OK:
            pathlist = dlg.GetPaths()
            filelist = []
            for path in pathlist:
                filelist.append(os.path.basename(str(path)))
            filelist.sort()
            os.chdir(os.path.dirname(str(pathlist[0])))
            path = os.path.dirname(str(pathlist[0]))
            return filelist, path
        else:
            return False
        dlg.Destroy()

    else:
        dlg = wx.FileDialog(None, message='Select a file ...', defaultDir=dir, wildcard=extfilter, style=wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            path = str(dlg.GetPath())
            datafile = os.path.basename(str(path))
            dir = os.path.dirname(str(path))
            return datafile, dir
        else:
            return False
        dlg.Destroy()


def wxSaveFile(dir=os.getcwd()):
    import wx

    extfilter = "HDF 5 Files (*.h5)|*.h5|All Files (*.*)|*.*"
    dlg = wx.FileDialog(None, message='Save as...', defaultDir=dir, wildcard=extfilter, style=wx.FD_SAVE)
    if dlg.ShowModal() == wx.ID_OK:
        path = str(dlg.GetPath())
        datafile = os.path.basename(str(path))
        dir = os.path.dirname(str(path))
        return datafile, dir
    else:
        return False
    dlg.Destroy()