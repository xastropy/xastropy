"""A collection of tools, tips, and tricks.

2009-07-20 22:36 IJC: Created

2010-10-28 11:53 IJMC: Updated documentation for Sphinx.

2011-06-15 09:34 IJMC: More functions have been added; cleaned documentation.
"""

import pdb
import numpy as np

def getfigs():
    """Return a list of all open matplotlib figures.

    No inputs or options."""
    from matplotlib._pylab_helpers import Gcf

    figs = [manager.canvas.figure for manager in Gcf.get_all_fig_managers()]
    figlist = [fig.number for fig in figs]

    return figlist

def nextfig():
    """Return one greater than the largest-numbered figure currently
    open. If no figures are open, return unity.

    No inputs or options."""
    # 2010-03-01 14:28 IJC: Created
    figlist = getfigs()
    if len(figlist)==0:
        return 1
    else:
        return max(figlist)+1
        
    return figlist

def printfigs(filename, figs=None, format=None, pdfmode='texexec', verbose=False, closefigs=False):
    """Print desired figures using designated 'format'. Concatenate PDFs.

    :Inputs:
        filename -- string.  prepended to all open figures

        figs -- int or list. 

          figures to access, then apply savefig to.  If None, print
                   all open figures; if -1, print current figure.

        format -- string or list of strings. 

          if 'pdf', all images are concatenated into one file (use
                   "pdfs" for individual pdf figure files)

        pdfmode -- string; 

           method of concatenating PDFs.  Either 'texexec' or 'gs' for
                   GhostScript

        closefigs -- bool

           If True, close each figure after printing it to disk.


    :NOTES:
      If no explicit path is passed and a subdirectory 'figures'
      exists in the current directory, the figures will be printed in
      'figures' instead.

    :EXAMPLE:
      ::

        from pylab import *
        figure(1); plot(arange(10), randn(10), 'ob')
        figure(2); plot(arange(15), randn(15), '-xr')
        printfigs('testing')
        !open testing.pdf

    """
    # 2009-07-20 23:10 IJC: Created; inspired by FGD.
    # 2009-09-08 13:54 IJC: Made it work with single-figure, non-list input.
    # 2010-02-02 11:50 IJC: Now it kills the 'logfile' detritus.
    # 2010-10-27 17:05 IJC: New texexec syntax is "result=...", not "result ..."
    # 2011-03-01 18:14 IJC: Added capability for multiple formats (in
    #                       a list).  Also, figure numbers are not
    #                       catted to the filename when saving a
    #                       single figure.
    # 2011-08-29 10:23 IJMC: Now don't try to concatenate single PDF figures.
    # 2012-11-01 11:41 IJMC: Slightly changed if-block for 'figs'.
    # 2014-05-03 15:04 IJMC: Added 'closefigs' flag.

    from pylab import savefig, figure, gcf, close
    from matplotlib._pylab_helpers import Gcf
    import os
    import pdb

    figlist = getfigs()
    if verbose: print "Available figure numbers>>" ,figlist
    if figs is None:
        figs = figlist
    elif figs is -1:
        figs = [gcf().number]
    else:
        if hasattr(figs, '__iter__'):
            figs = list(figs)
        else:
            figs = [figs]
    figlist = [val for val in figs if val in figlist]
    nfig = len(figlist)
    print "Figures to print>>",figlist

    if format==None:
        format = filename[-3::]
        filename = filename[0:len(filename)-4]

    if hasattr(format, 'capitalize'):
        format = [format]
        nformat = 1
    elif hasattr(format, '__iter__'):
        nformat = len(format)
    else:
        format = [str(format)]
        nformat = 1
        
    if len(figlist)==0:
        print "No open figures found; exiting."
        return
    
    for thisformat in format:
        fnamelist = []
        for ii in range(nfig):
            if nfig>1:
                fname = filename + str(figlist[ii]) 
            else:
                fname = filename 
            if thisformat=='pdf' and nfig>1:
                fname = fname + '_temp'
            if thisformat=='pdfs':
                fname = fname + '.pdf'
            else:
                fname = fname + '.' + thisformat
            figure(figlist[ii])
            savefig(fname )
            fnamelist.append(fname)
    
            if closefigs and thisformat==format[-1]:  # last time at this figure
                close(figlist[ii])

        if thisformat=='pdf':
            if nfig==1:
                savefig(fnamelist[0])
            else:  # we have to concatenate multiple PDF figures:
                bigfilename = filename + '.' + thisformat
                if os.path.isfile(bigfilename):
                    os.remove(bigfilename)
                if pdfmode=='gs':
                    execstr = 'gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile='
                    rmstr = ''
                elif pdfmode=='texexec':
                    execstr = 'texexec --pdfcopy --result='
                    rmstr = 'rm %s' % bigfilename.replace('pdf','log')
                execstr += bigfilename
                for fn in fnamelist:
                    execstr += ' ' + fn
                #pdb.set_trace()
                if verbose: print "GS exec call>>", execstr
                os.system(execstr)
                #pdb.set_trace()
                if len(rmstr)>0:
                    os.system(rmstr)
                for fn in fnamelist:
                    try:
                        os.remove(fn)
                    except:
                        pass
            

    return


def plotstyle(i, c=['b', 'g', 'r', 'c', 'm', 'y', 'k'], \
                  s=['.', 'x', 's', '^', '*', 'o', '+', 'v', 'p', 'D'], \
                  l=['-', '--', '-.', ':']):
    """Return plot properties to help distinguish many types of plot symbols.

    :INPUT:
       i -- int.  

    :OPTIONAL INPUT:
       c -- color, or list of colors accepted by pylab.plot

       s -- symbol, or list of symbols accepted by pylab.plot

       l -- linestyle, or list of linestyles accepted by pylab.plot

    :OUTPUT:
       tuple of (color, symbol, linestyle)

    :REQUIREMENTS:  :doc:`numpy`
       """
    # 2009-09-10 16:42 IJC: Created
    from numpy import tile, array

    if not c.__class__==list:
        c = list(c)
    if not s.__class__==list:
        s = list(s)
    if not l.__class__==list:
        l = list(l)
    nc = len(c)
    ns = len(s)
    nl = len(l)

    if not hasattr(i,'__iter__'):
        i = array([i])
    i = abs(array(i))

    nrepc = (max(i)/nc+1.).astype(int)
    nreps = (max(i)/ns+1.).astype(int)
    nrepl = (max(i)/nl+1.).astype(int)
    c = tile(c, nrepc)
    s = tile(s, nreps)
    l = tile(l, nrepl)
    if len(i)==1:
        ret = c[i][0], s[i][0], l[i][0]
    else:
        ret = list(c[i]),list(s[i]),list(l[i])
    return ret


def flatten(L, maxdepth=100):
    """Flatten a list. 

    Stolen from http://mail.python.org/pipermail/tutor/2001-January/002914.html"""
    # 2009-09-10 16:54 IJC: Input.
    if type(L) != type([]): return [L]
    if L == []: 
        return L
    else:
        maxdepth -= 1
        return flatten(L[0]) + flatten(L[1:], maxdepth=maxdepth)

def replaceall(seq, obj, rep):
    """Replace all instances of 'obj' with 'rep' in list 'seq'

    :INPUT:
       seq -- (list) list within which to find-and-replace elements

       obj -- target object to replace

       rep -- replacement object

    :EXAMPLE:
      ::

          import tools
          b = [2, ['spam', ['eggs', 5, dict(spam=3)]]]
          tools.replaceall(b, 'spam', 'bacon')
          print b

    :NOTES:
       -- Will fail if 'obj' is itself a list.

       -- Edits list in-place, so make a copy first if you want to
          retain the old version of your list.

       -- Has not been tested for extremely deep lists 

    :SEE ALSO:
       :func:`popall`   
       """
    #2009-09-11 10:22 IJC: Created
    n = len(seq)
    
    for ii in range(n):
        if seq[ii].__class__==list:
            replaceall(seq[ii], obj, rep)
        else:
            if seq[ii]==obj:
                seq[ii]=rep

    return

def popall(seq, obj):
    """Remove all instances of 'obj' from list 'seq'

    :INPUT:
       seq -- (list) list from which to pop elements

       obj -- target object to remove

    :EXAMPLE:
      ::

           import tools
           b = [3, 'spam', range(5)]
           tools.popall(b, 4)
           print b

    :NOTES:
       -- Will fail if 'obj' is itself a list.

       -- Edits list in-place, so make a copy first if you want to
          retain the old version of your list.

       -- Has not been tested for extremely deep lists 

    :SEE ALSO:
       :func:`replaceall`
       """
    #2009-09-11 10:22 IJC: Created
    n = len(seq)
    
    for ii in range(n):
        print ii,seq[ii]
        if seq[ii].__class__==list:
            popall(seq[ii], obj)

    doneYet = False
    while not doneYet:
        try:
            seq.remove(obj)
        except:
            doneYet = True

    return

def drawRectangle(x,y,width,height,**kw):
    """Draw a rectangle patch on the current, or specified, axes.

    :INPUT:
       x, y -- lower-left corner of rectangle

       width, height -- dimensions of rectangle

    :OPTIONAL INPUT:
       ax -- Axis to draw upon. if None, defaults to current axes.  

       dodraw -- if True, call 'draw()' function to immediately re-draw axes.

       **kw -- options passable to :func:`matplotlib.patches.Rectangle`

     :NOTE:  Axes will NOT auto-rescale after this is called.
      """
    # 2009-09-17 01:33 IJC: Created
    # 2014-03-01 13:51 IJMC: Added 'dodraw' option.
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    
    if kw.has_key('ax'):
        ax = kw.pop('ax')
    else:
        ax = plt.gca()
    p = mpatches.Rectangle((x,y), width, height, **kw)
    ax.add_patch(p)
    if kw.has_key('dodraw') and kw['dodraw']: plt.draw()

    return ax, p

def drawPolygon(xy,**kw):
    """Draw a rectangle patch on the current, or specified, axes.

    :INPUT:
       xy -- numpy array of coordinates, with shape Nx2.

    :OPTIONAL INPUT:
       ax -- Axis to draw upon. if None, defaults to current axes.  

       dodraw -- if True, call 'draw()' function to immediately re-draw axes.

       **kw -- options passable to :func:`matplotlib.patches.Polygon`

    :SEE ALSO:
       :func:`drawRectangle`
       
    :NOTE:  Axes will NOT auto-rescale after this is called.
      """
    # 2010-12-02 19:58 IJC: Created from drawRectangle
    # 2014-03-01 13:51 IJMC: Added 'dodraw' option.
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    
    if kw.has_key('ax'):
        ax = kw.pop('ax')
    else:
        ax = plt.gca()

    p = mpatches.Polygon(xy, **kw)
    ax.add_patch(p)
    if kw.has_key('dodraw') and kw['dodraw']: plt.draw()

    return ax, p

def drawCircle(x,y,radius,**kw):
    """Draw a circular patch on the current, or specified, axes.

    :INPUT:
       x, y -- center of circle

       radius -- radius of circle

    :OPTIONAL INPUT:
       ax -- Axis to draw upon. if None, defaults to current axes.  

       dodraw -- if True, call 'draw()' function to immediately re-draw axes.

       **kw -- options passable to :func:`matplotlib.patches.Circle`

    :NOTE:  Axes will NOT auto-rescale after this is called.
      """
    # 2011-01-28 16:03 IJC: Created
    # 2014-03-01 13:51 IJMC: Added 'dodraw' option.
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    
    if kw.has_key('ax'):
        ax = kw.pop('ax')
    else:
        ax = plt.gca()
    p = mpatches.Circle((x,y), radius, **kw)
    ax.add_patch(p)
    if kw.has_key('dodraw') and kw['dodraw']: plt.draw()

    return ax, p


def drawEllipse(x,y,width, height,**kw):
    """Draw an elliptical patch on the current, or specified, axes.

    :INPUT:
       x, y -- center of ellipse

       width -- width of ellipse

       height -- width of ellipse

    :OPTIONAL INPUT:
       ax -- Axis to draw upon. if None, defaults to current axes.  

       dodraw -- if True, call 'draw()' function to immediately re-draw axes.

       **kw -- options passable to :func:`matplotlib.patches.Ellipse`
               (angle, linewidth, fill, ...)

    :NOTE:  Axes will NOT auto-rescale after this is called.

    :SEE_ALSO:
       :func:`drawCircle`, :func:`drawRectangle`

      """
    # 2011-10-20 11:32 IJMC: Created
    # 2014-03-01 13:51 IJMC: Added 'dodraw' option.
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    
    if kw.has_key('ax'):
        ax = kw.pop('ax')
    else:
        ax = plt.gca()
    p = mpatches.Ellipse((x,y), width, height, **kw)
    ax.add_patch(p)
    if kw.has_key('dodraw') and kw['dodraw']: plt.draw()
    return ax, p


def errxy(x,y,xbins, xmode='mean', ymode='mean', xerr='minmax', yerr='sdom', clean=None, binfactor=None, verbose=False,returnstats=False, timing=False):
    """Bin down datasets in X and Y for errorbar plotting

    :INPUTS:
       x -- (array) independent variable data

       y -- (array) dependent variable data

       xbins -- (array) edges of bins, in x-space.  Only x-data
                between two bin edges will be used.  Thus if M bin
                edges are entered, (M-1) datapoints will be returned.
                If xbins==None, then no binning is done.

    :OPTIONAL INPUT:
       xmode/ymode -- (str) method to aggregate x/y data into datapoints:
              'mean' -- use numpy.mean
              'median' -- use numpy.median
              'sum' -- use numpy.sum
              None -- don't compute; return the empty list []

       xerr/yerr -- (str) method to aggregate x/y data into errorbars
              'std' -- sample standard deviation (numpy.std)
              'sdom' -- standard deviation on the mean; i.e., std/sqrt(N)
              'minmax' -- use full range of data in the bin
              None -- don't compute; return the empty list []

       binfactor -- (int) If not None, average over this many
              consecutive values instead of binning explicitly by
              time-based bins.  Can also be a sequence, telling the
              number of values over which to average.  E.g.,
              binfactor=[10,10,20] will bin over the first 10 points,
              the second 10 points, and the next 20 points.

       clean -- (dict) keyword options to clean y-data ONLY, via
                analysis.removeoutliers, with an additional "nsigma"
                keyword.  See removeoutliers for more information.
                E.g.:  clean=dict(nsigma=5,remove='both',niter=1)

    :OUTPUTS:     a tuple of four arrays to be passed to matplotlib.pyplot.errorbar:
       xx -- locations of the aggregated x-datapoint in each bin

       yy -- locations of the aggregated y-datapoint in each bin

       xerr -- x-errorbars

       yerr -- y-errorbars

    :EXAMPLE:
      ::

          x = hstack((arange(10), arange(20)+40))
          y = randn(len(x))
          xbins = [-1,15,70]
          xx,yy,xerr,yerr = errxy(x,y,xbins)
          plot(x,y, '.b')
          errorbar(xx,yy,xerr=xerr,yerr=yerr, fmt='or')
      
    :NOTES:

       To just bin down uncleaned data (i.e., no 'error' terms
          returned), set clean, xerr, yerr to None.  However, when
          computing all values (xerr and yerr not None) it is faster
          to set clean to some rediculous value, i.e.,
          clean=dict(niter=0, nsigma=9e99).  This probably means more
          optimization could be done.

       Be sure you call the errorbar function using the keywords xerr
          and yerr, since otherwise the default order of inputs to the
          function is (x,y,yerr,xerr).

       Data 'x' are determined to be in a bin with sides (L, R) when
          satisfying the condition (x>L) and (x<=R)

    :SEE ALSO:  matplotlib.pyplot.errorbar, :func:`analysis.removeoutliers`

    :REQUIREMENTS:  :doc:`numpy`, :doc:`analysis`
    """
    # 2009-09-29 20:07 IJC: Created w/mean-median and std-sdom-minmax.
    # 2009-12-14 16:01 IJC: xbins can be 'None' for no binning.
    # 2009-12-15 10:09 IJC: Added "binfactor" option.
    # 2009-12-22 09:56 IJC: "binfactor" can now be a sequence.
    # 2009-12-29 01:16 IJC: Fixed a bug with binfactor sequences.
    # 2010-04-29 09:59 IJC: Added 'returnstats' feature
    # 2010-10-19 16:25 IJC: Added 'sum' option for x-data
    # 2011-03-22 12:57 IJC: Added 'none' option for data and errors
    # 2012-03-20 16:33 IJMC: Fixed bug; xmode=='none' now works.
    # 2012-03-27 14:00 IJMC: Now using np.digitize -- speed boost.
    #                        Rewrote code to optimize (somewhat),
    #                        cleaned up 'import' statements.
    # 2012-04-08 15:57 IJMC: New speed boost from adopting
    #                        numpy.histogram-like implementation:
    #                        numpy.searchsorted, etc.

    import numpy as np
    from analysis import removeoutliers


    if timing:
        import time
        tic = time.time()

    def sdom(data):
        """Return standard deviation of the mean."""
        return np.std(data)/np.sqrt(data.size)

    def getcenter(data, cmode):
        """Get data center based on mode.  Helper function."""
        if cmode is None:
            ret = 0
        elif cmode=='mean':
            ret = np.mean(data)
        elif cmode=='median':
            ret = np.median(data)
        elif cmode=='sum':
            ret = np.sum(data)
        return ret

    def geterr(data, emode, cmode):  
        """Get errorbar. Helper function."""
        if emode is None:
            ret = []
        elif emode=='std':
            ret = np.std(data)
        elif emode=='sdom':
            ret = sdom(data)
        elif emode=='minmax':
            if len(data)==0:
                ret = [np.nan, np.nan]
            else:
                center = getcenter(data,cmode)
                ret = [center-min(data), max(data)-center]
        return ret

    def cleandata(data, clean, returnstats=False):
        """Clean data using removeoutliers. Helper function."""
        init_count = np.array(data).size

        if clean==None: # Don't clean at all!
            #clean = dict(nsigma=1000, niter=0)
            if returnstats:
                ret = data, (init_count, init_count)
            else:
                ret = data

        else:  # Clean the data somehow ('clean' must be a dict)
            if not clean.has_key('nsigma'):
                clean.update(dict(nsigma=99999))
            data = removeoutliers(data, **clean)
            if returnstats:
                ret = data, (init_count, np.array(data).size)
            else:
                ret = data

        return ret

    if timing:
        print "%1.3f sec since starting function; helpers defined" % (time.time() - tic)

    ####### Begin main function ##########
    sorted_index = np.argsort(x)
    x = np.array(x, copy=False)[sorted_index]
    y = np.array(y, copy=False)[sorted_index]
    #x = np.array(x,copy=True).ravel()
    #y = np.array(y,copy=True).ravel()
    xbins = np.array(xbins,copy=True).ravel()
    if xbins[0]==None and binfactor==None:
        if returnstats ==False:
            ret = x, y, np.ones(x.shape)*np.nan, np.ones(y.shape)*np.nan
        else:
            ret = x, y, np.ones(x.shape)*np.nan, np.ones(y.shape)*np.nan, (x.size, x.size)
        return ret

    if binfactor==None:  # used passed-in 'xbins'
        xbins = np.sort(xbins)
    elif hasattr(binfactor,'__iter__'): # use variable-sized bins
        binfactor = np.array(binfactor).copy()
        sortedx = np.sort(x)
        betweens = np.hstack((x.min()-1, 0.5*(sortedx[1::]+sortedx[0:len(x)-1]), x.max()+1))
        xbins = []
        counter = 0
        for ii in range(len(binfactor)):
            thisbin = betweens[counter]
            xbins.append(thisbin)
            counter += binfactor[ii]
        xbins.append(x.max() + 1)
    else: # bin down by the same factor throughout
        binfactor = int(binfactor)
        sortedx = np.sort(x)
        betweens = np.hstack((x.min()-1, 0.5*(sortedx[1::]+sortedx[0:len(x)-1]), x.max()+1))
        xbins = betweens[::binfactor]

    if timing:
        print "%1.3f sec since starting function; bins defined" % (time.time() - tic)

    nbins = len(xbins)-1

    arraynan = np.array([np.nan])

    exx = []
    eyy = []
    xx = np.zeros(nbins)
    yy = np.zeros(nbins)
    yy2 = np.zeros(nbins)

    init_count, final_count = y.size, 0
    if timing:
        setuptime = 0
        xdatatime = 0
        ydatatime = 0
        statstime = 0

    #import pylab as py
    #xxx = np.sort(x)

    if timing: tic1 = time.time()
    #inds = np.digitize(x, xbins)
    inds2 = [[x.searchsorted(xbins[ii], side='left'), \
                  x.searchsorted(xbins[ii+1], side='left')] for ii in range(nbins)]
    if timing: setuptime += (time.time() - tic1)
    #pdb.set_trace()
    #bin_means = [data[digitized == i].mean() for i in range(1, len(bins))]



    dox = xmode is not None 
    doy = ymode is not None 
    doex = xerr is not None 
    doey = yerr is not None 

    if clean is None:
        if timing: tic3 = time.time()
        if dox: exec ('xfunc = np.%s' % xmode) in locals()
        if doy: exec ('yfunc = np.%s' % ymode) in locals()
        for ii in range(nbins):
            #index = inds==(ii+1)
            if dox:
                #xx[ii] = xfunc(x[index])
                xx[ii] = xfunc(x[inds2[ii][0]:inds2[ii][1]])
            if doy:
                #yy[ii] = yfunc(y[index])
                yy[ii] = yfunc(y[inds2[ii][0]:inds2[ii][1]])
            if doex:
                #exx.append(geterr(x[index], xerr, xmode))
                exx.append(geterr(x[inds2[ii][0]:inds2[ii][1]], xerr, xmode))
            if doey:
                #eyy.append(geterr(y[index], yerr, ymode))
                eyy.append(geterr(y[inds2[ii][0]:inds2[ii][1]], yerr, ymode))

        if timing: statstime += (time.time() - tic3)
        #pdb.set_trace()
    else:
        for ii in range(nbins):
            if timing: tic1 = time.time()
            #index = inds==(ii+1)
            if timing: setuptime += (time.time() - tic1)

            if timing: tic2 = time.time()
            xdata = x[inds2[ii][0]:inds2[ii][1]]
            if timing: xdatatime += (time.time() - tic2)

            if timing: tic25 = time.time()
            if ymode is None and yerr is None:  # We're free to ignore the y-data:
                ydata = arraynan
            else:  # We have to compute something with the y-data:
                if clean is not None:
                    ydata, retstats = cleandata(y[inds2[ii][0]:inds2[ii][1]], clean, returnstats=True)
                    if returnstats:
                        final_count += retstats[1]
                else:  # We don't have to clean the data
                    ydata = y[inds2[ii][0]:inds2[ii][1]]
                    if returnstats:
                        final_count += ydata.size
            if timing: ydatatime += (time.time() - tic25)

            if timing: tic3 = time.time()
            xx[ii] = getcenter(xdata,xmode)
            if timing: tic4 = time.time()
            yy[ii] = getcenter(ydata,ymode)
            if timing: tic5 = time.time()
            exx.append(geterr(  xdata,xerr,xmode))
            if timing: tic6 = time.time()
            eyy.append(geterr(  ydata,yerr,ymode))
            if timing: tic7 = time.time()
            if timing: statstime += (time.time() - tic3)
            #exx[ii] = geterr(  xdata,xerr,xmode)
            #eyy[ii] = geterr(  ydata,yerr,ymode)

    if timing:
        print "%1.3f sec for setting up bins & indices..." % setuptime
        print "%1.3f sec for getting x data clean and ready." % xdatatime
        print "%1.3f sec for getting y data clean and ready." % ydatatime
        #print "%1.3f sec for computing x-data statistics." % (tic4-tic3)
        #print "%1.3f sec for computing y-data statistics." % (tic5-tic4)
        #print "%1.3f sec for computing x-error statistics." % (tic6-tic5)
        #print "%1.3f sec for computing y-error statistics." % (tic7-tic6)


        print "%1.3f sec for computing statistics........." % statstime

    if timing:
        print "%1.3f sec since starting function; uncertainties defined" % (time.time() - tic)



    #xx = array(xx)
    #yy = array(yy)
    exx = np.array(exx).transpose()  # b/c 2D if minmax option used
    eyy = np.array(eyy).transpose()  # b/c 2D if minmax option used

    #pdb.set_trace()
 
    if returnstats:
        ret= xx,yy,exx,eyy,(init_count, final_count)
    else:
        ret = xx,yy,exx,eyy

    #print 'tools: returnstats, len(ret)>>', returnstats, len(ret)
    if timing:
        print "%1.3f sec since starting function; returning" % (time.time() - tic)

    return ret 


def ploth(*args, **kw):
    """Plot 1D data in a histogram-like format.  If x-coordinates are
    specified, they refer to the centers of the histogram bars.

    Uses same format as matplotlib.pyplot.plot.  For example:
      ::

          ploth(x, y)         # plot x and y using solid linestyle (default)
          ploth(x, y, 'bo')   # plot x and y using blue circle markers w/no line
          ploth(y)            # plot y using x as index array 0..N-1
          ploth(y, 'r*--')    # ditto, but with red star corners and dashed line

    :OPTIONS:
      rot90 : bool
        If True, data will be plotted histogram-style vertically,
        rather than the standard horizontal plotting.

    :REQUIREMENTS: :doc:`numpy`, :doc:`analysis`
       """
    # 2009-09-17 09:26 IJC: Created
    # 2012-09-27 19:19 IJMC: Added 'rot90' keyword

    from numpy import arange, concatenate, vstack
    from pylab import plot

    if len(args)==1:
        y=args[0]
        ny = len(y)
        x = arange(ny)
        plotstr = '-'
    elif len(args)==2 and args[1].__class__==str:
        y=args[0]
        ny = len(y)
        x = arange(ny)
        plotstr = args[1]
    elif len(args)==2 and args[1].__class__<>str:
        x = args[0]
        y=args[1]
        ny = len(y)
        plotstr = '-'
    elif len(args)>=3:
        x = args[0]
        y=args[1]
        ny = len(y)
        plotstr = args[1]

    if kw.has_key('rot90') and kw['rot90']:
        temp = x
        x = y
        y = temp
        ny = len(y)
        nx = len(x)
        rot90 = kw.pop('rot90')
    else:
        rot90 = False
        
    x1= 0.5*(x[1::]+x[0:ny-1])
    xx = concatenate(([x[0]], vstack((x1,x1)).transpose().ravel(), [x[-1]]))
    yy = vstack((y,y)).transpose().ravel()
    if rot90:
        phandle = plot(xx,yy,plotstr,**kw)
    else:
        phandle = plot(xx,yy,plotstr,**kw)

    return phandle

def flatten(x, maxdepth=100):
    """flatten(sequence) -> list

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    :OPTIONAL INPUTS:
      maxdepth -- scalar
         number of layers deep to dig.  Seting to zero causes no flattening to occur.

    :Examples:
      ::

        >>> [1, 2, [3,4], (5,6)]
        [1, 2, [3, 4], (5, 6)]
        >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
        [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""

    # 2009-09-26 14:05 IJC: Taken from 
    #     http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
    # 2011-06-24 15:40 IJMC: Added maxdepth keyword

    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, basestring) and maxdepth>0:
            maxdepth -= 1
            result.extend(flatten(el, maxdepth=maxdepth))
        else:
            result.append(el)
    return result

def fconvolve(a, v, oversamp=2):
    """Returns the discrete, linear convolution of 1-D sequences a and v,
    using Fast Fourier Transforms.  Restrictions are: a and v must
    both be real-valued, and len(a)>len(v).

    :REQUIREMENTS: :doc:`analysis`, :doc:`numpy`
    """
    # 2009-10-29 11:00 IJC: Created

    from analysis import pad
    from numpy.fft import fft, ifft, fftshift
    from numpy import real, array

    a = array(a,copy=True).ravel()
    v = array(v,copy=True).ravel()
    na = len(a)
    nv = len(v)
    nfft = oversamp*na

    a2 = pad(a, 1, nfft)[0,:]
    v2 = pad(v, 1, nfft)[0,:]
    
    fa2 = fft(a2)
    fv2 = fft(v2)
    ret = real(fftshift(ifft(fa2 * fv2)))
    
    return pad(ret, 1, na).ravel()
    
    
    

def cplim(a1,a2):
    """Copy axis limits from one axes to another.

    :INPUTS:
        a1, a2 -- either (1) handles to axes objects, or (2) figure
        numbers.  If figures have subplots, you can refer to a
        particular subplot using decimal notation. So, 1.3
        would refer to subplot 3 of figure 1.

    :REQUIREMENTS: :doc:`matplotlib` (when this is written...)
        """
    # 2009-12-08 16:30 IJC: Had the idea...

    print "To be written -- and what a great day it will be."

    return

    
    
                  

def legc(leg,col='color'):
    """Color legend text to match linecolor.

    :Inputs:
      'leg' is a legend object.

      'col' sets the field of the leg.get_lines() objects to use to find
          the color.

    You may need to refresh the figure to see the changes."""
    # 2009-12-14 09:50 IJC: Created

    texts = leg.get_texts()
    lines = leg.get_lines()
    for label,line in zip(texts,lines):
        label.set_color(line.get_color())
    return leg

def keylist(filelist, keys):
    """Create an object based on FITS header keys extracted from a filelist.

    :Inputs:
      filelist -- sequence of strings representing filenames (for PyFITS)

      keys -- sequence of strings representing header keys

    #Keys not found in a file will result in the string value

    :REQUIREMENTS: :doc:`pyfits`, :doc:`spitzer`
    """
    import pyfits
    from spitzer import baseObject
    # 2010-01-24 15:23 IJC: Created
    # 2010-01-26 10:27 IJC: Solved a pernicious little problem: always
    #                       call object creation w/parentheses!
    obj = baseObject()
    for k in keys:
        exec('obj.%s=[]'%k)
    for f in filelist:
        h = pyfits.getheader(f)
        for k in keys:
            exec("obj.%s.append(h['%s'])" %(k,k) )

    return obj

def plotc(x,y,z,**kw):
    """Plot x,y data with an evolving z component by changing its
    color using a matplotlib colormap 'cm' object.

    Will bomb if z elements are non-finite.
    
    :OPTIONS:
      map : str
        colormap to use for z

      zmin, zmax : floats
         maximum/minimum values of z for colorscale

      sizenotcolor : bool
         If True, 'z' specifies marker size not color.

      others : various
         Any options passable to matplotlib's :func:`plot` 

    :SEE ALSO: :func:`matplotlib.colormaps`

    :REQUIREMENTS: :doc:`pylab`
    """
    # 2010-02-08 16:47 IJC: Created
    # 2011-09-07 10:51 IJMC: Added zmin, zmax options. And use
    #                        variable-number of keywords now.
    # 2013-10-09 17:16 IJMC: Added sizenotcolor option.

    from pylab import plot, cm, array
    
    defaults = dict(map='Blues', zmin=None, zmax=None, linestyle='None', marker='o', sizenotcolor=False)

    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]

    map = kw.pop('map')
    zmin = kw.pop('zmin')
    zmax = kw.pop('zmax')
    sizenotcolor = kw.pop('sizenotcolor')

    try:
        cmap = eval('cm.'+map)
    except:
        print "Colormap %s not found -- exiting!" % map
        return -1

    z2 = array(z,copy=True)
    if not sizenotcolor:
        if zmin is None:
            zmin = z2.min()
        if zmax is None:
            zmax = z2.max()

        z3 = (z2-zmin) / (zmax-zmin)
        z4 = z3 * cmap.N
        z2 = cmap(z4.astype(int))

    plotlist = []
    print kw
    for xx,yy,param in zip(x,y,z2):
        if sizenotcolor:
            kw['ms'] = param
        else:
            kw['color'] = param
        plotlist.append(plot([xx],[yy],**kw))

    return plotlist
    


def hist2d(x,y, bins=None):
    """Compute 2-d histogram data for specified bins.

    :INPUT:
       x

       y

    :OPTIONAL INPUT:
      bins:  a two-tuple containing one of the following:
         (nx,ny) -- tuple, number of bins in each direction

         (xlims, ylims) -- tuple of sequences (x- and y-bin edges)

    :OUTPUT:
       A 3-tuple consisting of:
         xbins -- the centers of the x bins

         ybins -- the centers of the y bins

         hist -- The 2D histogram array

    :SEE ALSO: :func:`numpy.histogram2d`

    :REQUIREMENTS: :doc:`numpy`
    """

    # 2010-02-25 17:26 IJC: Created from my old Hess Diagram fcn.
    # 2010-03-04 13:59 IJC: Fixed a typo -- now actually returns bin centers

    from numpy import arange, array, zeros, isfinite, linspace

    x = array(x).ravel()
    y = array(y).ravel()
    x = x[isfinite(x)]
    y = y[isfinite(y)]

    if bins==None:
        bins = [20,20]
    if hasattr(bins,'__iter__'):
        if len(bins)<2:
            print "bins must have len>2"
            return -1
    else:
        print "bins must have len>2"
        return -1


    # Test X bins
    if hasattr(bins[0],'__iter__'):  # sequence of limits
        nx = len(bins[0])-1
        xlims = bins[0]
    else:
        nx = bins[0]
        xlims = linspace(x.min(), x.max(), nx+1)

    # Test Y bins
    if hasattr(bins[1],'__iter__'):  # sequence of limits
        ny = len(bins[1])-1
        ylims = bins[1]
    else:
        ny = bins[1]
        ylims = linspace(y.min(), y.max(), ny+1)

    xcen = zeros(nx,float)
    ycen = zeros(ny,float)
    hist = zeros((ny,nx),int)
    for ii in range(nx):
        xindex = (x>xlims[ii]) * (x<xlims[ii+1])
        xcen[ii] = 0.5*(xlims[ii]+xlims[ii+1])
        for jj in range(ny):
            ycen[jj] = 0.5*(ylims[jj]+ylims[jj+1])
            yindex = (y>ylims[jj]) * (y<ylims[jj+1])
            index = xindex * yindex
            hist[jj,ii] = index.sum()

    return (xcen, ycen, hist)
    
    

def plotcorrs(params, labs=None, tit=None, xrot=0, yrot=0, cmap=None,figsize=None,plotregion=[.1,.1,.8,.8], n=6, nbins=None, clim=None, docontour=False, contourcolor='k', newfig=True):
    """ Plot correlation coefficient matrix in one big, huge, honkin'
     figure.  Color-code the figure based on the correlation
     coefficients between parameters.

     :INPUTS:
       params -- (M x N array)  M instantiations of a set of N parameters.

     :OPTIONS:
       labs -- (list) labels for each of the N parameters

       tit -- (str) title for figure

       xrot/yrot -- (float) axis label rotation, in degrees

       cmap -- (matplotlib cm) -- colormap for color-coding.

       figsize -- (2-list) -- width and height of figure

       plotregion -- (4-list) -- (left, bottom, width, height) of plotted region in each figure

       n -- (int) -- number of subplots across each figure

       nbins : int
         Bin the data into this many bins, and show 2D histograms instead of points.

       clim : None
         Colorscale limits for normalized 2D histograms (where hist.sum() = 1.0)

       docontour : bool
         Whether to plot contours, or do an 'imshow'

       newfig : bool
         Whether to generate a new figure, or plot in the current axes.

       contourcolor 
         Color of contour line, if "docontour" is set to a list of confidence intervals.

     :REQUIREMENTS: :doc:`pylab`, :doc:`nsdata`

     :NOTES:  
       Based on the concept by Nymyer and Harrington at U. Central Florida

       Beware of plotting two many points, and clogging your system!
    """
    # 2010-05-27 09:27 IJC: Created
    # 2010-08-26 14:49 IJC: Added test for higher-rank dimensional
    # 2010-08-27 09:21 IJC: Moved 'tit' title text lower
    # 2011-11-03 12:03 IJMC: Added 'nbins' option
    # 2012-03-30 09:04 IJMC: Moved axes labels to upper-right, rather than lower-left.
    # 2013-08-19 13:17 IJMC: Added 'docontour' option.
    # 2013-10-09 14:59 IJMC: Added 'newfig' option.

    import pylab as py
    import nsdata as ns
    import kdestats as kde
    n = int(n)

    n, m = params.shape
    if n>=m:
        npts0 = n
        nparam = m
    else:
        npts0 = m
        nparam = n
        params = params.copy().transpose()

    if nbins is not None:
        nbins = int(nbins)
        hist_bins = [py.linspace(min(params[:,ii]), max(params[:,ii]), nbins+1) for ii in range(nparam)]
        hist_cmap = cmap

    nind = params.shape[1]
    if labs is None:
        labs = ['']*nind
    if figsize is None:
        figsize = [9,9]
    nperpage = min(n,nind-1)


    nfigs = py.ceil(float(nind-1.)/nperpage).astype(int)
    #print "nind, nperpage, nfigs>>",nind, nperpage, nfigs
    subx0,suby0, xwid,ywid = plotregion
    subdx = xwid/(nperpage) # was nind-1.
    subdy = ywid/(nperpage) # was nind-1.
    #print "subx0,suby0,subdx,subdy>>",[subx0,suby0,subdx,subdy]


    oldaxlinewidth = py.rcParams['axes.linewidth']
    if nind>40:
        py.rcParams['axes.linewidth'] = 0

    figs = []
    allsubplots = []
    # Iterate over figure columns
    for kk2 in range(nfigs):
        # Iterate over figure rows
        for kk1 in range(nfigs):
            if newfig:
                f=py.figure(nextfig(),figsize)
            else:
                f = py.gcf()
            subplots = []
            jj0 = 0
            #Iterate over subplot columns:
            for jj in range(nperpage*kk2,min(nperpage*(kk2+1),nind)):
                # Set the vertical panel offset:
                if kk1==kk2:  # a figure on the diagonal
                    ii0 = jj0+1
                elif kk1>kk2: # a figure below the diagonal
                    ii0 = 1
                #Iterate over subplots rows:
                for ii in range(max(jj+1,nperpage*kk1+1), min(nperpage*(kk1+1)+1,nind)):
                    #print '(kk2,kk1,jj,jj0,ii,ii0): (%i,%i,%i,%i,%i,%i)'%(kk2,kk1,jj,jj0,ii,ii0), [subx0+subdx*jj0,suby0+subdy*(nperpage-ii0),subdx,subdy]
                    s = py.axes([subx0+subdx*jj0,suby0+subdy*(nperpage-ii0),subdx,subdy])
                    param_doesnt_vary = params[:,jj].std()==0 or params[:,ii].std()==0 or \
                        (py.np.abs(params[:,jj].std()/py.median(params[:,jj])) < 1e-9) or \
                        (py.np.abs(params[:,ii].std()/py.median(params[:,ii])) < 1e-9)
                    if nbins is None or param_doesnt_vary:
                        py.plot(params[:,jj],params[:,ii],',k')
                        #pdb.set_trace()
                    else:
                        #pdb.set_trace()
                        thishist = py.histogram2d(params[:,jj], params[:,ii], \
                                                      bins=[hist_bins[jj], hist_bins[ii]])
                        if docontour:
                            xplot = 0.5*(thishist[1][1:] + thishist[1][0:-1])
                            yplot = 0.5*(thishist[2][1:] + thishist[2][0:-1])
                            if hasattr(docontour, '__iter__'):
                                clev = [kde.confmap(1.0*thishist[0]/npts0, thisDC) for thisDC in docontour]
                                py.contour(xplot, yplot, 1.0*thishist[0].transpose()/npts0, clev, colors=contourcolor, linewidths=2)
                            else:
                                py.contourf(xplot, yplot, 1.0*thishist[0].transpose()/npts0, cmap=hist_cmap)
                            h_axis = py.xlim() + py.ylim()
                        else:
                            ns.imshow(1.0*thishist[0].transpose()/npts0, x=thishist[1], y=thishist[2], cmap=hist_cmap)
                            h_axis = py.xlim() + py.ylim()[::-1]
                        #pdb.set_trace()
                            py.axis(h_axis)
                        if clim is None:
                            py.clim([0., thishist[0].ravel().max()*1.0/npts0])
                        else:
                            py.clim(clim)

                    if jj0>0:  # 
                        s.set_yticklabels('');
                    else:
                        #py.ylabel(labs[ii], rotation=yrot)
                        pass
                    if newfig: s.set_yticks(s.get_yticks()[1:-1]);
                    if ii0 == (jj0+1):
                        s.get_xaxis().set_label_position('top')
                        py.xlabel(labs[jj])
                        s.get_yaxis().set_label_position('right')
                        py.ylabel(labs[ii], rotation='horizontal')
                    if ii0<(nperpage-1) and ii<(nind-1):
                        s.set_xticklabels('');
                    else:
                        #py.xlabel(labs[jj],rotation=xrot)
                        s.get_xaxis().set_major_formatter(py.FormatStrFormatter('%01.2f'));
                    if newfig: s.set_xticks(s.get_xticks()[1:-1]);
                    if nperpage>10:
                        s.set_xticklabels('');
                        s.set_yticklabels('');
                    if nperpage>50:
                        s.set_xticks([])
                        s.set_yticks([])
                    else:
                        [obj.set_rotation(90.) for obj in s.get_xticklabels()] ; 
                    if cmap is not None:
                        s.set_axis_bgcolor(cmap(.3+.7*abs(py.corrcoef(params[:,jj],params[:,ii])[0,1])))
                    #py.title('(kk2,kk1,jj,jj0,ii,ii0): (%i,%i,%i,%i,%i,%i)'%(kk2,kk1,jj,jj0,ii,ii0))
                    if nbins is not None and (not param_doesnt_vary):
                        py.axis(h_axis)
                    subplots.append(s)
                    ii0 += 1

                jj0+=1

            figs.append(f)
            allsubplots.append(subplots)

    if tit is not None:
        f.text(.5,.9,tit,fontsize=24,horizontalalignment='center')


    py.draw()
    py.rcParams['axes.linewidth'] = oldaxlinewidth

    for ff,ss in zip(figs,allsubplots):
        if len(ss)==0:
            py.close(ff)
            

    return figs, allsubplots


def pparams(params, npts=None, figsize=[15,10],newfig=True, labs=None):
    """Take a set of parameters and plot them.  Assume that the larger
    of the array's two dimensions is N (the number of instantiations)
    and the smaller is M (the number of parameters). 

    If npts is not None, then pick only every (N/npts)th point.  
    if newfig is False, plot into the current figure.

    :REQUIREMENTS: :doc:`numpy`, :doc:`pylab`
    """
    # 2012-02-17 09:45 IJMC: Added labs option

    from numpy import sqrt
    from pylab import figure, plot, subplot, xticks, gcf, title

    n, m = params.shape

    if n>=m:
        npts0 = n
        nparam = m
    else:
        npts0 = m
        nparam = n
        params = params.copy().transpose()


    ndown = sqrt(nparam).astype(int)
    nacross = int(1.0*nparam/ndown)
    if ndown*nacross<nparam: 
        nacross +=1
        
    if npts is None:
        npts = npts0
    sampling = max(int(npts0/npts), 1)


    if newfig==True:
        fig = figure(nextfig(), figsize)
    else:
        fig = gcf()

    for ii in range(nparam):
        subplot(ndown,nacross, ii+1)
        plot(params[::sampling,ii]); xticks([])
        if labs is not None:
            title(labs[ii])

    return fig


def hparams(params, nbins=10, figsize=[15,10], normed=False,newfig=True, labs=None, cumulative=False, minorticks=True, plotconf=None, plotmid=False, color=None):
    """Take a set of parameters and histogram them.  Assume that the larger
    of the array's two dimensions is N (the number of instantiations)
    and the smaller is M (the number of parameters). 

    :Options:
      nbins sets the number of bins

      normed sets the histogram normalization

      if newfig is False, plot into the current figure.

      labs is a list of string labels

      cumulative plots normalized cumulative distribution function

      minorticks displays minor tick marks

      plotconf displays confidence levels at the desired (fractional)
          threshold.  E.g., plotconf=0.683 displays the one-sigma
          confidence limits.
      
      plotmid plots nothing (if False), the median (if True), or the
          specified values (if a sequence).  Defaults to median if
          plotconf is True.

      color sets the plotting color.

    :Requirements:  :doc:`numpy`, :doc:`pylab`
    """
    # 2010-11-30 12:11 IJC: Now xticks are rotated 90 degrees.
    # 2012-01-22 16:35 IJMC: Added 'labs' option
    # 2012-03-26 13:09 IJMC: Added 'cumulative' and 'minorticks' options
    # 2012-05-10 16:30 IJMC: Added plotconf and plotmid options.
    # 2013-03-11 08:45 IJMC: Added color option.
    from numpy import sqrt
    from pylab import figure, gcf, subplot, draw
    from analysis import dumbconf

    n, m = params.shape

    if n>=m:
        npts0 = n
        nparam = m
    else:
        npts0 = m
        nparam = n
        params = params.copy().transpose()

    if cumulative:
        normed=True

    if (plotmid is True) or ((plotconf is not None) and (not hasattr(plotmid, '__iter__'))):
        plotmid = np.median(params, axis=0)

    ndown = sqrt(nparam).astype(int)
    nacross = int(1.0*nparam/ndown)
    if ndown*nacross<nparam: 
        nacross +=1
        
    if newfig==True:
        fig = figure(nextfig(), figsize)
    else:
        fig = gcf()
        
    histoptions = dict(histtype='step', cumulative=cumulative, normed=normed)
    if color is not None:
        histoptions['color'] = color

    for ii in range(nparam):
        sax = subplot(ndown,nacross, ii+1)
        sax.hist(params[:,ii],nbins, **histoptions)
        if cumulative:
            sax.set_ylim([0,1])
        else:
            sax.set_ylim([0, sax.get_ylim()[1]])
            sax.set_yticklabels([])
        if minorticks:
           sax. minorticks_on()
        [tick.set_rotation(90.) for tick in sax.get_xaxis().get_ticklabels()]
        if labs is not None:
            sax.text(.05, .9, labs[ii], verticalalignment='top', horizontalalignment='left', transform=sax.transAxes)

        if plotmid is not False:
            sax.plot([plotmid[ii]]*2, sax.get_ylim(), '-k')
        if plotconf is not None:
            x0 = plotmid[ii]
            xup = -x0 + dumbconf(params[:,ii], (1.-plotconf)/2., mid=x0, type='lower')[0]
            xdn =  x0 - dumbconf(params[:,ii], (1.-plotconf)/2., mid=x0, type='upper')[0]
            if xup==0 or xdn==0:
                ndigit = 1
            else:
                ndigit = np.abs(int(np.floor(min(np.log10(np.abs(xup)), np.log10(np.abs(xdn)))) - 1))
            precstr = '%' + ('1.%i' % ndigit) + 'e'
            sax.plot([x0-xdn]*2, sax.get_ylim(), '--k')
            sax.plot([x0+xup]*2, sax.get_ylim(), '--k')
            sax.set_ylabel((precstr + '\n(+' + precstr + ' / -' + precstr + ')\n') % (x0, xup, xdn), horizontalalignment='center')


    draw()

    return fig


def getparams(params, chisq=None, conf=[.683]):
    """Find confidence levels and optimal parameters.

    If chisq is None:
        Take the median of each parameter, and compute the upper and
        lower confidence limits using the parameter distributions and
        this model.
    Else:
        Take a set of fit parameters and associated chi-squared values.
        Take as the optimum model that model with the MEDIAN chi-squared
        value.  Compute upper and lower confidence limits using the
        parameter distributions and this optimum model.

    :INPUTS:
        params : N-D numpy array
            model parameters
        chisq : 1-D numpy array
            chi-squared values.  
        conf : sequence or float
            confidence levels to report

    :OUTPUTS:
    

    :SEE ALSO: 
        :func:`analysis.dumbconf`

    :REQUIREMENTS:  :doc:`numpy`
    """
    # 2010-07-20 11:37 IJC: Created
    # 2010-10-28 12:01 IJC: Updated documentation
    from numpy import sqrt, median, array, ceil,sort,nonzero,floor

    if len(params.shape)==1:
        params = params.reshape(len(params),1)
    n, m = params.shape
    if chisq is None:
        if m>n:
            params = array(params,copy=True).transpose()
        nmod = params.shape[0]
    else:
        nmod = len(chisq)
        if n==nmod:
            pass
        else:
            params = array(params,copy=True).transpose()

    # Don't allow an even number of parameter sets
    if nmod/2==nmod/2.: 
        params = array(params,copy=True)[1::]#[0:nmod-1,:]
        if chisq is None:
            pass
        else:
            chisq = array(chisq,copy=True)[1::]#[0:nmod-1]
        nmod -= 1

    nparam = params.shape[1]

    ret = []
    if chisq is None:
        medmodel = median(params,0)
    else:
        medmodel = params[chisq==median(chisq),:].ravel()
    for ii in range(nparam):
        thisret = [medmodel[ii]]
        for thisconf in conf:
            thisconf = abs(thisconf)
            sorted_param = sort(params[:,ii])
            mid_index = nonzero(sorted_param==medmodel[ii])[0].mean()
            n_offset = nmod*thisconf/2.
            upper_index = ceil(min(mid_index+n_offset,nmod-1)).astype(int)
            lower_index = floor(max(mid_index-n_offset,0)).astype(int)
            upper = sorted_param[upper_index]-medmodel[ii]
            lower = sorted_param[lower_index]-medmodel[ii]
            thisret += [upper, lower]
        ret.append(thisret)

    return ret


def combinations(input_list):
    """Return all possible combinations of the elements in an input
    sequence.  The last returned element will be the empty list.

    E.g., combinations([0,1]) returns [[0, 1], [0], [1], []]

    Taken from the internet:
      http://desk.stinkpot.org:8080/tricks/index.php/2008/04/get-all-possible-combinations-of-a-lists-elements-in-python/

    :Requirements: :doc:`copy`
    """
    # 2010-10-09 17:06 IJC: Created from internet
    import copy

    swap_list_list = [[]]
    for swap in input_list:
        temp_lists = []
        for list_stub in swap_list_list:
            this_list = copy.copy(list_stub)
            this_list.append(swap)
            temp_lists.append(this_list)
            temp_lists.append(list_stub)
        swap_list_list = temp_lists

    return swap_list_list

def pcmodelxcorr(pcaout, data, model, npcs=6, nblock=1000, xl=50, modstr='model', titstr=''):
    """Plot cross-correlations between projection of principal components
    onto data and a model.
    
    :INPUTS:
        pcaout -- output from pcsa.pca, but its only important component
                  is that pcaout[2] is an array of PC eigenvectors; the
                  strongest vector would be pcaout[2][:,-1], etc.

        data -- numpy array.  Should be shape N x M -- N observations of M
                variables, arranged in L blocks, and should have been
                mean-subtracted prior to PCA (i.e., data -= data.mean(0))

        model -- numpy array.  Should be shape M.

        npcs -- int.  Number of principal components to cross-correlate
                with.

        nblock -- int.  NUmber of channels to use at a time in
                  correlations.  Must be an integral divisor of
                  data.shape[1]

        xl -- int.  +/- X-limits to display in plots.

    :Requirements: :doc:`pylab`, :doc:`numpy`, :doc:`pcsa`

    Usage is pretty specific to echelle-data
    """
    # 2010-10-15 11:10 IJC: Created
    import pylab as py
    from numpy import abs
    from pcsa import pca_project
    import pdb

    ps = [pca_project(pcaout[2][:,-ii], data-data.mean(0))[0] \
              for ii in range(1,npcs+1)]

    nobs, ndat = data.shape
    nord = ndat/nblock

    figs = []
    maxlag = nblock/2
    lags = py.arange(-maxlag+1, maxlag)
    for ii in range(npcs):
        figs.append(py.figure())
        x_mod = [py.xcorr(pcp-pcp.mean(), mod-mod.mean(), maxlags=maxlag-1)[1]  \
                     for pcp, mod in zip(ps[ii].reshape(nord,nblock), \
                                             model.reshape(nord, nblock))]
        py.close()
        pdb.set_trace()
        xc_gmean = abs(py.array(x_mod).prod(0))**(1./nord)
        xc_amean = py.array(x_mod).mean(0)

        py.figure()
        py.subplot(211)
        py.semilogy(lags, xc_gmean)
        py.ylabel('geometric mean'  )
        py.xlabel('lags')
        py.xlim([-50,50])
        py.ylim([xc_gmean[abs(lags)<xl].min(), xc_gmean.max()])
        py.title(titstr + '\n' + modstr + ' and PC #%i' % ii)
        py.subplot(212)
        py.plot(lags, xc_amean)
        py.ylabel('arithmetic mean' )
        py.xlabel('lags')
        py.xlim([-xl,xl])

    return figs


def dcf(t, x, y, zerolag=True, nbin=11, binwidth=None, bins=None, prebin=None):
    """Compute the Discrete Correlation Function for unevenly sampled data.

    If your data are evenly sampled, just use :func:`numpy.correlate`!

    :INPUTS:
       t -- (1D sequence) - time sampling of input data.

       x, y -- (1D sequences) - input data.

       Note that t, x, y should all have the same lengths!

    :OPTIONAL INPUT:

       zerolag -- (bool) - whether to compute DCF for zero-lag datapoints.

       nbin -- (int) - number of computed correlation/lag values to average over

       binwidth -- (float) - width of bins to average over, in units
                             of t.  If set, this overrides nbin.

       bins -- (sequence) - edges of bins to use in averaging DCF
                            values.  If set, this overrides both nbin
                            and binwidth.
                            
       prebin -- (scalar) - factor by which to bin down initial data
                            and lag pairs.  This translates into a
                            speed boost of about this same factor.
                            

    :RETURNS:
       meanlags -- average lag in each averaged bin

       rdcf -- DCF value in each averaged bin

       rdcf_err -- uncertainty in the DCF value for each averaged bin

    :SEE ALSO:  :func:`numpy.correlate`
    
    :REQUIREMENTS: :doc:`analysis`, :doc:`numpy`
    
       """
    # 2010-11-12 IJC: Created
    # 2010-11-17 09:25 IJC: Added 'bins' option
    # 2011-03-22 13:54 IJC: Added 'prebin' option
    # 2012-03-21 13:25 IJMC: Switched "import *" to "import array, " ... etc.

    import analysis as an
    #from numpy import array, meshgrid, nonzero, arange, argsort, isfinite
    import pdb

    t = np.array(t)
    x = np.array(x, copy=True)
    y = np.array(y, copy=True)
    nx, ny = x.size, y.size
    npts = max(nx, ny)
    
    x -= x.mean()
    y -= y.mean()
    xx, yy = np.meshgrid(x, y)

    sx = x.std()
    sy = y.std()

    # Generate data-pair indices:
    if zerolag:
        pairs1, pairs2 = np.nonzero(np.ones((npts, npts), bool))
    else:
        xind, yind = np.meshgrid(np.arange(npts), np.arange(npts))
        pairs1, pairs2 = np.nonzero(xind<>yind)
        del xind
        del yind

    # Compute lags:
    tt, tt2 = np.meshgrid(t, t)
    lags = (tt-tt2)[pairs1, pairs2]
    del tt
    del tt2

    uij = (xx * yy)[pairs1, pairs2] / (sx * sy)
    del xx
    del yy

    tind = np.argsort(lags)

    lags = lags[tind]
    uij = uij[tind]
    del tind

    #pdb.set_trace()

    # The regular "DCF" uncertainty is just the standard deviation of the mean:
    #if bins is not None:
    #    meanlags, rdcf, meanlags_width, rdcf_err = \
    #              errxy(lags, uij, bins, xerr='minmax', yerr='sdom')
    #elif binwidth is None:
    #    meanlags, rdcf, meanlags_width, rdcf_err = \
    #              errxy(lags, uij, None, xerr='minmax', yerr='sdom', binfactor=nbin)
    #else:
    #    bins = arange(lags.min(), lags.max() + binwidth, binwidth)
    #    meanlags, rdcf, meanlags_width, rdcf_err = \
    #              errxy(lags, uij, bins, xerr='minmax', yerr='sdom')

    if prebin>1:
        lags = an.binarray(lags, 2)
        uij = an.binarray(uij, 2)

    if bins is not None:
        meanlags, rdcf, meanlags_width, rdcf_err = \
                  errxy(lags, uij, bins, xerr=None, yerr=None)
    elif binwidth is None:
        meanlags, rdcf, meanlags_width, rdcf_err = \
                  errxy(lags, uij, None, xerr=None, yerr=None, binfactor=nbin)
    else:
        bins = np.arange(lags.min(), lags.max() + binwidth, binwidth)
        meanlags, rdcf, meanlags_width, rdcf_err = \
                  errxy(lags, uij, bins, xerr=None, yerr=None)
        

    finite_ind = np.isfinite(meanlags) * np.isfinite(rdcf) #* \
#                   isfinite(rdcf_err)
    meanlags = meanlags[finite_ind]
    rdcf = rdcf[finite_ind]
    #rdcf_err = rdcf_err[finite_ind]

    return meanlags, rdcf, rdcf_err



def getfilelist(path='.', includes=[], excludes=[]):
    """Return a list of filenames meeting certain criteria.

    :INPUTS:
      path -- (str) path to directory to be scanned

      includes -- (list of strs) -- all strings in this list must be
                     present in a filename to be returned

      excludes -- (list of strs) -- any string in this list will
                     prevent a filename from being returned

    """
    # 2011-01-30 22:20 IJC: Created

    import os

    f = os.popen('ls %s/' % path)
    filenames = [os.path.split(line)[1].strip() for line in f.readlines()]
    f.close()

    filtered_filenames = []
    if len(includes)>0:
        for fn in filenames:
            string_found = True
            for incstr in includes:
                if fn.find(incstr)==-1:
                    string_found = False
            if string_found:
                filtered_filenames.append(fn)

    if len(excludes)>0:
        returned_filenames = []
        for fn in filtered_filenames:
            file_excluded = False
            for excstr in excludes:
                if fn.find(excstr)>-1:
                    file_excluded = True
            if not file_excluded:
                returned_filenames.append(fn)
    else:
        returned_filenames = filtered_filenames

    return returned_filenames

def loadpickle(filename):
    """Load a pickle from a given filename.  If it can't be loaded by
    pickle, return -1 -- otherwise return the pickled object.

    E.g., 
       data = tools.loadpickle(filename)"""
    # 2011-02-10 15:45 IJC: Created
    import pickle

    good = True

    try:
        f = open(filename, 'r')
    except:
        print "Could not open file: %s" % filename
        good = False

    if good:
        try:
            ret = pickle.load(f)
        except:
            print "Could not load pickle from %s" % filename
            good = False

    try:
        f.close()
    except:
        print "Could not close file %s" % filename
        good = False

    if good:
        pass
    else:
        ret = -1

    return ret
    

def savepickle(obj, filename):
    """Save a pickle to a given filename.  If it can't be saved by
    pickle, return -1 -- otherwise return the file object.

    To save multiple objects in one file, use (e.g.) a dict:

       tools.savepickle(dict(a=[1,2], b='eggs'), filename)
       """
    # 2011-05-21 11:22 IJMC: Created from loadpickle.
    # 2011-05-28 09:36 IJMC: Added dict example
    import pickle

    good = True

    try:
        f = open(filename, 'wb')
    except:
        print "Could not open file: %s : for writing." % filename
        good = False

    if good:
        try:
            ret = pickle.dump(obj, f)
        except:
            print "Could not write object to pickle file: %s" % filename
            good = False

    try:
        f.close()
    except:
        print "Could not close pickle file %s" % filename
        good = False

    if good:
        pass
    else:
        f = -1

    return f
    

def dict2obj(dic):
    """Take an input Dict, and turn it into an object with fields
    corresponding to the dict's keys."""
    # 2011-02-17 09:41 IJC: Created

    from spitzer import baseObject
    
    ret = baseObject()

    if not isinstance(dic, dict):
        print "Input was not a Python dict!  Exiting."
    else:
        for key in dic.keys():
            exec('ret.%s = dic["%s"]' % (key, key))

    return ret

        


#def loadspectrum(filename, lambdascale=1., limits=None, skiprows=None, lamdacol=0, datacol=1 ):
#    """Load a spectrum from a FITS or ASCII file, and return the
#    wavelength and data in two vectors.#
#
#    :INPUTS:
#
#       filename (str) -- filename to load.  If extension is "FITS"
#                         load using pyfits, otherwise try to load as
#                         an ASCII file.
#
#       lambdacol (int) -- column number



def find_peaks(vec, sep=0, thresh=None):
    """
    Find all large values in input vector that are separated by at least
    'wid' pixels.  

    :INPUTS:

       vec (sequence) -- 1D vector

       sep (scalar) -- minimum separation of returned peaks 

       thresh (scalar) -- ignore all peaks lower than this value.

    :EXAMPLE:

       import pylab as py
       import tools
       x = py.linspace(0, 10, 100)   # Generate fake time index
       y = py.sin(6.28*x/10.) + py.sin(6.28*x/2.)   # Generate fake data
       peakvals, peaklocs = tools.find_peaks(y, sep=10)   # Find the peaks
       py.plot(x, y, '-', x[peaklocs], peakvals, 'or')   # Plot them

    :RETURNS:

       peakvals, peakindices
    """
    # 2011-03-22 14:54 IJC: Created
    # 2012-08-09 22:51 IJMC: Added thresh option.

    import numpy as np
    from pylab import find
    import pdb

    if thresh is None:
        thresh = -np.inf

    vec = np.array(vec)
    npix = len(vec)
    sep = np.floor(sep).astype(int)

    available_index = np.ones(npix, bool)
    peakvals    = np.zeros(npix, float)
    peakindices = np.zeros(npix, float)
    npks = 0
    inds = np.arange(npix, dtype=int)

    while available_index.any():
        #pdb.set_trace()
        this_peak = vec[available_index].max()
        this_location = inds[available_index][vec[inds[available_index]] == vec[inds[available_index]].max()][0] 
        available_index[max(0, this_location - sep):min(npix, this_location + sep)+1] = False
        if this_peak >= thresh:
            peakvals[npks]    = this_peak
            peakindices[npks] = this_location
            npks += 1
        
    peakvals = peakvals[0:npks]
    peakindices = peakindices[0:npks].astype(int)

    return peakvals, peakindices

def addobj(obj1, obj2, exclude='TBD'):
    """Combine fields in two objects with the same attributes.  A
    handy function!

    :INPUTS:
       obj1, obj2 : objects of the same kind

    :RETURNS:
       obj12 : new object, with fields of 1 and 2 combined.
         
         OR

       -1, if the objects have absolutely no attributes in common.

    :NOTES:
       Methods/attributes beginning with an underscore (_) are not combined.
       """
    # 2011-05-25 10:11 IJMC: Created
    # 2011-06-01 15:15 IJMC: Fixed a bug in concatenate(common_dim)

    import copy
    import numpy

    # Helper function:
    def combine_atts(a1, a2):
        """Combine two objects, as best you can!  

        If they are of different classes, then return None.

        If they are non-iterable, combined them in a 2-tuple."""
        combined = False
        if a1.__class__ <> a2.__class__:
            ret = None
            combined = True

        if ((not hasattr(a1, '__iter__')) or isinstance(a1, list) or isinstance(a1, str)) \
                and (not combined):
            try:
                ret = a1 + a2
            except:
                ret = (a1, a2)
            combined = True

        if isinstance(a1, numpy.ndarray) and (not combined):
            sh1, sh2 = a1.shape, a2.shape
            if sh1 == sh2:  # Exactly same dimensions; I don't know
                            # which dimension to combine along, so
                            # just create a new one.
                ret = numpy.array([a1, a2])
                combined = True
            elif len(sh1) <> len(sh2): # Wholly different shapes;
                                       # can't combine.
                ret = (a1, a2)
                combined = True
            elif len(set(sh1).difference(sh2)) > 1: # Too many disparate dimensions.
                ret = (a1, a2)
                combined = True
            else: # There must be exactly 1 non-common dimension.
                if len(sh1)==0 or len(sh1)==1:  # Scalar (0D) or Vector (1D):
                    ret = numpy.concatenate((a1.ravel(), a2.ravel()))
                    combined = True
                else:  #  Array (>1 D)
                    common_dim = (numpy.arange(len(sh1)) * \
                         (True - (numpy.array(sh1) == set(sh1).intersection(sh2).pop()))).sum()
                    ret = numpy.concatenate((a1, a2), common_dim)
                    #if common_dim==0:
                    #    ret = numpy.hstack((a1, a2))
                    #    combined = True
                    #elif common_dim==1:
                    #    ret = numpy.vstack((a1, a2))
                    #    combined = True
                    combined = True

        if not combined:
            ret = (a1, a2)

        return ret

    try:
        newobj = copy.deepcopy(obj1)
    except:
        try:
            newobj = copy.deepcopy(obj2)
        except:
            print "copy.deepcopy() failed... could not copy either input object"
            return -1

    # Get all attributes:
    att1all = dir(obj1)
    att2all = dir(obj2)
    
    # Exclude everything beginning with an underscore:
    att1 = []
    att2 = []
    for att in att1all:
        if att[0] == '_':
            pass
        else:
            att1.append(att)
    for att in att2all:
        if att[0] == '_':
            pass
        else:
            att2.append(att)
    
    # Get attributes common to both objects:
    common_atts = set(att1).intersection(att2)
    abandoned_atts = (set(att1).union(att2)).difference(common_atts)

    if len(common_atts)==0:
        return -1

    # Remove all non-common attributes:
    for att in abandoned_atts:
        if hasattr(newobj, att):
            setattr(newobj, att, None)

    
    # Combine all common attributes:
    for att in common_atts:
        if (not hasattr(newobj, att)) or (not hasattr(obj1, att)) or \
                (not hasattr(obj2, att)): # Test for coding mistakes:
             print "Couldn't find attribute '%s'; something went wrong!" % att
        else:
            setattr(newobj, att, combine_atts(getattr(obj1, att), getattr(obj2, att)))

    return newobj


def multifunc(params, func1, func2, npar1, npar2=None, args1=(), args2=()):
    """Multiply two functions together.

    :EXAMPLE:
      ::
    
       import numpy as np
       multifunc([.785, .785], np.cos, np.sin, 1)  # returns cos(pi/4) * sin(pi/4) ~ 0.5

    :EXAMPLE:
      ::

       import pylab as py
       import tools
       x = 2*py.pi * py.linspace(-5, 5, 500)
       y = tools.multifunc([x/8., x], np.sin, np.sin, 1).ravel()
       py.plot(x, y)
       """
    # 2011-06-03 11:26 IJC: Created
    # 2012-06-10 16:08 IJMC: Now reflects preferred use of multifunc_general.

    print "Multifunc() is deprecated; use multifunc_general() instead."

    # Set the parameters for each function:
    if npar2 is None:
        npar2 = len(params) - npar1
    return multifunc_general(params, (func1, func2), (npar1, npar2), (args1, args2))


def multifunc_general(params, funcs, npar=None, args=None):
    """Multiply results of several functions together.

    :INPUTS:
      params : sequence
        Concatenated sequence of parameters (first arguments) for each
        of the functions in 'funcs'.

      funcs : function or tuple
        Single function, or a tuple of several functions to call.  

      npar : tuple
        If more than one function is used, npar must be a tuple
        specifying the number of parameters passes to each function
        (as its first input).  E.g., if npar = (2, 3) then the result
        will be: funcs[0](params[0:2]) * funcs[1](params[2:5])
       
      args : tuple
        If more than one function is used, args must be a tuple
        specifying the additional arguments to be passed to each
        function. E.g.: funcs[0](params[0:2], *args[0])

    :EXAMPLE:
      ::
    
       import numpy as np
       multifunc_general([.785, .785], (np.cos, np.sin), (1, 1))  # returns cos(pi/4) * sin(pi/4) ~ 0.5

    :EXAMPLE:
      ::

       import pylab as py
       import tools
       x = 2*py.pi * py.linspace(-10, 10, 500)
       y = tools.multifunc_general([x/8., x], (np.sin, np.sin), (1, 1)).ravel()
       py.plot(x, y)
       """
    # 2012-06-10 15:54 IJMC: Created; modified from multifunc


    if not hasattr(funcs, '__iter__'):
        ret = funcs(params, *args)
    else:
        ret = 1.
        nfunc = len(funcs)
        if npar is None:
            print "Multiple functions input; npar ought not to be None!"
            return -1

        if args is None:
            args = ((),) * nfunc

        i0 = 0
        for ii in range(nfunc):
            these_params = params[i0:i0+npar[ii]]
            i0 += npar[ii]
            ret *= funcs[ii](these_params, *args[ii])

    return ret


def sumfunc_general(params, funcs, npar=None, args=None):
    """Add results of several functions together.

    :INPUTS:
      params : sequence
        Concatenated sequence of parameters (first arguments) for each
        of the functions in 'funcs'.

      funcs : function or tuple
        Single function, or a tuple of several functions to call.  

      npar : tuple
        If more than one function is used, npar must be a tuple
        specifying the number of parameters passes to each function
        (as its first input).  E.g., if npar = (2, 3) then the result
        will be: funcs[0](params[0:2]) * funcs[1](params[2:5])
       
      args : tuple
        If more than one function is used, args must be a tuple
        specifying the additional arguments to be passed to each
        function. E.g.: funcs[0](params[0:2], *args[0])

    :EXAMPLE:
      ::
    
       import numpy as np
       sumfunc_general([.785, .785], (np.cos, np.sin), (1, 1))  

    :EXAMPLE:
      ::

       import pylab as py
       import tools
       x = 2*py.pi * py.linspace(-10, 10, 500)
       y = tools.sumfunc_general([x/8., x], (np.sin, np.sin), (1, 1)).ravel()
       py.plot(x, y)
       """
    # 2012-06-10 15:54 IJMC: Created; modified from multifunc_general


    if not hasattr(funcs, '__iter__'):
        ret = funcs(params, *args)
    else:
        ret = 1.
        nfunc = len(funcs)
        if npar is None:
            print "Multiple functions input; npar ought not to be None!"
            return -1

        if args is None:
            args = ((),) * nfunc

        i0 = 0
        for ii in range(nfunc):
            these_params = params[i0:i0+npar[ii]]
            i0 += npar[ii]
            ret += funcs[ii](these_params, *args[ii])

    return ret



def sumfunc(params, func1, func2, npar1, npar2=None, args1=(), args2=()):
    """Add two functions together.

    :EXAMPLE:
      ::
    
       import numpy as np
       sumfunc([.785, .785], np.cos, np.sin, 1)  # returns cos(pi/4) + sin(pi/4) ~ XXX

    :EXAMPLE:
      ::

       import pylab as py
       import tools
       x = 2*py.pi * py.linspace(-5, 5, 500)
       y = tools.sumfunc([x/8., x], np.sin, np.sin, 1).ravel()
       py.plot(x, y)
       """
    # 2011-10-27 11:50 IJMC: Created from multifunc

    # Set the parameters for each function:
    if npar2 is None:
        npar2 = len(params) - npar1
    params1 = params[0:npar1]
    params2 = params[npar1:(npar1 + npar2)]

    return func1(params1, *args1) + func2(params2, *args2)


def imrot(*varargin):
    """
     Rotate a SQUARE image by theta (in radians)

     :SYNTAX:
       res = imrot (img, theta)
    
    
    :NOTES: 
       This is VERY SLOW and only uses nearest-neighbor interpolation.

       Image must be square; size remains the same

    """
    # 2011-06-08 15:48 IJMC: Adapted to Python, but without interp2.
    # 2006/11/30 IJC: Added interpolation option
    # Written by Joseph Green at the Jet Propulsion Laboratory
    
    from analysis import pad
    from pylab import array, arange, meshgrid, cos, sin, isnan, dot, sqrt, find, np, dot, vstack
    import pdb


    nargin = len(varargin)

    if nargin==1:
        return varargin[0]
    elif nargin==2:
        img = varargin[0]
        theta = varargin[1]
        interpMethod = 'nn'
    else:
        img = varargin[0]
        theta = varargin[1]
        interpMethod = varargin[2];

    npix = max(img.shape);
    img  = pad(img,npix);
    xs = arange(-npix/2, npix/2)  + 0.

    x, y = meshgrid(xs,xs)

    M = array([[cos(theta), sin(theta)], [-sin(theta), cos(theta)]])

    xr=x.copy()
    yr=y.copy()

    xxr, yyr = dot(M, vstack((x.ravel(), y.ravel())))
    #pdb.set_trace()

    res = 0*img
    xlim = xs.min(), xs.max()

    for m in range(npix):
        for n in range(npix):
            if xr[m,n] < xlim[0] or xr[m,n] > xlim[1] or \
                    yr[m,n] < xlim[0] or yr[m,n] > xlim[1]:
                res[m,n] = 0.
            else:
                distances = sqrt((x - xr[m,n])**2 + (y - yr[m,n])**2).ravel()
                #distances = np.abs((x - xr[m,n]) + 1j*(y - yr[m,n]))
                mindist = distances.min()
                #locations = find(distances==mindist)
                if mindist==0:
                    res[m,n] = (img.ravel()[distances==mindist]) #.mean()
                else:
                    res[m,n] = 0.

                
    return res


def readPSplotdata(filename, spacedelimiter=' ', pairdelimiter=','):
    """Read in a raw PostScript plot datafile and output Numpy arrays.

    :INPUTS:
      filename : str
         filename of the data file.  

      spacedelimiter : str
         delimiter between data pairs; by default a space (" ")

      pairdelimiter : str
         delimiter within data pairs; by default a comma (",")

    :FILE_FORMAT:
      The datafile should have the format:
       "m x0,y0 dx1,dy1 dx2,dy2 ..."

    :OUTPUTS:
      (x, y) -- tuple of 1D arrays
      """
    # 2011-08-24 09:20 IJMC: Created
    from numpy import array

    # Read file:
    f = open(filename, 'r')
    rawdata = f.readlines()
    f.close()

    # Check for validity

    # Get datapoints
    datavals = []
    for line in rawdata:
        datavals = datavals + line.split(spacedelimiter)
        
    # 
    foundinit = False
    for dataval in datavals:
        if dataval.find(pairdelimiter)>-1:
            xval, yval = map(float, dataval.split(pairdelimiter))
            if not foundinit: # starting point
                foundinit = True
                x, y = [xval], [yval]
            else: # differential 
                x.append(x[-1] + xval)
                y.append(y[-1] + yval)

    return array(x), array(y)


def bnu(T, lam):
    """Planck function in frequency.

    :INPUTS:
      T : scalar or array
        temperature in Kelvin

      lam : scalar or array
        wavelength in microns [but intensity will be per Hz]

    Value returned is in cgs units: erg/s/cm^2/Hz/sr
    """
    # 2011-11-04 10:47 IJMC: Added to Tools
    # 2013-02-19 19:42 IJMC: Updated to use precise value of c.
    from numpy import exp

    c = 29979245800. # cm/s
    nu = c/(lam/1e4)
    h = 6.626068e-27 # cgs units
    k =  1.3806503e-16 # cgs units
    expo = h*nu/(k*T)
    nuoverc = 1./ (lam/1e4)
    return ((2*h*nuoverc**2 * nu)) / (exp(expo)-1)

def blam(T, lam):
    """Planck function in wavelength.

    :INPUTS:
      T : scalar or array
        temperature in Kelvin

      lam : scalar or array
        wavelength in microns

    Value returned is in (nearly) cgs units: erg/s/cm^2/micron/sr
    """
    # 2012-06-12 19:56 IJMC: Created.
    # 2013-02-19 19:42 IJMC: Updated to use precise value of c.

    c = 29979245800. # cm/s
    nu = c/(lam/1e4) # Hz

    # Take advantage of the fact that (lam F_lam) = (nu F_nu):

    return bnu(T, lam) * (nu / lam) 


def planet2temp(rprs, fpfs, teff, lam=24., gamma=0.8, ntrials=10000):
    """Convert planet/star size and flux contrasts to brightness temperatures.

    :INPUTS:
      rprs : 2-tuple, either
        [0] : (value, dist) where dist is the sampled posterior
              distribution of planet/star radius ratio (e.g., from
              MCMC output).  OR:
        [1] : (value, uncertainty)

      fpfs : 2-tuple, either
        [0] : (value, dist) where dist is the sampled posterior
              distribution of planet/star flux ratio (e.g., from
              MCMC output).  OR:
        [1] : (value, uncertainty)

      teff : 2-tuple, either
        [0] : (value, dist) where dist is the sampled posterior
              distribution of stellar effective temperatures (e.g., from
              MCMC output).  OR:
        [1] : (value, uncertainty)

      lam : float
        wavelength of observations, in microns

      gamma : float
        factor to account for the fact that a star's infrared flux is
        lower than predicted by a blackbody of a given effective
        temperature.  Set to 0.8 for 24 microns.

      ntrials : int
        Number of times to resample distributions.

     :REQUIREMENTS:
       :doc:`scipy.optimize`, :doc:`numpy`
       """
    # 2011-11-04 11:27 IJMC: Created

    from numpy import random
    from scipy import optimize, array

    def err_temp_bb(temp, bnu_t, lam):
        return bnu_t - bnu(temp, lam)

    if not hasattr(rprs, '__iter__'):
        rprs = (rprs, 0)
    if not hasattr(fpfs, '__iter__'):
        fpfs = (fpfs, 0)
    if not hasattr(teff, '__iter__'):
        teff = (teff, 0)

    # Create resampled distributions:
    npts_r = len(rprs)
    npts_f = len(fpfs)
    npts_t = len(teff)
    resample_r = hasattr(rprs[1], '__iter__')
    resample_f = hasattr(fpfs[1], '__iter__')
    resample_t = hasattr(teff[1], '__iter__')
    rprs_0 = rprs[0]
    fpfs_0 = fpfs[0]
    teff_0 = teff[0]

    if not resample_r:
        rprs = random.normal(rprs[0], rprs[1], ntrials)
    else:
        rprs = array(rprs[1], copy=False)[random.uniform(0, len(rprs[1]), ntrials).astype(int)]
    if not resample_f:
        fpfs = random.normal(fpfs[0], fpfs[1], ntrials)
    else:
        fpfs = array(fpfs[1], copy=False)[random.uniform(0, len(fpfs[1]), ntrials).astype(int)]
    if not resample_t:
        teff = random.normal(teff[0], teff[1], ntrials)
    else:
        teff = array(teff[1], copy=False)[random.uniform(0, len(teff[1]), ntrials).astype(int)]

    planet_flux  = fpfs * bnu(teff*gamma, lam) / (rprs**2)
    planet_flux_0  = fpfs_0 * bnu(teff_0*gamma, lam) / (rprs_0**2)

    planet_temps = array([optimize.fsolve(err_temp_bb, 1000, \
                                              args=(thisflux, lam)) \
                              for thisflux in planet_flux]).ravel()
    planet_temp_0 = optimize.fsolve(err_temp_bb, 1000, \
                                        args=(planet_flux_0, lam))

    return (planet_temp_0, planet_temps)


def erf_approx(z, N):
    """Weideman 1994's approximate complex error function.

    :INPUTS:
      z : real or complex float

      N : number of terms to use.

    :NOTES:
      returns w(z) = exp(-z**2) erfc(-1j*z)
      """
    # 2011-11-14 17:01 IJMC: Created

    import numpy as np

    M = 2 * N
    M2 = 2*M # number of sampling points
    k = np.arange(-M+1, M).reshape(M2 - 1, 1)
    L = np.sqrt(N / np.sqrt(2))
    
    theta = k * np.pi/M
    t = L * np.tan(theta/2)  

    # Function to be transformed:
    f = np.exp(-t**2) * (L**2 + t**2)
    f = np.concatenate(([0], f))


def ee_psf(psf, energy=0.8, center=None):
    """%% Determine the diameter in pixels for a given Encircled Energy
    %
    % [d, ee] = ee_psf(psf, energy)
    % [d, ee] = ee_psf(psf, energy, center);
    %
    % INPUTS: psf   - array representing a point spread function
    %         energy- encircled energy percentage (default = 0.8).
    %                   Can also be a vector of values, in which case the
    %                   outputs 'd' and 'encircledenergy' are also vectors.
    % OPTIONAL INPUT:
    %         center- [x y] coordinates of the center of the psf.  If not
    %                   passed, defaults to the coordinates of the maximum-
    %                   valued point of the psf array.
    %
    % OUTPUTS: d    - diameter of circle enclosing 'energy' [pix] at the
    %                   corresponding desired encircled energy value
    %          ee   - encircled energy at computed 'd' at the
    %                   corresponding desired encircled energy value
    %function [d, encircledenergy] = ee_psf(psf, energy, center)
    """
    # 2012-01-06 15:48 IJMC: Converted from Matlab to Python;
    #                        originally written at JPL.

    from analysis import pad
    import numpy as np

    psf = np.array(psf, copy=False)
    npix = max(psf.shape)
    psf = pad(psf,npix);

    if center is None:
         xc, yc = np.nonzero(psf==psf.ravel().max())
         center = [yc, xc];

    if energy is None:
        energy = np.array([0.8])
    elif not hasattr(energy, '__iter__'):
        energy = np.array([energy])
    else:
        energy = np.array(energy, copy=False)

    encircledenergy = np.zeros(len(energy))  # initialize encircled energy array
    d  = np.zeros(len(energy))  # initialize diameter-pixel array

    xc, yc = center[0:2]
    ee = 0
    eeold = -1     
    xs = np.arange(npix) - xc
    ys = np.arange(npix) - yc

    [x,y] = np.meshgrid(xs, ys)
    r = np.abs(x + 1j*y)
    energy_tot = psf.sum()

    for inc_energy in range(len(energy)):
         rad = 0
         ee = 0
         dfactor = 4       # these two lines determine how it searches
         dpix = dfactor**np.floor(np.log(npix)/np.log(dfactor) - 1)

         while  (ee <> eeold):
             while (ee < energy[inc_energy]) and (ee<>eeold):
                 rad = rad+dpix
                 emasked = (r < rad)*psf
                 eeold = ee
                 ee = emasked.sum()/energy_tot

             #disp(['Inner: Rad=' num2str(rad) ' ,Encircled energy = ' num2str(ee)])

             rad = rad - dpix
             dpix = dpix/dfactor
             rad = rad + dpix
             emasked = (r < rad)*psf
             eeold = ee
             ee = emasked.sum()/energy_tot
           #disp(['Outer: Rad=' num2str(rad) ', encircled energy = ' num2str(ee)])

         d[inc_energy] = 2*rad
         encircledenergy[inc_energy] = ee

    return d, encircledenergy



def extinct_cardelli(lam_um, RV=3.1, warn=True):
    """Compute the Cardelli et al. 1989 A_lam/A_V extinction.

    :INTPUTS:
      lam_um : float or Numpy array
        wavelength desired, in microns

      RV : float
        R_V extinction parameter

      warn : bool
        If True, print a warning if wavelength is outside the valid range.

    :OUTPUTS:
      extinct : float or Numpy array
        extinction -- A_lambda/A_V

    :NOTES:
      This is only valid for wavelengths in the range 0.3 - 3.3 microns!
    """
    # 2012-02-24 23:06 IJMC: Created
    from numpy import array

    if not hasattr(lam_um, '__iter__'):
        lam_um = array([lam_um], copy=True)
    else:
        lam_um = array(lam_um, copy=True)

    
    ind0 = lam_um < (1./1.1)  # visible
    ind1 = lam_um >= (1./1.1)  # near-infrared
    ind2 = (lam_um < 0.25) + (lam_um > 3.5) 
    extinct = 0*lam_um

    if warn and ind2.any():
        print "Wavelengths found outside valid range (0.3-3 microns).  Results are suspect."

    x = 1./lam_um
    if ind0.any():
        y = x[ind0] - 1.82
        extinct[ind0] += 1.
        a = [0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999]
        b = [1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260, -2.09002]
        ypow = y.copy()
        for ii in range(7):
            extinct[ind0] += (a[ii] + b[ii] / RV) * ypow
            ypow *= y

    if ind1.any():
        extinct[ind1] = x[ind1]**1.61 * (0.574 - 0.527 / RV)
    
    return extinct


def shift_image(im, dx=0, dy=0, buffer=0):
    """Shift a 2D image by the specified number of integer pixels.
    
    :INPUTS:
      im : 2D Numpy array
        numpy array of desired size

      dx : int
        number of pixels to shift in x direction

      dy : int
        number of pixels to shift in y direction

      buffer : scalar
        Default value for unfilled pixels

        """
    # 2012-02-25 02:46 IJMC: Created
    
    from numpy import array, abs

    im = array(im, copy=False)

    nx, ny = im.shape
    newim = 0*im + buffer
    if abs(dx)<nx and abs(dy)<ny:
        xnew1 = max(0, dx)
        xnew2 = min(nx, nx+dx)
        xorg1 = max(0, -dx)
        xorg2 = min(nx, nx-dx)
        ynew1 = max(0, dy)
        ynew2 = min(ny, ny+dy)
        yorg1 = max(0, -dy)
        yorg2 = min(ny, ny-dy)
        newim[xnew1:xnew2,ynew1:ynew2] = im[xorg1:xorg2,yorg1:yorg2]

    return newim


def resamp(frame, resamp, retcoords=False):
    """Resample a 2D array by a given factor, using bilinear interpolation.

    :INPUTS:
      frame : 2D Numpy Array
        data array to be resampled

      resamp : positive scalar
        resampling factor (typically an integer greater than 1)

      retcoords: bool
        If True, return tuple (frame2, x2, y2)

     :NOTES:
       Result needs to be scaled by (1/resamp^2) to conserve flux
        """
    # 2012-02-25 07:21 IJMC: Created
    # 2012-02-26 14:19 IJMC: Added retcoords option

    from numpy import array, arange
    from scipy import interpolate

    # Parse inputs:
    resamp = float(resamp)
    frame = array(frame, copy=False)
    nx0, ny0 = frame.shape

    nx = ((nx0 - 1)*resamp + 1.)  # Avoid resampling at pixel locations
    ny = ((ny0 - 1)*resamp + 1.)  #   outside the original boundaries.
       
    xx0 = range(nx0)
    yy0 = range(ny0)
    x1,y1 = arange(nx)/resamp, arange(ny)/resamp
    rectspline = interpolate.fitpack2.RectBivariateSpline(xx0, yy0, frame, kx=1, ky=1, s=0)
    frame2 = rectspline(x1, y1)#/resamp/resamp

    if retcoords:
        ret = frame2, x1, y1
    else:
        ret = frame2

    return ret



def plotlikelihood_2d(L, x=None, y=None, conf=[0.683], figsize=[8, 6], contourargs=dict(colors='k'), posteriorargs=dict(color='k'), limitargs=dict(color='k', linestyle=':'), xlabel=None, ylabel=None, buffers=[.1, .1, .1, .1]):
    """Plot contours and histograms for 2D likelihood array.

    :INPUTS:
      L : 2d Numpy array

        Likelihood values, not necessarily normalized.  (Remember that
        L = exp[-chi^2 / 2.] )

      x : sequence
        Values along the first dimension of L.  Thus len(x) must equal
        L.size[0]

      y : sequence
        Values along the second dimension of L.  Thus len(x) must equal
        L.size[1]

      figsize : 2-sequence
        Size of figure to be created, in inches.

      conf : scalar or sequence
        Confidence intervals to plot

      contourargs : dict
        Keyword arguments to be passed to matplotlib.contour

      posteriorargs : dict
        Keyword arguments to be passed to matplotlib.plot for
        posterior distributions

      limitargs : dict
        Keyword arguments to be passed to matplotlib.plot for
        1D confidence limits

      xlabel : str
      
      ylabel : str

      buffers : 4-sequence
        fractional buffer width around edges: [left, right, bottom, top]
    """
    # 2012-06-19 11:47 IJMC: Created
    from kdestats import confmap
    import pylab as py

    if not hasattr(conf, '__iter__'):
        conf = [conf]

    nconf = len(conf)

    # Define array indices:
    nx, ny = L.shape
    x0 = np.arange(nx)
    y0 = np.arange(ny)

    # Define axis values:
    if x is None:
        x = x0.copy()
    if y is None:
        y = y0.copy()

    # Define axis labels:
    if xlabel is None:
        xlabel = ''
    if ylabel is None:
        ylabel = ''

    # Define marginalized posterior probability distributions:
    yppd = L.sum(0)
    xppd = L.sum(1)

    # Define confidence levels:
    conflevels = confmap(L, conf)
    xconf = [[np.interp(val, np.cumsum(xppd/xppd.sum()), x) for val in [(1.-cl)/2., 1. - (1.-cl)/2.]] for cl in conf]
    yconf = [[np.interp(val, np.cumsum(yppd/yppd.sum()), y) for val in [(1.-cl)/2., 1. - (1.-cl)/2.]] for cl in conf]


    # Prepare for plotting:
    lbuffer, rbuffer, bbuffer, tbuffer = buffers


    eff_width = 1. - lbuffer - rbuffer
    eff_height = 1. - tbuffer - bbuffer
    contour_frac = 0.6

    # Begin plotting:
    py.figure(nextfig(), figsize)
    cax = py.subplot(2, 2, 3, position=[lbuffer, bbuffer, eff_width*contour_frac, eff_height*contour_frac])
    py.contour(x, y, L.transpose(), conflevels, **contourargs)
    py.minorticks_on()
    py.xlabel(xlabel)
    py.ylabel(ylabel)
    axlim = py.axis()
    for limits in yconf:
        py.plot(axlim[0:2], [limits[0]]*2, **limitargs)
        py.plot(axlim[0:2], [limits[1]]*2, **limitargs)
    for limits in xconf:
        py.plot([limits[0]]*2, axlim[2:4], **limitargs)
        py.plot([limits[1]]*2, axlim[2:4], **limitargs)


    hax = py.subplot(2, 2, 1, position=[lbuffer, bbuffer + eff_height*contour_frac + 1e-5, eff_width*contour_frac, eff_height*(1. - contour_frac)])
    hax.plot(x, xppd, **posteriorargs)
    x_ylim = py.ylim()
    for limits in xconf:
        py.plot([limits[0]]*2, x_ylim, **limitargs)
        py.plot([limits[1]]*2, x_ylim, **limitargs)
    py.ylim(x_ylim)
    py.xlim(axlim[0:2])
    yt = py.yticks()[0]
    py.xticks(py.xticks()[0], [])
    py.yticks(yt[yt>0], [])
    py.ylim(x_ylim)
    py.xlim(axlim[0:2])
    py.minorticks_on()
    
    hay = py.subplot(2, 2, 4, position=[lbuffer + eff_width*contour_frac + 1e-5, bbuffer, eff_width*(1. - contour_frac), eff_height*contour_frac])
    hay.plot(yppd, y, **posteriorargs)
    y_xlim = py.xlim()
    for limits in yconf:
        py.plot(y_xlim, [limits[0]]*2, **limitargs)
        py.plot(y_xlim, [limits[1]]*2, **limitargs)
    py.ylim(axlim[2::])
    py.xlim(y_xlim)
    xt = py.xticks()[0]
    py.xticks(xt[xt>0], [])
    py.yticks(py.yticks()[0], [])
    py.xlim(y_xlim)
    py.ylim(axlim[2::])
    py.minorticks_on()    
    
    return cax, hax, hay
    

def areaOverlappingCircles(r1, r2, d):
    """Return overlapping area of two circles.  From Wolfram Mathworld.

    r1, r2 are the radii of the circles.

    d is the distance between their centers.
    """
    # 2012-07-30 15:40 IJMC: Created, from Wolfram Mathworld

    if r1> r2:
        temp = r1
        r1 = r2
        r2 = temp

    r1s = r1*r1
    r2s = r2*r2

    if not hasattr(d, '__iter__'):
        d = np.array([d])
    else:
        d = np.array(d, copy=False)

    valid_index = (d < (r1+r2)) * (d > (r2-r1))
    fulloverlap_index = d <= (r2-r1)

    ds = d[valid_index]**2


    area = np.zeros(d.shape)
    area[fulloverlap_index] = np.pi*r1s
    area[valid_index] = \
        r1s * np.arccos((ds + r1s - r2s) / (2*d[valid_index]*r1)) + \
        r2s * np.arccos((ds + r2s - r1s) / (2*d[valid_index]*r2)) - \
        0.5 * np.sqrt((-d[valid_index]+r1+r2)*(d[valid_index]+r1-r2)*(d[valid_index]-r1+r2)*(d[valid_index]+r1+r2))

    return area


def findRectangles(a, minsepy=None, minsepx=None, edgepad=10):
    """Find corner coordinates of approximate rectangle shapes in an array.

    :INPUTS:
      a : 2D Numpy array
        
      minsep : scalar
        Minimum separation from one upper or left-hand border to the
        next. (cf. :func:`find_peaks`)

      edgepad : int
        Pad the array with this many zeros on all sides (to find
        rectangles abutting array edges)
      """
    # 2012-08-09 23:51 IJMC: Created
    # 2012-08-17 14:53 IJMC: Added 'edgepad' option.
    # 2012-10-21 22:46 IJMC: Fixed 'x_borders' when edgepad is called.

    from scipy.signal import medfilt2d
    from analysis import pad

    dtype = a.dtype
    if dtype not in (np.float32, np.float64):
        a = a.astype(np.float32)

    #sm = medfilt2d(a, 5)
    sm = pad(medfilt2d(a, 5), a.shape[0]+edgepad*2, a.shape[1]+edgepad*2)
    smlox = np.hstack((np.diff(sm, axis=1), np.zeros((sm.shape[0],1))))
    smloy = np.vstack((np.diff(sm, axis=0), np.zeros(sm.shape[0])))

    boxcar = np.ones(5)/5.
    smloys = np.convolve(smloy.sum(1), boxcar, 'same')
    smloxs = np.convolve(smlox.sum(0), boxcar, 'same')

    y_lowers = np.sort(find_peaks(smloys, minsepy, thresh=10)[1]) - edgepad
    y_uppers = np.sort(find_peaks(-smloys, minsepy, thresh=10)[1]) - edgepad
    nreg = y_uppers.size
    corners = np.zeros((nreg, 4), dtype=int)
    corners[:,2] = y_lowers + 1
    corners[:,3] = y_uppers + 1
    for ii in range(nreg):
         subreg = smlox[y_lowers[ii]:y_uppers[ii]]
         xedges = np.abs(np.convolve(subreg.sum(0), boxcar, 'same'))
         x_borders = np.sort(find_peaks(xedges, minsepx)[1][0:2])
         corners[ii,0:2] = x_borders - edgepad+1

    return corners



def extractSubregion(fitsfilename, corners=None, dx=None, dy=None, kw='subreg', retall=False):
    """Extract a specified rectangular subregion from a FITS file.
    
    :INPUTS:
      fitsfilename : str
        Name of the (2D) FITS file.

      corners : str, 4-sequence
        if sequence: [x0, x1, y0, y1], corners of subregion.

        if str: header keyword containing this sequence. 

        In either case, the extracted subregion (when dx=dy=0) will be:
          data[corners[2]:corners[3], corners[0]:corners[1]]

      dx : None, 2-sequence
        If sequence: [x0, x1] will become [x0-dx[0], x1+dx[1]]

      dy : None, 2-sequence
        If sequence: [y0, y1] will become [y0-dy[0], y1+dy[1]]

      kw : None, str
        If str: this header keyword will be updated with the new
        corners (possibly modified by dx, dy)

    :OUTPUTS:
      (subregion_data, [fitsfile_header, corners_used])

      If the specified header keyword is not found, or the specified
      corners return an error, then this function will crash inelegantly.

    :NOTE:
      WCS headers will not be updated, so be careful when using this
      routine for imaging data!
    """
    # 2012-08-28 15:25 IJMC: Created
    # 2013-01-20 14:50 IJMC: Clarified documentation of 'corners' input.
    import pyfits

    if dx is None:
        dx = 0
    if dy is None:
        dy = 0
    if not hasattr(dx, '__iter__'):
        dx = [dx]
    if not hasattr(dy, '__iter__'):
        dy = [dy]
    if len(dx) < 2:
        dx = [dx[0], dx[0]]
    if len(dy) < 2:
        dy = [dy[0], dy[0]]



    if not isinstance(fitsfilename, np.ndarray):
        data = pyfits.getdata(fitsfilename)
        header = pyfits.getheader(fitsfilename)
        if isinstance(corners, str):
            corners = scrapeints(header[corners])
    else:
        data = fitsfilename
        header = pyfits.Header()

    ny, nx = data.shape


    newcorners = [max(0, corners[0]-dx[0]), min(nx, corners[1]+dx[1]), max(0, corners[2]-dy[0]), min(ny, corners[3]+dy[1])]
    subreg = data[newcorners[2]:newcorners[3], newcorners[0]:newcorners[1]]
    #pdb.set_trace()

    header.update('dx', str(dx))
    header.update('dy', str(dy))
    if kw is not None:
        header.update(kw, str(newcorners))

    if retall:
        ret = subreg, header, newcorners
    else:
        ret = subreg

    return ret


def scrapeints(string):
    """Extract a series of integers from a string.  Slow but robust.

    readints('[124, 56|abcdsfad4589.2]')

    will return:

    [124, 56, 4589, 2]
    """
    # 2012-08-28 16:19 IJMC: Created

    numbers = []
    nchar = len(string)
    thisnumber = ''
    for n, char in enumerate(string):
        try:
            val = int(char)
            anumber = True
        except:
            anumber = False

        if anumber:
            thisnumber = thisnumber + char
        else: #if (not anumber) or n==(nchar-1):
            if len(thisnumber)>0:
                numbers.append(int(thisnumber))
                thisnumber = ''
    return numbers


def array_or_filename(input, kw_getdata=dict(), kw_array=dict(), noneoutput=None):
    """If input is a Numpy array, return it.  If it is of type str,
    use Pyfits to read in the file and return it.  Keyword options are
    available for the calls to pyfits.getdata and numpy.array.  If
    input is None, return noneoutput."""
    # 2012-09-03 11:43 IJMC: Created
    from pyfits import getdata

    if input is None:
        output = noneoutput
    else:
        if isinstance(input, str):
            output = getdata(input, **kw_getdata)
        else:
            output = np.array(input, **kw_array)

    return output


def roundparams(value, error, seconderror=None, extraplaces=1):
    """Round a value and its associated error(s) to 2 decimal places.

    :INPUTS:
      value, error : scalars

    :OUTPUTS:
      value, error

    :EXAMPLE:
      ::

        import tools
        tools.roundparams(1.23456, .09876)

    :SEE_ALSO:
      :func:`roundvals`
    """

    # 2012-09-27 17:47 IJMC: Created
    if seconderror is not None:
        pow = np.abs(np.floor(np.log10(0.5*(np.abs(error)+np.abs(seconderror)))))
    else:
        pow = np.abs(np.floor(np.log10(np.abs(error))))
    value = np.round(value, decimals=pow+extraplaces)
    val1 = np.round(error, decimals=pow+extraplaces)
    ret = value, val1
    fstrp = '%+'+ ('1.%i' % (pow+1)) + 'f'
    fstrm = '%+'+ ('1.%i' % (pow+1)) + 'f'
    if seconderror is not None:
        val2 = np.round(seconderror, decimals=pow+1)
        ret = ret + (val2,)
        fstr =  '$%s^{%s}_{%s}$' % (fstrp,fstrp,fstrm)
        retstr = fstr % (value, val1, val2)
    else:
        fstr =  '$%s \\pm %s$' % (fstrp,fstrp)
        retstr = fstr % (value, val1)

    return retstr


def sample_2dcdf(cdf, x, y, nsamp=1, verbose=False):
    """Sample a 2D Probability Distribution Function (2d-histogram)

    :INPUTS:
      CDF : 2D NumPy array
        Two-dimensional (N x M) probability distribution function
        (histogram) from which you wish to draw samples. This need not
        be in any particular normalized form -- the only condition is
        that the value in each cell is proportional to the probability
        of that cell.

      x : 1D NumPy array
        N+1 Values indexing the cells of the CDF (one per row), evenly spaced

      y : 1D NumPy array
        M+1 Values indexing the cells of the CDF (one per column), evenly spaced

      nsamp : int
        Number of samples to be drawn.

    :NOTES:
      Note that this approach uses simple, dumb nearest-neighbor
      interpolation when drawing candidate samples.  Exceedingly
      granular CDFs may suffer.  

      Note that this code isn't optimized; try it with small samples
      before you use it on large ones!

    """
    # 2012-09-30 21:44 IJMC: Created

    # Initialize:
    cdf = np.array(cdf, copy=False)
    x = np.array(x, copy=False)
    y = np.array(y, copy=False)
    ret1 = np.zeros(nsamp, dtype=y.dtype)
    ret2 = np.zeros(nsamp, dtype=x.dtype)
    xcenters = 0.5 * (x[1:] + x[0:-1])
    ycenters = 0.5 * (y[1:] + y[0:-1])


    #NOrmalize:
    prob = 1.0 * cdf / cdf.max()

    success_estimate = 1.0 * prob.sum() / prob.size
    n_maxsamp = int(2*nsamp/success_estimate)
    if verbose:
        print "Success estimate is %1.3e" % success_estimate
        print "Will need to draw approximately %i samples to generate desired number (%i)." % (n_maxsamp, nsamp)

    # Generate initial random samples to draw from:
    sampx = np.random.uniform(low=x.min(), high=x.max(), size=n_maxsamp)
    sampy = np.random.uniform(low=y.min(), high=y.max(), size=n_maxsamp)
    accept_prob = np.random.uniform(size=n_maxsamp)

    nsampled = 0
    ii = 0
    while nsampled < nsamp:
        # Locate the candidate sample in the CDF array: 
        ind_x = (np.abs(xcenters-sampx[ii])==np.abs(xcenters-sampx[ii]).min())#.nonzero()
        ind_y = (np.abs(ycenters-sampy[ii])==np.abs(ycenters-sampy[ii]).min())#.nonzero()

        # Accept or reject the candidate sample:
        if prob[ind_x, ind_y] > accept_prob[ii]:
            ret1[nsampled] = sampx[ii]
            ret2[nsampled] = sampy[ii]
            nsampled += 1
        ii += 1

        # If we ran out of random samples, draw more:
        if ii >= n_maxsamp:
            sampx = np.random.uniform(low=x.min(), high=x.max(), size=n_maxsamp)
            sampy = np.random.uniform(low=y.min(), high=y.max(), size=n_maxsamp)
            accept_prob = np.random.uniform(size=n_maxsamp)
            ii = 0
            
    return ret1, ret2


def sample_1dcdf(pdf, x, nsamp=1, verbose=False, absnorm=False):
    """Sample a 1D Posterior Distribution Function (1d-histogram)

    :INPUTS:
      PDF : 1D NumPy array
        Distribution function (histogram) from which you wish to draw
        samples. This need not be in any particular normalized form --
        the only condition is that the value in each cell is
        proportional to the probability of that cell.

      x : 1D NumPy array
        N Values indexing the cells of the CDF (one per row)

      nsamp : int
        Number of samples to be drawn.

      absnorm : bool
        If True, normalize pdf so that it integrates to 1.0

    """
    # 2012-10-28 13:13 IJMC: Created.
    # 2012-11-20 16:18 IJMC: Fixed a small bug.

    # Initialize:
    argind = np.argsort(x)
    x = np.array(x, copy=False)[argind]
    pdf = np.array(pdf, copy=False)[argind]

    ret = np.zeros(nsamp, dtype=x.dtype)

    #NOrmalize:
    #if absnorm:
    prob = 1.0 * pdf / pdf.max()
    #else:
    #    prob = pdf

    success_estimate = 1.0 * prob.sum() / prob.size
    n_maxsamp = int(2*nsamp/success_estimate)
    if verbose:
        print "Success estimate is %1.3e" % success_estimate
        print "Will need to draw approximately %i samples to generate desired number (%i)." % (n_maxsamp, nsamp)

    # Generate initial random samples to draw from:

    nsampled = 0
    ii = 0
    
    while nsampled < nsamp:
        sampx = np.random.uniform(low=x.min(), high=x.max(), size=n_maxsamp)
        accept_prob = np.interp(sampx, x, prob, left=0.0, right=0.0)
        accept_ind = np.random.uniform(low=0., high=1., size=n_maxsamp) < accept_prob

        n_accepted = accept_ind.sum()
        ind0 = nsampled
        ind1 = min(nsampled+n_accepted, nsamp)
        ret[ind0:ind1] = sampx[accept_ind][0:ind1-ind0]
        nsampled += n_accepted
        #pdb.set_trace()
        #if verbose:
        #    print "%i sampled accepted" % n_accepted
        #    print nsampled, n_accepted, nsamp
            
    return ret


def textfig(textlist, **kw):
    """Generate a figure from a list of strings.

    :INPUTS:
      textlist : list
         List of text strings, one per entry.

    :OPTIONS:
      any options valid for passing to matplotlib.text
      
    :RETURNS:
      (fig, ax) -- handles to Figure and Axes that were created
      """
    # 2013-03-10 13:27 IJMC: Created
    import pylab as py

    # Handle input options:
    defaults = dict(horizontalalignment='left', fontsize=9, family='Courier', weight='bold')
    if kw.has_key('figoptions'):
        figoptions = kw.pop('figoptions')
    else:
        figoptions = dict()

    if kw.has_key('axoptions'):
        axoptions = kw.pop('axoptions')
    else:
        axoptions = dict(position=[.02, .02, .96, .9])

    for key in defaults.keys():
        if not kw.has_key(key):
            kw[key] = defaults[key]

    # Generate figure and axis:
    fig = py.figure(nextfig(), **figoptions)
    ax = py.axes(**axoptions)

    nlines = len(textlist)
    vertpos = py.linspace(.95, .05, nlines)
    py.xticks([])
    py.yticks([])

    # Plot the text:
    for line, pos in zip(textlist, vertpos):
        py.text(.05, pos, line, transform=ax.transAxes, **kw)

    return fig, ax

def roundvals(input, ndigit=2, strout=True):
    """Round all input values to the smallest number of digits used.

     :INPUTS:
      input : scalar, or 1D list or array
        values to be rounded

      ndigit : int
        Number of significant digits to be retained

      strout : bool
        Return text strings?  If not, return a 1D NumPy array.
        
     :EXAMPLE:
       ::
       
         import tools
         tools.roundvals([235, -4, 0.045380])

     :SEE_ALSO:
       :func:`roundparams`

        """
    # 2013-03-11 08:29 IJMC: Created
    # 2013-04-20 09:35 IJMC: Now, don't be thrown off by zero values.

    if not hasattr(input, '__iter__'):
        input = [input]
        scalarInput = True
    else:
        scalarInput = False

    input = np.array(input, copy=False)
    #np.abs(input).min()
    ndig = np.abs(np.floor(np.log10(np.abs(input[np.abs(input)>0])).min()).astype(int)) + (ndigit-1)
    if strout:
        precstr = '%' + ('1.%i' % ndig) + 'f'
        ret = [precstr % val for val in input]
    else:
        ret = np.round(input, ndig)
    
    if scalarInput:
        ret = ret[0]

    return ret


def gelman_rubin(chains):
    """Compute the Gelman-Rubin convergence metric for MCMC chains.

    :INPUTS:
      chains : 2D NumPy array
        A stack of all MCMC chains to be compared, created with
        something like :func:`numpy.vstack`.  The chains must all be
        the same length, and they must have more links than the total
        number of chains.

        OR

      chains : 3D NumPy array
        N chains of L links for P parameters (e.g., the 'chain'
        attribute of an emcee.sampler object), of shape NxLxP.

    :OUTPUTS:
      R metric. If this is 'close to 1.0,' then your chains have
      converged to the same distribution.  The definition of 'close'
      could be 1.2, 1.1, 1.01... it's up to you!

    :REFERENCE:
      Eq. 20 of Gelman & Rubin 1992, Statistical Sciences, Vol. 7, p. 457

      http://www.star.le.ac.uk/~sav2/idl/rhat.pro
      """
    # 2014-01-23 11:14 IJMC: Created by copying from the IDL code
    #                        located at
    #                        http://www.star.le.ac.uk/~sav2/idl/rhat.pro


    #---------------------------------------------------------
    #Main routine

    chains = np.array(chains, copy=False)
    if 0 in chains.shape:
        print "Input array cannot have any dimensions of size zero! Returning 9e99."
        return np.array([9e99])

    if chains.ndim==2:
        if (chains.shape[0] > chains.shape[1]):
            dimension = 0
        else:
            dimension = 1
        nn = chains.shape[dimension]

        #mean value of each chain
        mean_j = chains.mean(axis=dimension)

        #variance within each chain
        var_j = chains.var(axis=dimension)

        #now compute B and W from Gelman & Rubin
        B = nn*mean_j.var()
        W = var_j.mean()

        #compute Gelman-Rubin R^2
        R2 = ( W*(nn-1)/nn + B/nn ) / W
        R  = np.array([np.sqrt(R2)])

        if verbose:
            print '-- Mean of each chain: j=0,1,2,...'
            print mean_j
            print '-- Variance of each chain: j=0,1,2,...'
            print var_j
            print '-- B/N, W, R^2, R, N'
            print (B/nn),W,R2,R,nn

    elif chains.ndim==3:
        nn = chains.shape[1]
        mean_j = chains.mean(axis=1)
        var_j = chains.var(axis=1)
        B = nn * mean_j.var(axis=0)
        W = var_j.mean(axis=0)
        R2 = ( W*(nn-1)/nn + B/nn ) / W
        R  = np.sqrt(R2)

    else:
        print "Must input a 2D or 3D array! Returning 9e99."
        R = np.array([9e99])

    return R


def get_emcee_start(bestparams, variations, nwalkers, maxchisq, args, homein=True, retchisq=False, depth=np.inf):
    """Get starting positions for EmCee walkers.

    :INPUTS:
      bestparams : sequence (1D NumPy array)
        Optimal parameters for your fitting function (length N)

      variations : 1D or 2D NumPy array
        If 1D, this should be length N and new trial positions will be
        generated using numpy.random.normal(bestparams,
        variations). Thus all values should be greater than zero!

        If 2D, this should be size (N x N) and we treat it like a
        covariance matrix; new trial positions will be generated using
        numpy.random.multivariate_normal(bestparams, variations). 

      nwalkers : int
        Number of positions to be chosen.

      maxchisq : int
        Maximum "chi-squared" value for a test position to be
        accepted.  In fact, these values are computed with
        :func:`phasecurves.errfunc` as errfunc(test_position, *args)
        and so various priors, penalty factors, etc. can also be
        passed in as keywords.

      args : tuple
        Arguments to be passed to :func:`phasecurves.errfunc` for
        computing 'chi-squared' values.

      homein : bool
        If True, "home-in" on improved fitting solutions. In the
        unlikely event that a randomly computed test position returns
        a better chi-squared than your specified best parameters,
        reset the calculation to start from the new, improved set of
        parameters.

      retchisq : bool
        If True, return the tuple (positions, chisq_at_positions)
        
     :BAD_EXAMPLE:
      ::

        pos0 = tools.get_emcee_start(whitelight_bestfit[0], np.abs(whitelight_bestfit[0])/1000., nwalkers, 10*nobs, mcargs)
        """
    # 2013-05-01 11:18 IJMC: Created
    
    #get_emcee_start(bestparams, variations, nwalkers, maxchisq, args):
    from phasecurves import errfunc

    best_chisq = errfunc(bestparams, *args)
    if best_chisq >= maxchisq:
        print "Specified parameter 'maxchisq' is larger than the chi-squared value for the specified best parameters. Try increasing maxchisq."
        return -1

    npar = len(bestparams)
    if variations.ndim==2:
        usecovar = True
    else:
        usecovar = False

    pos0 = np.zeros((nwalkers, npar), dtype=float)
    chisq = np.zeros(nwalkers, dtype=float)
    npos = 0
    while npos < nwalkers:
        if usecovar:
            testpos = np.random.multivariate_normal(bestparams, variations)
        else:
            testpos = np.random.normal(bestparams, variations)
        testchi = errfunc(testpos, *args)
        if np.isfinite(testchi) and (testchi < best_chisq) and homein and depth>0:
            return get_emcee_start(testpos, variations, nwalkers, maxchisq, args, homein=homein, retchisq=retchisq, depth=depth-1)
        elif testchi < maxchisq:
            pos0[npos] = testpos
            chisq[npos] = testchi
            npos += 1

    if retchisq:
        ret = pos0, chisq
    else:
        ret = pos0
    return ret


def findFrac(validValues, thisValue, retinds=False):
  """
  Helper tool for simply linear interpolation.

   :INPUTS:
     validValues : sequence
       List of valid values

     thisValue : scalar
       Value of interest.

   :OUTPUTS:
     (TwoClosestValues, relativeFractions)
  """
  # 2013-08-12 09:48 IJMC: Created
  validValues = np.array(validValues, copy=False)
  if thisValue<=validValues.min():
    ret = ([validValues.min()], [1.0])
  elif thisValue>=validValues.max():
    ret = ([validValues.max()], [1.0])
  else:
    dif = np.abs(validValues - thisValue)
    inds = np.argsort(dif)[0:2]
    val1, val2 = validValues[inds]
    c1 = np.abs(1.0*(val2-thisValue)/(val2-val1))
    c2 = np.abs(1.0*(val1-thisValue)/(val2-val1))
    ret = ([val1, val2], [c1, c2])

  if retinds:
      inds = [np.nonzero(validValues==val)[0][0] for val in ret[0]]
      ret = (inds, ret[1])
  return ret



def feps_interpol(x, y, a, linear=True):
    """
    Wrapper script for NumPy interpolation. Culls duplicate values and
    puts x into a monotonically increasing grid.

    :INPUTS:
      x : NumPy array
        1D sequence of values defining the grid coordinates at which
        the input values 'y' are defined.

      y : NumPy array
        1D sequence of values.

      a : NumPy array
        Values of 'x' at which to interpolate from the values of 'y'.

    :EXAMPLE:
     ::

      import numpy as np
      
      x = np.linspace(-3,3,61)
      v = np.sin(x)
      u = np.array([-2.50, -2.25, -1.85, -1.55, -1.20, -0.85, -0.50, -0.10, \
          0, 0.75, 0.85, 1.05, 1.45, 1.85, 2.00, 2.25, 2.75 ])
      b = feps_interpol(x,v,u)
      

    :NOTES:
      Converted from IDL code from J. Bouwman. Documentation was:
      ;(SH Feb 26 1999)
      ;We need to make the grid mononic therefore spline needs to do
      ;some cleaning before execution
    """
    # 2013-12-01 23:47 IJMC: Translated from IDL.

    u, idx = np.unique(x, return_index=True)
    xt = x[idx]
    yt = y[idx]
    idx2 = np.argsort(xt)
    xt = xt[idx2]
    yt = yt[idx2]
    idx3 = np.argsort(a)
    at = a[idx3]
    if linear is not True:
        print "nonlinear mode not implemented!"
    else:
        bt = np.interp(at, xt, yt)
    return bt


def removeoutliers(data, nsigma, remove='both', center='mean', niter=99999, retind=False, verbose=False):
    """Strip outliers from a dataset, iterating until converged.

    :INPUT:
      data -- 1D numpy array.  data from which to remove outliers.

      nsigma -- positive number.  limit defining outliers: number of
                standard deviations from center of data.

    :OPTIONAL INPUTS:               
      remove -- ('min'|'max'|'both') respectively removes outliers
                 below, above, or on both sides of the limits set by
                 nsigma.

      center -- ('mean'|'median'|value) -- set central value, or
                 method to compute it.

      niter -- number of iterations before exit; defaults to Inf,
               which can occasionally result in empty arrays returned

      retind -- (bool) whether to return index of good values as
                second part of a 2-tuple.

    :EXAMPLE: 
       ::

           from numpy import hist, linspace, randn
           from analysis import removeoutliers
           data = randn(1000)
           hbins = linspace(-5,5,50)
           d2 = removeoutliers(data, 1.5, niter=1)
           hist(data, hbins)
           hist(d2, hbins)

       """
    # 2009-09-04 13:24 IJC: Created
    # 2009-09-24 17:34 IJC: Added 'retind' feature.  Tricky, but nice!
    # 2009-10-01 10:40 IJC: Added check for stdev==0
    # 2009-12-08 15:42 IJC: Added check for isfinite

    from numpy import median, ones, isfinite

    def getcen(data, method):
        "Get central value of a 1D array (helper function)"
        if method.__class__==str:
            if method=='median':
                cen = median(data)
            else:
                cen = data.mean()
        else:
            cen = method
        return cen

    def getgoodindex(data, nsigma, center, stdev, remove):
        "Get number of outliers (helper function!)"
        if stdev==0:
            distance = data*0.0
        else:
            distance = (data-center)/stdev
        if remove=='min':
            goodind = distance>-nsigma
        elif remove=='max':
            goodind = distance<nsigma
        else:
            goodind = abs(distance)<=nsigma
        return goodind

    data = data.ravel().copy()

    ndat0 = len(data)
    ndat = len(data)
    iter=0
    goodind = ones(data.shape,bool)
    goodind *= isfinite(data)
    while ((ndat0<>ndat) or (iter==0)) and (iter<niter) and (ndat>0) :
        ndat0 = len(data[goodind])
        cen = getcen(data[goodind], center)
        stdev = data[goodind].std()
        thisgoodind = getgoodindex(data[goodind], nsigma, cen, stdev, remove)
        goodind[find(goodind)] = thisgoodind
        if verbose:
            print "cen>>",cen
            print "std>>",stdev
        ndat = len(data[goodind])
        iter +=1
        if verbose:
            print ndat0, ndat
    if retind:
        ret = data[goodind], goodind
    else:
        ret = data[goodind]
    return ret

def egaussian(p, x, y, e=None):
    """ Compute the deviation between the values y and the gaussian defined by p, x:

    p is a three- or four-component array, list, or tuple.

    Returns:   y - p3 - p0/(p1*sqrt(2pi)) * exp(-(x-p2)**2 / (2*p1**2))

    if an error array, e (typ. one-sigma) is entered, the returned value is divided by e.

    SEE ALSO:  :func:`gaussian`"""
    # 2008-09-11 15:19 IJC: Created
    # 2009-09-02 15:20 IJC: Added weighted case
    # 2011-05-18 11:46 IJMC: Moved to analysis.
    from numpy import ones

    if e==None:
        e=ones(x.shape)
    fixval(e,y.max()*1e10)

    z = (y - gaussian(p, x))/e
    fixval(z,0)

    return z

def gaussian(p, x):
    """ Compute a gaussian distribution at the points x.

        p is a three- or four-component array, list, or tuple:

        y =  [p3 +] p0/(p1*sqrt(2pi)) * exp(-(x-p2)**2 / (2*p1**2))

        p[0] -- Area of the gaussian
        p[1] -- one-sigma dispersion
        p[2] -- central offset (mean location)
        p[3] -- optional constant, vertical offset

        NOTE: FWHM = 2*sqrt(2*ln(2)) * p1  ~ 2.3548*p1

        SEE ALSO:  :func:`egaussian`"""
    #2008-09-11 15:11 IJC: Created for LINEPROFILE
    # 2011-05-18 11:46 IJC: Moved to analysis.
    # 2013-04-11 12:03 IJMC: Tried to speed things up slightly via copy=False
    # 2013-05-06 21:42 IJMC: Tried to speed things up a little more.

    if not isinstance(x, np.ndarray):
        x = array(x, dtype=float, copy=False)

    if len(p)==3:
        p = array(p, copy=True)
        p = concatenate((p, [0]))
    #elif len(p)==4:
    #    p = array(p, copy=False)

    return  p[3] + p[0]/(p[1]*sqrt(2*pi)) * exp(-(x-p[2])**2 / (2*p[1]**2))

def stdr(x, nsigma=3, niter=99999, finite=True, verbose=False, axis=None):
    """Return the standard deviation of an array after removing outliers.
    
    :INPUTS:
      x -- (array) data set to find std of

    :OPTIONAL INPUT:
      nsigma -- (float) number of standard deviations for clipping
      niter -- number of iterations.
      finite -- if True, remove all non-finite elements (e.g. Inf, NaN)
      axis -- (int) axis along which to compute the mean.

    :EXAMPLE:
      ::

          from numpy import *
          from analysis import stdr
          x = concatenate((randn(200),[1000]))
          print std(x), stdr(x, nsigma=3)
          x = concatenate((x,[nan,inf]))
          print std(x), stdr(x, nsigma=3)

    SEE ALSO: :func:`meanr`, :func:`medianr`, :func:`removeoutliers`, 
              :func:`numpy.isfinite`
    """
    # 2010-02-16 14:57 IJC: Created from mear
    # 2010-07-01 14:06 IJC: ADded support for higher dimensions
    from numpy import array, isfinite, zeros, swapaxes

    x = array(x, copy=True)
    xshape = x.shape
    ndim =  len(xshape)
    if ndim==0:
        return x

    if  ndim==1 or axis is None:
        # "1D" array
        x = x.ravel()
        if finite:
            x = x[isfinite(x)]
        x = removeoutliers(x, nsigma, niter=Inf, verbose=verbose)
        return x.std()
    else:
        newshape = list(xshape)
        oldDimension = newshape.pop(axis) 
        ret = zeros(newshape, float)

        # Prevent us from taking the action along the axis of primary incidices:
        if axis==0: 
            x=swapaxes(x,0,1)

        if axis>1:
            nextaxis = axis-1
        else:
            nextaxis = 0

        for ii in range(newshape[0]):
            #print 'x.shape>>',x.shape, 'newshape>>',newshape, 'x[ii].shape>>',x[ii].shape, 'ret[ii].shape>>',ret[ii].shape,'ii>>',ii
            ret[ii] = stdr(x[ii], nsigma=nsigma,niter=niter,finite=finite,\
                                 verbose=verbose, axis=nextaxis)
        return ret

