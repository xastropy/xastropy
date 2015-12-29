""" Module to build DLA summary files and spectra for pyigm
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import glob, os, sys, io, json

from astropy import units as u
from astropy.table import Table, Column, QTable

from linetools.spectra import io as lsio
from linetools.abund import ions as ltai

from pyigm.surveys.dlasurvey import DLASurvey

from xastropy.xutils import fits as xxf
from xastropy.xutils import xdebug as xdb

def neeleman13():
    """ Build a summary file for the Neeleman+13 sample
    """
    prefix = 'H100'
    outpath = os.getenv('DROPBOX_DIR')+'/Public/DLA/'+prefix+'/'
    dlasurvey = DLASurvey.from_flist('Lists/Neeleman13.lst',
                                     tree=os.environ.get('DLA'))
    dlasurvey.ref = 'Neeleman+13'
    #
    dlasurvey.fill_ions(use_Nfile=True)
    mk_json_ions(dlasurvey, prefix, outpath+prefix+'_DLA_ions.json')

    # Summary file and spectra
    mk_summary(dlasurvey, prefix, outpath+prefix+'_DLA.fits',
               specpath=outpath+'/Spectra/',
               htmlfil=outpath+prefix+'_DLA.html')

##
def fits_idx(idx):
    """Converts an index to instrument
    """
    fits_list = dict(zip(list(map((lambda x: 2**x),range(6))),
                         ['HIRES','ESI','UVES','XX','MIKEb','MIKEr']))
    try:
        return fits_list[idx]
    except:
        return 'Unknown'

##
def survey_name(prefix, idla):
    """
    Parameters
    ----------
    prefix : str
      Survey prefix
    idla : DLASystem

    Returns
    -------
    jname : str
      Name in J2000

    """
    jname = prefix+'_J{:s}{:s}_z{:0.3f}'.format(
        idla.coord.ra.to_string(unit=u.hour,pad=True,sep='',precision=2),
        idla.coord.dec.to_string(pad=True,alwayssign=True,sep='',precision=1),
        idla.zabs)
    return jname


####
def mk_1dspec(idla, outpath=None, name=None, clobber=False):
    """ Collate and rename the spectra
    Parameters
    ----------
    idla : DLASystem
    name : str, optional
    clobber : bool, optional

    """
    #
    if outpath is None:
        outpath = 'Spectra/'
    if name is None:
        raise ValueError('Not setup for this')

    # Spectra files
    spec_dict = idla._clmdict['fits_files']
    all_spec = []
    for key in spec_dict.keys():
        instr= fits_idx(key)
        if instr == 'XX':
            xdb.set_trace()
            if 'PROGETTI' in spec_dict[key]:
                instr = 'MAGE'
            elif 'HIRES' in spec_dict[key]:
                instr = 'HIRES'
            else:
                xdb.set_trace()
        # Generate new filename
        spec_fil = name+'_'+instr+'.fits'
        # Copy over?
        tmp = glob.glob(outpath+spec_fil+"*")
        if (len(tmp) == 0) | (clobber):
            # Read
            spec = lsio.readspec(spec_dict[key])
            # Write
            spec.write_to_fits(outpath+spec_fil, clobber=True, add_wave=True)
        # Append
        all_spec.append(str(spec_fil))
    # Return
    return all_spec

####
def mk_summary(dlas, prefix, outfil, specpath=None, htmlfil=None):
    """ Loops through the DLA list and generates a Table

    Also pushes the 1D spectra into the folder

    Parameters
    ----------
    dlas : DLASurvey
    prefix : str
    outfil : str
      Name of the output FITS summary file
    htmlfil : str, optional

    Returns
    -------
    """
    #
    if htmlfil is None:
        htmlfil = 'tmp.html'

    # # Constructing
    # QSO, RA/DEC
    cqso = Column(dlas.qso, name='QSO')
    ra = dlas.coord.ra.degree[0]
    dec = dlas.coord.dec.degree[0]
    jname = []
    for abs_sys in dlas._abs_sys:
        jname.append(survey_name(prefix, abs_sys))

    cjname = Column(jname, name='Name')
    cra = Column(ra, name='RA', unit=u.degree)
    cdec = Column(dec, name='DEC', unit=u.degree)
    czem = Column(dlas.zem, name='Z_QSO')

    # Begin the Table
    dla_table = QTable( [cjname, cqso, cra, cdec, czem] )

    # LLS properties
    czabs = Column(dlas.zabs, name='ZABS')
    cNHI = Column(dlas.NHI, name='logNHI')
    csigNHI = Column(dlas.sig_NHI, name='sig(logNHI)')

    # Add to Table
    dla_table.add_columns([czabs, cNHI, csigNHI])

    # Spectra files
    all_sfiles = []
    for jj,ills in enumerate(dlas._abs_sys):
        sub_spec = mk_1dspec(ills, name=cjname[jj], outpath=specpath)
        # Pad
        while len(sub_spec) < 5:
            sub_spec.append(str('NULL'))
        # Append
        all_sfiles.append(sub_spec)

    cspec = Column(np.array(all_sfiles), name='SPEC_FILES')
    dla_table.add_column( cspec )

    # Sort
    dla_table.sort('RA')

    # Write
    print('Writing {:s}'.format(outfil))
    xxf.table_to_fits(dla_table,outfil)
    print('Writing {:s}'.format(htmlfil))
    Table(dla_table).write(htmlfil)

    return dla_table

#
def mk_json_ions(dlas, prefix, outfil):
    """ Generate a JSON table of the Ion database
    Parameters
    ----------
    dlas : DLASurvey
    prefix : str
    outfil : str
      Output JSON file
    """
    # Sort
    ra = dlas.coord.ra.degree[0]
    srt = np.argsort(np.array(ra))

    all_ions = {}
    # Loop on DLA
    for jj, isrt in enumerate(srt):
        idla = dlas._abs_sys[isrt]
        # Astropy Table
        ion_tab = idla._ionN
        # Convert key to standard names
        new_dict = {}
        for row in ion_tab:
            Zion = (row['Z'], row['ion'])
            # Skip HI
            if Zion == (1,1):
                continue
            # Get name
            new_key = ltai.ion_name(Zion)
            # Fine structure?
            if row['Ej'] > 0.:
                new_key = new_key+'*'
            new_dict[new_key] = dict(zip(row.dtype.names, row))
        # Write to all_ions
        name = survey_name(prefix, idla)
        all_ions[name] = new_dict

    # Write
    print('Writing {:s}'.format(outfil))
    with io.open(outfil, 'w', encoding='utf-8') as f:
        f.write(unicode(json.dumps(all_ions, sort_keys=True, indent=4,
            separators=(',', ': '))))

    # Return
    return all_ions


# ##################################################
if __name__ == '__main__':

    if len(sys.argv) == 1:

        flg_smpl = 0
        flg_smpl += 2**0  # Neeleman+13

        if (flg_smpl % 2**1) >= 2**0:
            neeleman13()
    else: # Run something
        pass