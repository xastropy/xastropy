def get_package_data():
    # Installs the testing data files. Unable to get package_data
    # to deal with a directory hierarchy of files, so just explicitly list.
    return {'xastropy.igm.abs_sys.tests': ['files/*.fits', 'files/*.dat',
                                   'files/*.all']}
