#!/usr/bin/env python

from subprocess import check_output, call
import os

ggoFiles = {
    'Dump2XYZ':            'src/applications/dump2Xyz/Dump2XYZ.ggo',
    'DynamicProps':        'src/applications/dynamicProps/DynamicProps.ggo',
    'Hydro':               'src/applications/hydrodynamics/Hydro.ggo',
    'icosahedralBuilder':  'src/applications/nanoparticleBuilder/icosahedralBuilder.ggo',
    'nanoparticleBuilder': 'src/applications/nanoparticleBuilder/nanoparticleBuilder.ggo',
    'nanorod_pentBuilder': 'src/applications/nanoparticleBuilder/nanorod_pentBuilder.ggo',
    'nanorodBuilder':      'src/applications/nanoparticleBuilder/nanorodBuilder.ggo',
    'omd2omd':             'src/applications/omd2omd/omd2omd.ggo',
    'randomBuilder':       'src/applications/randomBuilder/randomBuilder.ggo',
    'recenter':            'src/applications/recenter/recenter.ggo',
    'SequentialProps':     'src/applications/sequentialProps/SequentialProps.ggo',
    'simpleBuilder':       'src/applications/simpleBuilder/simpleBuilder.ggo',
    'StaticProps':         'src/applications/staticProps/StaticProps.ggo',
    'thermalizer':         'src/applications/thermalizer/thermalizer.ggo'
}
        
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


ggoProgram = which('gengetopt')
gitProgram = which('git')

topLevel = check_output([gitProgram, "rev-parse", "--show-toplevel"]).strip()


for ggo, ggoFile in ggoFiles.items():
    theDir = os.path.dirname(topLevel + os.sep + ggoFile)
    theFile = os.path.basename(topLevel + os.sep + ggoFile)
    os.chdir(theDir)

    print "Generating gengetopt parser from %s in %s" % (theFile, os.getcwd())

    call([ggoProgram, "-C", "--input=" + theFile])

    
