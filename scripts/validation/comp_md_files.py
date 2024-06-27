#!/usr/bin/env python

import comp_md

#print comp_md.validate_md_time("/Users/snjoroge/Documents/openmd/development/samples/alkane/butane.omd", "/Users/snjoroge/Documents/openmd/development/samples/alkane/butane.omd")
#print comp_md.file_is_md("/home/sasa/Documents/development/samples/alkane/butane.omd")
#print comp_md.compare("../samples/argon/500.stat", "../samples/argon/500_v.stat")
comp_md.setupDirectories()
comp_md.scanForMdFiles()
comp_md.runMdFiles()
comp_md.cleanUp()
#print comp_md.compareEor("../samples/argon/500.eor", "../samples/argon/500_v.eor")
