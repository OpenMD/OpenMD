These are sample files used to test the bond order parameter
code in StaticProps.  Perfect BCC, FCC, HCP, and Icosahedral
clusters have known analytic values for the bond order parameter
and the files in this directory can be used to make sure the
correct values are being computed.  The central atom in each
cluster is a copper atom, so the proper way to run the tests would
be a command like:

  $(OOPSE_HOME)/bin/StaticProps --bo -i bcc.md --rcut=9 --sele1="select Cu"
