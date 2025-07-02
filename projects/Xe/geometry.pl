# Geometries can be specified either manually in the array 'geometries'
# by copying what is between start copy and end copy for each geometry
# or automatically as shown below the hash array %geometry
#
# Note that 'length_unit' should be set to the actual length unit used for values in 'atoms' array,
# later (in &generate_geometries) $model{'r_unit'} is set also to this unit for consistency reasons
#
# Note that there's a restriction on the orientation of the molecules in the scripts (although there isn't one in the codes themselves)
# Your molecular target should be oriented as follows for the corresponding point groups ('any orientation' still requires the axes to
# be along the cartesian axes):
# D2h: any orientation
# C2v: C2 axis along Z        
# C2h: C2 axis along Z        
# D2:  any orientation       
# C2:  C2 axis along Z 
# Cs:  Molecule in the YZ plane
# Ci:  any orientation
# C1:  any orientation 

%geometry = (

  'suffix', ".Xe_atom",   # string added to model directory to distinguish runs
                              # with different geometry settings

  'geometry_labels',  "   Xe atom   ",         # labels used on the first line of output files
                                          # it should correspond to numbers given in 'geometries'->'description'
  'correct_cm',   1,                      # correct the center of mass to be at the origin
  'length_unit',  0,                      # 0 - atomic units, 1 - Angstroms   

  'geometries', [
       # start copy
#      { 'description', "   1.81     1.81     104.48", # string to use in output files to describe this particular geometry, can be anything
#        'gnuplot_desc', "R1 = 1.81, R2 = 1.81, Theta = 104.48",    # used in gnuplot files for keys
#        # specify ALL atoms (even redundant with respect to symmetry elements) 
#        'atoms', [ [ "O",           0.00,          0.00,          0.00 ],
#                   [ "H",           0.00,   1.430954704,  -1.108363043 ],
#                   [ "H",           0.00,  -1.430954704,  -1.108363043 ] ]
#      },
#      # end copy
  ],

  # the following option can be used e.g. to continue an interupted run
  # all geometries are generated but codes will run only for specified geometries
  'start_at_geometry', 1,                # codes will run only for geometries with index >= than this number
  'stop_at_geometry',  0,               # this can be used to stop at certain geometry 
                                         # if zero then codes will run for all geometries
);

# Here is an example of automatic generation of geometries.
# Specifically the symmetric stretch of the water molecule is generated.
# If you want to use this then comment first the part in %geometry between # start copy ... # end copy
# and change 'suffix' and 'geometry_labels' in %geometry accordingly, 
# e.g. as ".sym_stretch.theta104.48.r1.6-0.2-2.6" and "     O-H     "
#
# Xenon atom geometry - single atom at origin
push(@{$geometry{'geometries'}},
     { 'description', "   Xe atom  ",
       'gnuplot_desc', "Xe atom",
       'atoms', [ [ "Xe", 0.0, 0.0, 0.0 ] ]
     }
);