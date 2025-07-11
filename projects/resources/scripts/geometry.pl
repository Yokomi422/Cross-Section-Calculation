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

  'suffix', ".equilibrium",   # string added to model directory to distinguish runs
                              # with different geometry settings

  'geometry_labels',  "   R1       R2       Theta    ",         # labels used on the first line of output files
                                          # it should correspond to numbers given in 'geometries'->'description'
  'correct_cm',   1,                      # correct the center of mass to be at the origin
  'length_unit',  0,                      # 0 - atomic units, 1 - Angstroms   

  'geometries', [
       # start copy
       { 'description', "   1.81     1.81     104.48", # string to use in output files to describe this particular geometry, can be anything
         'gnuplot_desc', "R1 = 1.81, R2 = 1.81, Theta = 104.48",    # used in gnuplot files for keys
         # specify here all atoms, use Cartesian coordinates
         'atoms', [ [ "O",           0.00,          0.00,          0.00 ],
                    [ "H",           0.00,   1.430954704,  -1.108363043 ],
                    [ "H",           0.00,  -1.430954704,  -1.108363043 ] ]
       },
       # end copy
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
#my $theta = 104.48;
#my $theta_rad = 3.14159265359 * $theta / 360.0;
#for (my $rrr = 1.6; $rrr <= 2.61; $rrr += 0.2) {
#  push(@{$geometry{'geometries'}},
#       { 'description', sprintf("%8.5f", $rrr)."  ".sprintf("%8.5f", $rrr)."  ".sprintf("%6.2f", $theta),
#         'gnuplot_desc', "R = ".sprintf("%5.2f", $rrr).", theta = ".sprintf("%6.2f", $theta), # used in gnuplot files for keys
#         'atoms', [ [ "O", 0.0,  0.0,  0.0 ],
#                    [ "H", 0.0,  $rrr * sin($theta_rad), -$rrr * cos($theta_rad) ],
#                    [ "H", 0.0, -$rrr * sin($theta_rad), -$rrr * cos($theta_rad) ] ]
#       }
#  );
#}

# Here is an example of automatic generation of geometries 
# relative to a given (usually equilibrium) geometry specified above in %geometry.
# This example generates new geometries by stretching both OH bonds of the H2O molecule
# symmetrically by -0.4, -0.2, 0.0, 0.1, 0.2, 0.3 and 0.4 a.u.
# Note that this subroutine replaces even the reference geometry specified in %geometry,
# thus one should include in 'ranges' even stretching by 0.0 a.u.
# Modify 'geometry_labels' in %geometry accordingly, e.g. one could use " Change of O-H "
#
#&set_geometries_by_stretching_bonds(\%geometry,
#  {'ranges', [ [-0.4, -0.2, 0.2],
#               [ 0.0,  0.4, 0.1] ],
#   'shift', { '1-2', [2],             # this means: move the atom 2 (the first  hydrogen) in the direction of the bond O-H1
#              '1-3', [3] },           # this means: move the atom 3 (the second hydrogen) in the direction of the bond O-H2
#   'gnuplot_prefix', "Delta O-H = "
#  }
#);

