# Settings related to a specific installation of UK R-matrix codes

# Path to executables can be spicified directly in %dirs or later if you are using several computers
# (see switch($run{'computer'}) below
# Use ${bs} (see dirfile.pm) in relative paths instead of \ or / for portability
# Some directories are determined later automatically
# If the full path must be used then it is not necessary to use ${bs}

%dirs = (
 'bin_in',    "/home/pj25000080/ku50001398/ukrmol+/ukrmol-in-3.2/build/bin",
 'bin_out',   "/home/pj25000080/ku50001398/ukrmol+/ukrmol-out-3.2/build/bin",
  'molpro',    "/home/app/Molpro/2024.1.0_mpipr/bin/",
  'psi4',      "",
  'molcas',    "",
  'basis',     "/home/pj25000080/ku50001398/ukrmol+/projects/resources/basis.sets",
  'templates', "/home/pj25000080/ku50001398/ukrmol+/projects/resources/input.templates",
  'use_tmple_output",
  'libs',      "../resources/lib",
  'output',    "sample_output",
)
