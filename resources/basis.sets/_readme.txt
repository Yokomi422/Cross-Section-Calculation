The basis sets (except swmol3.continuuum... basis sets) in this directory 
were downloaded as Molpro format from EMSL Basis Set Exchange, 
https://bse.pnl.gov/bse/portal and they are ready to use. 

If you need another basis, download it from EMSL in Molpro format and 
copy it in the basis.sets directory on your computer.

Names for some basis were slightly changed because of filename conventions:

6-31Gxx    is actually 6-31G**
6-311Gxx   is actually 6-311G**
6-311ppGxx is actually 6-311++G**

etc.

Continuum basis sets, such as swmol3.continuum.r10.L4, were generated using numcbas and gtobas codes
and their format was used for old quantum chemistry code 'swmol3' which is no longer supported.
Input for the newer scatci_integrals code will be generated automatically from these files by the scripts.
An example of bash scripts to generate such basis are in tools/continuum.basis.


