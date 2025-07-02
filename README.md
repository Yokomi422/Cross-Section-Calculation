# Programs setup
if you want to run ukrmol+ scirpts, please setup as following
```bash
# move to ukrmol+ directory
cd /home/pj25000080/ku50001398/ukrmol+
```

Before you execute ukrmol scripts you have to copy some files into your working directory.
1. /scripts
2. /lib

After you have set up these above, to run ukrmol scripts please run the command below. 

before running the scripts, you have to set properties of your target molecule on model.pl, geometry.pl
```bash
perl main.pl dirs.pl config.pl model.pl geometry.pl
```

# error handling in main.pl
## Can't locate ForkManager.pm
```bash
an't locate ForkManager.pm in @INC (you may need to install the ForkManager module) (@INC contains: /home/pj2│ 21 
5000080/ku50001398/ukrmol+/projects/He/../lib /home/pj25000080/ku50001398/perl5/lib/perl5/5.26.3/x86_64-linux-│ 22 
thread-multi /home/pj25000080/ku50001398/perl5/lib/perl5/5.26.3
6_64-linux-thread-multi /home/pj25000080/ku50001398/perl5/lib/perl5 /usr/local/lib64/perl5 /usr/local/share/pe│~                                                                                                                               
rl5 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5) at main.pl li│~                                                                                                                               
ne 23.    
```
main.pl:23 
```
- use ForkManager;
+ use Parallel::ForkManager;
```

## Cat't locate ukrmollib.pm
`
`bash
Can't locate ukrmollib.pm in @INC (you may need to install the ukrmollib module) (@INC contains: /home/pj25000│
080/ku50001398/ukrmol+/projects/He/../lib /home/pj25000080/ku50001398/perl5/lib/perl5/5.26.3/x86_64-linux-thre│
ad-multi /home/pj25000080/ku50001398/perl5/lib/perl5/5.26.3 /home/pj25000080/ku50001398/perl5/lib/perl5/x86_64│
-linux-thread-multi /home/pj25000080/ku50001398/perl5/lib/perl5 /usr/local/lib64/perl5 /usr/local/share/perl5 │
/usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5) at main.pl line 2│
4.                                                                                                            │
BEGIN failed--compilation aborted at main.pl line 24. 
```
somehow, cannot find modules on the same directory, so change the path in the main.pl directory.
- about line 21, ../lib -> ../resources/lib

## cannot make directory no permissions
```perl
use File::Path qw(make_path);
```
- not effective
- make_dir -> make_path

change path
```perl
$dirs{'output'} = $dirs{'cwd'} . '/output' unless defined $dirs{'output'} && $dirs{'output'} ne '';
nfo'} = 'both';warn ">>> output dir = $dirs{'output'}\n";
warn ">>> dirs{output} = '$dirs{output}'\n";
warn ">>> model{directory} = '$model{directory}'\n";
warn ">>> about to make dirs{model} = '$dirs{model}'\n";


## Died at /home/pj25000080/ku50001398/ukrmol+/projects/He/../resources/lib/ukrmollib.pm line 628.

in which_continuum_basis_file_to_use subroutine, 
reason: basis.sets' path is inaccurate

fix below
```bash
unless ( defined $dirs{'basis'} && $dirs{'basis'} ne '' ) {
  $dirs{'basis'} = "$dirs{'cwd'}${bs}../resources/basis.sets";
}
```
## Could not open he.molden: No such file or directory at /home/pj25000080/ku50001398/ukrmol+/projects/He/../resources/lib/MultiSpace.pm line 1450.



## display print_info
```perl
$run{'print_info'} = 'both';
```
## sh: /molpro: No such file or directory
default path setting should be done

## Warning: no template file for molpro.inp !

## molpro input must be in pwd
input file means target.molpro.inp, not inputs directory
mpi program is not working, maybe because multi theaeds are on unknown pwds
-> no mpi

* fehler on processor = maybe multi process issues?

in ukrmollib.pl not relative path, not . included
```
my $input = basename($r_par->{'data'}->{'inputfile'});
my $output = basename($r_par->{'data'}->{'outputfile'});
```

## run_code
- not include .
- not variable, but constant string

# important functions
## run_code
