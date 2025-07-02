### What is this?

This folder contains all launching scripts for various programs like Gaussian etc.
It does not contain system-wide scripts from Polach (like `M12_SP`).

All scripts should work on all clusters in our laboratory.
This is achieved via script `SetEnvironment.sh`,
which determines the cluster and exports correct $PATHs and such.

### How to use it?
In your scripts, that you submit to queue,
use the following to get environmental variables from the `$USER`
    #$ -V
    SetEnvironment.sh MOLPRO [version]  

SetEnvironment.sh must be in `$PATH` 
Otherwise, the system won't find the scripts.


### Architecture of the scripts

Each program should have two separate scripts:
1. The script that actually executes the program, e.g. TERA
    It should:
    i)  set necessary variables such as scratch dirs etc.
    ii) run the program using input given from user

2. A launching script, e.g. launchGAUSS
   This script takes arguments from the user, checks them and creates
   a very simple auxiliary launching script, which is subsequently submitted.

Such separation has a big advantage if you need to run hundreds of calculations.
You can then easily pile them up into few jobs, i.e. you can have several G09 calculations in one job.
This strategy is used e.g. in script ../SPECTRA/LAUNCH/RecalcGeometries.sh.


### SCRATCH directories
They should look like this (if possible):
    export SCRDIR="/scratch/$USER/orca_$1_${JOB_ID}"

i.e. `program_jobname_jobid`

### Determining the cluster
Use the following code to determine the cluster name:

```bash
if [[ "$node" =~ ^s[0-9]+$|as67-1 ]];then
   cluster=as67
elif [[ "$node" =~ ^a[0-9]+$|403-a324-01 ]];then
   cluster=a324
elif [[ "$node" =~ ^n[0-9]+$|403-as67-01  ]];then
   cluster=as67gpu
else
   echo "I did not recognize any of the PHOTOX clusters. Please check the script $0"
   echo "Exiting..."
   exit 1
fi

if [[ "$cluster" -eq "as67" ]];then
   code for as67
else if [[ "$cluster" -eq "a324" ]];then
   code for a324
fi
```
