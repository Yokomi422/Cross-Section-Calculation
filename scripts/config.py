from dataclasses import dataclass


@dataclass
class UKRmolConfig:
    """
    Configuration related to config.pl
    """

    suffix: str = ""
    print_info = "file"
    add_files_to_backup: list[str] | None = None

    # Running options
    # Which codes will you run? Choose one or the other
    molpro: int = 1
    psi4: int = 0
    # What type of calculation/process?
    scattering: int = 1  # run all programs, if you want to run target only, set to 0
    photoionization: int = (
        0  # claculation photoionization cross sections of the (N+1)-electron system instead of electron scattering cross sections
    )
    rmt_interface: int = (
        0  # from the RMT molecular input file (rmt_interface and photoionization are mutually exclusive options)
    )
    # Ocassionally, you may want to turn these parts of the calculation too
    skip_radden: int = 1  # skip calculating radian densities
    # Options to keep or delete files
    gather_data: int = 1  # gather eigenphase sums, cross section, enegies...
    clean: int = 1  # removing fort* etc (except moints with molecuar integrals)
    remove_moints: int = 1  # remove moints with molecular integrals
    # Options to use existing input data
    use_templates: int = (
        1  # set this to 0, if you need to modify generated inputs manually and than rerun everything
    )
    # but do not change the filesnames1
    # MPI lancher options
    mpi_integrals: str = (
        ""  # MPI launcher command for scatci_integrals (e.g. mpirun -ilp64 -n 1)
    )
    mpi_scatci: str = (
        ""  # MPI launcher command for mpi-sctci (e.g. "mpirun -np 32; if empty, legacy scatci will be used")
    )
    mpi_rsolve: str = (
        ""  # MPI launcher command for mpi-rsolve (e.g. "mpirun -np 32; if empty, serial rsolve ill be used")
    )

    # Run-time options for SCATCI_INTEGRALS
    buffer_size: int = (
        5000  # Size of temporary arrays for integral transformation in MiB
    )
    delta_rl: float = (
        0.25  # Length, in Bohr, of the elementary radial quadrature needed for evaluation of the mixed BTO/GTO integrals
    )
    transform_alg: int = (
        0  # Choice of the integral transformation algorithm: 0 = auto, 1 = sparse (optimal for BTO and mixed BTO/GTO continuum), other = dense
    )
    # # The sparse transformation is not available in distributed (MPI) mode.

    # Saving options
    # These are files that can be plotted
    save_eigenph: int = (
        1  # set 1 to copy fort.10XX into a file like eigenph.singlet.Ag (eigenphase sums)
    )
    save_xsec: int = (
        1  # set 1 to copy fort.12XX into a file like xsec.singlet.Ag (cross sections)
    )
    # These are files that can be used to start an outer region calculation
    save_channels: int = (
        0  # set 1 to copy fort.10XX into a file like eigenph.singlet.Ag (eigenphase sums)
    )
    save_rmat_amp: int = (
        0  # set 1 to copy fort.12XX into a file like xsec.singlet.Ag (cross sections)
    )

    save_Kmatrix: int = 0
    save_Tmatrix: int = 0
    keep_inputs: int = 1
    keep_inputs: int = 1

    parallel_geom: int = 1
    parallel_symm: int = 1

    bound: int = 0
    use_cdenprop: int = 0
    ight: int = 2
    ighs: int = 2

    use_saved_ramps: int = 0
    run_eigenp: int = 1
    run_matrix: int = 1
    run_ixsecs: int = 1
    run_reson: int = 0
    run_time_delay: int = 0
    dipelm_smooth: int = 0

    def _format(self, value) -> str:
        """Format Python value into Perl format"""
        if isinstance(value, str):
            return f'"{value}"'
        elif isinstance(value, (int, float)):
            return str(value)
        elif isinstance(value, list):
            formatted_items = [self._format(item) for item in value]
            return f"[{', '.join(formatted_items)}]"
        else:
            return str(value)

    def _comment(self, key: str) -> str | None:
        """Add appropriate comment for each configuration key."""
        comments = {
            "suffix": "# string added to model directory to distinguish different runs which have the same model\n                         # it can be useful e.g. if you run only target calculations etc.",
            "print_info": '# "screen" or "file" or "both" or "none" - option whether info messages will be printed on the screen or into files in the subdirectory logs',
            "add_files_to_backup": '# all files given on the command line as "perl main.pl config.run.pl ..." are copied into directory script\n                                           # if you want to backup additional files, specify them here',
            "molpro": "# do SCF calculations using Molpro",
            "psi4": "# do SCF calculations using Psi4, a free alternative to Molpro",
            "scattering": "# run all programs, if you want to run target only, set to 0",
            "photoionization": "# calculate photoionization cross sections of the (N+1)-electron system instead of electron scattering cross sections",
            "rmt_interface": "# form the RMT molecular input file (rmt_interface and photoionization are mutually exclusive options)",
            "skip_radden": "# skip calculating radial densities",
            "gather_data": "# gather eigenphase sums, cross sections, energies ...",
            "clean": "# removing fort.* etc (except moints with molecular integrals)",
            "remove_moints": "# remove moints with molecular integrals",
            "use_templates": "# set this to 0, if you need to modify generated inputs manually and than rerun everything\n                         # but do not change the filenames !",
            "mpi_integrals": '# MPI launcher command for scatci_integrals (e.g. "mpirun -ilp64 -n 1")',
            "mpi_scatci": '# MPI launcher command for mpi-scatci (e.g. "mpirun -np 32"; if empty, legacy scatci will be used)',
            "mpi_rsolve": '# MPI launcher command for mpi-rsolve (e.g. "mpirun -np 32"; if empty, serial rsolve will be used)',
            "buffer_size": "# Size of temporary arrays for integral transformation in MiB",
            "delta_r1": "# Length, in Bohr, of the elementary radial quadrature needed for evaluation of the mixed BTO/GTO integrals",
            "transform_alg": "# Choice of the integral transformation algorithm: 0 = auto, 1 = sparse (optimal for BTO and mixed BTO/GTO continuum), other = dense\n                         # The sparse transformation is not available in distributed (MPI) mode.",
            "save_eigenph": "# set 1 to copy fort.10XX into a file like eigenph.singlet.Ag (eigenphase sums)",
            "save_xsec": "# set 1 to copy fort.12XX into a file like xsec.singlet.Ag (cross sections)",
            "save_channels": "# set 1 to copy fort.10  into a file like channels.singlet.Ag",
            "save_rmat_amp": "# set 1 to copy fort.21  into a file like ramps.singlet.Ag",
            "save_Kmatrix": "# set 1 to copy fort.9XX into a file like K-matrix.singlet.Ag",
            "save_Tmatrix": "# set 1 to copy fort.11XX into a file like T-matrix.singlet.Ag",
            "keep_inputs": "# set 1 to keep input files for UK R-matrix codes",
            "keep_outputs": "# set 1 to keep output files of UK R-matrix codes",
            "parallel_geom": "# number of geometries to be invoked in the same time",
            "parallel_symm": "# number of symmetries to be invoked in the same time (values <= 0 trigger alternative workflow using mpi-scatci)",
            "bound": "# set this option to 1 if you want to calculate bound states",
            "use_cdenprop": "# Force use of cdenprop for target calculation.",
            "ight": "# Set a value for igh. ight = for target-scatci, ighs = for scattering-scatci.\n  # 2 = Auto select (default), -1 = Arpack, 0 = Davidson, 1 = Givens-Householder",
            "ighs": "",
            "use_saved_ramps": "# If set to 1 the amplitudes and channel data saved from a previous run (using the options 'save_channels' and 'save_rmat_amp')\n  # will be used instead of running SWINTERF to generate them. This is useful in case the inner region data (e.g.fort.25) have\n  # not been saved only a rerun of the outer region is needed e.g. with a different energy grid or using a different set of programs.",
            "run_eigenp": "# calculate eigenphase sums for all symmetries",
            "run_tmatrx": "# calculate T-matrices for all symmetries",
            "run_ixsecs": "# calculate cross sections for all symmetries",
            "run_reson": "# calculate resonance fits for all symmetries",
            "run_time_delay": "# calculate time delays (requires a program that is not part of UKRmol+)",
            "dipelm_smooth": "# set: (0) for raw data, (1) for smoothed data or (2) for both.",
        }
        return comments.get(key, None)

    def transform(self):
        """
        Transform Python dataclass config to Perl config.
        """
        perl_lines = ["# Running and saving options", "", "%run = (", ""]

        # Group configurations by sections
        sections = [
            (
                "",
                ["suffix", "print_info", "add_files_to_backup"],
            ),
            ("Running options", ["molpro", "psi4"]),
            (
                "What type of calculation/process?",
                ["scattering", "photoionization", "rmt_interface"],
            ),
            (
                "Occasionally, you may want to turn these parts of the calculation too",
                ["skip_radden"],
            ),
            (
                "Options to keep or delete files",
                ["gather_data", "clean", "remove_moints"],
            ),
            ("Options to use existing input data", ["use_templates"]),
            ("MPI launcher options", ["mpi_integrals", "mpi_scatci", "mpi_rsolve"]),
            (
                "Run-time options for SCATCI_INTEGRALS",
                ["buffer_size", "delta_r1", "transform_alg"],
            ),
            (
                "Saving options",
                [
                    "save_eigenph",
                    "save_xsec",
                    "save_channels",
                    "save_rmat_amp",
                    "save_Kmatrix",
                    "save_Tmatrix",
                ],
            ),
            ("Keep files", ["keep_inputs", "keep_outputs"]),
            ("Parallelization", ["parallel_geom", "parallel_symm"]),
            ("For very special run to get bound states", ["bound"]),
            ("Expert settings", ["use_cdenprop", "ight", "ighs", "use_saved_ramps"]),
            (
                "Which programs to run",
                [
                    "run_eigenp",
                    "run_tmatrx",
                    "run_ixsecs",
                    "run_reson",
                    "run_time_delay",
                ],
            ),
            ("Photoionization", ["dipelm_smooth"]),
        ]

        for section, keys in sections:
            perl_lines.append(f"  # {section}")

            for key in keys:
                if hasattr(self, key):
                    value = getattr(self, key)
                    formatted = self._format(value)
                    comment = self._comment(value)

                    if comment is None:
                        perl_lines.append(f"  '{key}', {' ' * (18 - len(key))}{formatted}")
                    else:
                        perl_lines.append(f"  '{key}', {' ' * (18 - len(key))}{formatted}, {comment}")

        
            perl_lines.append("")

        perl_lines.append(");")
        return "\n".join(perl_lines)
    
config = UKRmolConfig()
perl_config = config.transform()

with open("config.pl", "w") as f:
    f.write(perl_config)


