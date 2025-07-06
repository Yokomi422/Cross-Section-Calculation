from dataclasses import dataclass


@dataclass
class UKRmolConfig:
    """
    config.plに関連する設定
    """

    suffix: str = ""
    print_info = "file"
    add_files_to_backup: list[str] | None = None

    # 実行オプション
    # どのコードを実行しますか？どちらかを選択
    molpro: int = 1
    psi4: int = 0
    # どのタイプの計算/プロセス？
    scattering: int = 1  # run all programs, if you want to run target only, set to 0
    photoionization: int = (
        0  # 電子散乱断面積の代わりに(N+1)電子系の光イオン化断面積を計算
    )
    rmt_interface: int = (
        0  # RMT分子入力ファイルから(rmt_interfaceとphotoionizationは排他的オプション)
    )
    # 時々、計算のこれらの部分も有効にしたい場合があります
    skip_radden: int = 1  # 径方密度の計算をスキップ
    # ファイルを保持または削除するオプション
    gather_data: int = 1  # 固有位相和、断面積、エネルギーなどを収集
    clean: int = 1  # fort*などを削除（分子積分のmointを除く）
    remove_moints: int = 1  # 分子積分のmointを削除
    # 既存の入力データを使用するオプション
    use_templates: int = (
        1  # 生成された入力を手動で修正してからすべてを再実行する必要がある場合は0に設定
    )
    # ただしファイル名は変更しないでください！
    # MPIランチャーオプション
    mpi_integrals: str = (
        ""  # scatci_integrals用MPIランチャーコマンド (例: mpirun -ilp64 -n 1)
    )
    mpi_scatci: str = (
        ""  # mpi-scatci用MPIランチャーコマンド (例: "mpirun -np 32"; 空の場合は旧式scatciを使用)
    )
    mpi_rsolve: str = (
        ""  # mpi-rsolve用MPIランチャーコマンド (例: "mpirun -np 32"; 空の場合は連続rsolveを使用)
    )

    # SCATCI_INTEGRALSの実行時オプション
    buffer_size: int = (
        5000  # 積分変換用一時配列のサイズ(MiB単位)
    )
    delta_r1: float = (
        0.25  # 混合BTO/GTO積分の評価に必要な基本径方求積の長さ（Bohr単位）
    )
    transform_alg: int = (
        0  # 積分変換アルゴリズムの選択: 0=自動, 1=税密(BTOと混合BTO/GTO連続体に最適), その他=密
    )
    # # 税密変換は分散(MPI)モードでは利用できません。

    # 保存オプション
    # これらはプロット可能なファイルです
    save_eigenph: int = (
        1  # fort.10XXをeigenph.singlet.Agのようなファイルにコピーする場合は1に設定(固有位相和)
    )
    save_xsec: int = (
        1  # fort.12XXをxsec.singlet.Agのようなファイルにコピーする場合は1に設定(断面積)
    )
    # これらは外部領域計算を開始するために使用できるファイルです
    save_channels: int = (
        0  # fort.10XXをeigenph.singlet.Agのようなファイルにコピーする場合は1に設定(固有位相和)
    )
    save_rmat_amp: int = (
        0  # fort.12XXをxsec.singlet.Agのようなファイルにコピーする場合は1に設定(断面積)
    )

    save_Kmatrix: int = 0
    save_Tmatrix: int = 0
    keep_inputs: int = 1
    keep_outputs: int = 1

    parallel_geom: int = 1
    parallel_symm: int = 1

    bound: int = 0
    use_cdenprop: int = 0
    ight: int = 2
    ighs: int = 2

    use_saved_ramps: int = 0
    run_eigenp: int = 1
    run_tmatrix: int = 1
    run_ixsecs: int = 1
    run_reson: int = 0
    run_time_delay: int = 0
    dipelm_smooth: int = 0

    def _format(self, value) -> str:
        """Python値をPerl形式にフォーマット"""
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
        """各設定キーに適切なコメントを追加"""
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
            "run_tmatrix": "# calculate T-matrices for all symmetries",
            "run_ixsecs": "# calculate cross sections for all symmetries",
            "run_reson": "# calculate resonance fits for all symmetries",
            "run_time_delay": "# calculate time delays (requires a program that is not part of UKRmol+)",
            "dipelm_smooth": "# set: (0) for raw data, (1) for smoothed data or (2) for both.",
        }
        return comments.get(key, None)

    def transform(self):
        """
        Pythonデータクラス設定をPerl設定に変換。
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
                    "run_tmatrix",
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
                    comment = self._comment(key)

                    if comment is None:
                        perl_lines.append(f"  '{key}', {' ' * (18 - len(key))}{formatted}")
                    else:
                        perl_lines.append(f"  '{key}', {' ' * (18 - len(key))}{formatted}, {comment}")

        
            perl_lines.append("")

        perl_lines.append(");")
        return "\n".join(perl_lines)
    
