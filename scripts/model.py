from dataclasses import dataclass
from typing import  Any, Optional


@dataclass
class UKRmolModel:
    """
    UKRmol+計算用のモデル設定
    model.plに対応
    """
    
    # ディレクトリと基本情報
    directory: str = ""
    molecule: str = "Ne"
    atoms: list[str] | None = None
    nelectrons: int = 10
    symmetry: str = "C1"
    
    # 単位
    r_unit: int = 0  # 0 - atomic units, 1 - Angstroms
    e_unit: int = 2  # 0 - atomic units, 1 - Rydbergs, 2 - eV
    x_unit: int = 1  # 1 - atomic units, 2 - Angstroms^2
    
    # 基底とモデル
    basis: str = "cc-pVTZ"
    orbitals: str = "natural"
    charge_of: str = "target"
    select_orb_by: str = "molden"
    model: str = "CAS"
    
    # 軌道空間
    nfrozen: int = 1
    nactive: int = 4
    nvirtual: int = 3
    nreference: int = 5
    
    # 対称性別軌道配列
    frozen_orbs: list[int] | None = None
    active_orbs: list[int] | None = None
    virtual_orbs: list[int] | None = None
    reference_orbs: list[int] | None = None
    
    # MASアプローチ
    use_MASSCF: int = 0
    use_MAS: int = 0
    qchem_MAS: Optional[Any] = None
    qchem_constraints: Optional[Any] = None
    MAS: list[Any] | None= None
    constraints: Optional[Any] | None = None
    l2_MAS: list[Any] | list = None
    l2_constraints: Optional[Any] = None
    
    # ターゲット状態
    ncasscf_states: dict[str, list[int]] | None = None
    ntarget_states: dict[str, list[int]] | None = None
    ntarget_states_used: str = "2"
    
    # 削除闾値
    delthres: list[str] = None
    
    # PCOオプション
    use_PCO: int = 0
    reduce_PCO_CAS: int = 1
    maxl_PCO: int = 2
    PCO_alpha0: list[float] = None
    PCO_beta: list[float] = None
    num_PCOs: list[int] = None
    PCO_gto_thrs: list[float] = None
    PCO_delthres: list[str] = None
    
    # 散乱状態
    scattering_states: dict[str, list[int]] = None
    
    # 連続基底
    use_GTO: int = 1
    radius_GTO: int = 10
    maxl_GTO: int = 4
    use_BTO: int = 0
    start_BTO: float = 8.0
    order_BTO: int = 9
    no_of_BTO: int = 12
    maxl_BTO: int = 4
    maxl_legendre_1el: int = 70
    maxl_legendre_2el: int = 55
    rmatrix_radius: float = 10.0
    
    # 伝搬
    max_multipole: int = 2
    raf: float = 70.0
    
    # エネルギーグリッド
    nescat: str = "20, 400"
    einc: str = "0.00015, 0.0005, 0.001, 0.05"
    
    # 断面積
    maxi: str = "1"
    maxf: str = "0"
    
    # 散乱タイプ
    positron_flag: int = 0
    
    # 光イオン化
    initialsym: str = "1"
    first_Ip: float = 0.0
    
    # 出力オプション
    norbitals_to_print: int = 10
    nstates_to_print: int = 10
    
    # メモリ制御
    lndo: int = 10000000
    
    def __post_init__(self):
        if self.atoms is None:
            self.atoms = ["Ne"]
        if self.frozen_orbs is None:
            self.frozen_orbs = [1, 0, 0, 0, 0, 0, 0, 0]
        if self.active_orbs is None:
            self.active_orbs = [4, 0, 0, 0, 0, 0, 0, 0]
        if self.virtual_orbs is None:
            self.virtual_orbs = [0, 0, 0, 0, 0, 0, 0, 0]
        if self.reference_orbs is None:
            self.reference_orbs = [0, 0, 0, 0, 0, 0, 0, 0]
        if self.MAS is None:
            self.MAS = [8, [8, 8], 'active']
        if self.l2_MAS is None:
            self.l2_MAS = [
                [8, [8, 9], 'active'],
                [5, [0, 1], 'active']
            ]
        if self.ncasscf_states is None:
            self.ncasscf_states = {
                'singlet': [1, 0, 0, 0, 0, 0, 0, 0],
                'triplet': [0, 0, 0, 0, 0, 0, 0, 0]
            }
        if self.ntarget_states is None:
            self.ntarget_states = {
                'singlet': [2, 0, 0, 0, 0, 0, 0, 0],
                'triplet': [0, 0, 0, 0, 0, 0, 0, 0]
            }
        if self.delthres is None:
            self.delthres = ["1.0D-07"] * 8
        if self.PCO_alpha0 is None:
            self.PCO_alpha0 = [0.15, 0.15, 0.15]
        if self.PCO_beta is None:
            self.PCO_beta = [1.3, 1.3, 1.3]
        if self.num_PCOs is None:
            self.num_PCOs = [1, 1, 1]
        if self.PCO_gto_thrs is None:
            self.PCO_gto_thrs = [-1.0, -1.0, -1.0]
        if self.PCO_delthres is None:
            self.PCO_delthres = ["1.0D-06"] * 8
        if self.scattering_states is None:
            self.scattering_states = {
                'doublet': [1, 0, 0, 0, 0, 0, 0, 0],
                'quartet': [0, 0, 0, 0, 0, 0, 0, 0]
            }

    def _format(self, value) -> str:
        """Python値をPerl形式にフォーマット"""
        if isinstance(value, str):
            return f'"{value}"'
        elif isinstance(value, (int, float)):
            return str(value)
        elif isinstance(value, list):
            if len(value) == 0:
                return "[]"
            formatted_items = [self._format(item) for item in value]
            return f"[{', '.join(formatted_items)}]"
        elif isinstance(value, dict):
            items = []
            for k, v in value.items():
                items.append(f"'{k}' => {self._format(v)}")
            return f"{{{', '.join(items)}}}"
        elif value is None:
            return "undef"
        else:
            return str(value)

    def transform(self):
        """
        PythonデータクラスをPerlハッシュ形式に変換
        """
        perl_lines = [
            "# Settings which are necessary for UK R-matrix codes input files",
            "#          except geometries (set in %geometry in a different file)",
            "",
            "%model = (",
            "",
            f"  'directory',     {self._format(self.directory)}, # directory (relative path) which should correspond to a given model",
            "                       # if it is empty, it will be set later according to the model settings",
            "",
            "  # Molecule",
            "",
            f"  'molecule',    {self._format(self.molecule)},               # Used only for directory and descriptions",
            f"  'atoms',        {self._format(self.atoms)},             # All atoms must be specified",
            f"  'nelectrons',    {self.nelectrons},               # number of target electrons",
            f"  'symmetry',     {self._format(self.symmetry)},              # Point group symmetry",
            "",
            "  # Units",
            "",
            f"  'r_unit',        {self.r_unit},                 # distance unit: 0 - atomic units, 1 - Angstroms",
            f"  'e_unit',        {self.e_unit},                 # energy unit:   0 - atomic units, 1 - Rydbergs, 2 - eV",
            f"  'x_unit',        {self.x_unit},                 # cross section: 1 - atomic units, 2 - Angstroms^2",
            "",
            "  # Model - each model will have its own directory from these settings",
            f"  'basis',         {self._format(self.basis)},         # basis set",
            f"  'orbitals',      {self._format(self.orbitals)},         # which orbitals to use",
            f"  'charge_of',      {self._format(self.charge_of)},         # use orbitals for N-electron or (N+1)-electron system",
            f"  'select_orb_by',  {self._format(self.select_orb_by)},         # molden = use orbitals as ordered in Molden file",
            f"  'model',             {self._format(self.model)},         # model type",
            f"  'nfrozen',       {self.nfrozen},                 # number of frozen target orbitals",
            f"  'nactive',       {self.nactive},                 # number of active target orbitals",
            f"  'nvirtual',      {self.nvirtual},                 # number of virtual orbitals",
            f"  'nreference',    {self.nreference},                 # number of orbitals used for searching reference orbitals",
            "",
            "  # Orbital arrays per symmetry",
            f"  'frozen_orbs',   {self._format(self.frozen_orbs)}, # which orbitals for each symmetry to use as frozen",
            f"  'active_orbs',   {self._format(self.active_orbs)}, # which orbitals for each symmetry to use as active",
            f"  'virtual_orbs',  {self._format(self.virtual_orbs)}, # which orbitals for each symmetry to use as virtual",
            f"  'reference_orbs',{self._format(self.reference_orbs)}, # which orbitals for each symmetry to use for reference",
            "",
            "  # Multiple Active Spaces (MAS) approach",
            f"  'use_MASSCF', {self.use_MASSCF},  # use MASSCF approach",
            f"  'use_MAS',    {self.use_MAS},  # use MAS approach",
            f"  'qchem_MAS', {self._format(self.qchem_MAS)}, # quantum chemistry MAS",
            f"  'qchem_constraints', {self._format(self.qchem_constraints)}, # quantum chemistry constraints",
            f"  'MAS', {self._format(self.MAS)},",
            f"  'constraints', {self._format(self.constraints)}, # additional constraints",
            f"  'l2_MAS', {self._format(self.l2_MAS)},",
            f"  'l2_constraints', {self._format(self.l2_constraints)}, # l2 constraints",
            "",
            "  # Target states",
            f"  'ncasscf_states', {self._format(self.ncasscf_states)},",
            f"  'ntarget_states', {self._format(self.ntarget_states)},",
            f"  'ntarget_states_used', {self._format(self.ntarget_states_used)},",
            "",
            "  # Deletion thresholds",
            f"  'delthres',      {self._format(self.delthres)},",
            "",
            "  # PCO options",
            f"  'use_PCO',               {self.use_PCO},",
            f"  'reduce_PCO_CAS',        {self.reduce_PCO_CAS},",
            f"  'maxl_PCO',              {self.maxl_PCO},",
            f"  'PCO_alpha0',   {self._format(self.PCO_alpha0)},",
            f"  'PCO_beta',     {self._format(self.PCO_beta)},",
            f"  'num_PCOs',     {self._format(self.num_PCOs)},",
            f"  'PCO_gto_thrs', {self._format(self.PCO_gto_thrs)},",
            f"  'PCO_delthres',      {self._format(self.PCO_delthres)},",
            "",
            "  # Scattering states",
            f"  'scattering_states', {self._format(self.scattering_states)},",
            "",
            "  # Continuum basis set",
            f"  'use_GTO',             {self.use_GTO},",
            f"  'radius_GTO',         {self.radius_GTO},",
            f"  'maxl_GTO',            {self.maxl_GTO},",
            f"  'use_BTO',             {self.use_BTO},",
            f"  'start_BTO',         {self.start_BTO},",
            f"  'order_BTO',           {self.order_BTO},",
            f"  'no_of_BTO',          {self.no_of_BTO},",
            f"  'maxl_BTO',            {self.maxl_BTO},",
            f"  'maxl_legendre_1el',  {self.maxl_legendre_1el},",
            f"  'maxl_legendre_2el',  {self.maxl_legendre_2el},",
            f"  'rmatrix_radius',   {self.rmatrix_radius},",
            "",
            " # Propagation step",
            f"  'max_multipole',       {self.max_multipole},",
            f"  'raf',              {self.raf},",
            "",
            " # Energy grid",
            f"  'nescat',      {self._format(self.nescat)},",
            f"  'einc',     {self._format(self.einc)},",
            "",
            " # Initial and final states",
            f"  'maxi',              {self._format(self.maxi)},",
            f"  'maxf',              {self._format(self.maxf)},",
            "",
            "  # Scattering type",
            f"  'positron_flag',       {self.positron_flag},",
            "",
            "  # Photoionization settings",
            f"  'initialsym',        {self._format(self.initialsym)},",
            f"  'first_Ip',          {self.first_Ip},",
            "",
            "  # Printing options",
            f"  'norbitals_to_print', {self.norbitals_to_print},",
            f"  'nstates_to_print',   {self.nstates_to_print},",
            "",
            "  # Memory control",
            f"  'lndo',          {self.lndo},",
            ");"
        ]
        
        return "\n".join(perl_lines)
