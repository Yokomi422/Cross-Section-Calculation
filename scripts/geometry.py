from dataclasses import dataclass
from typing import List, Dict, Any, Optional
try:
    from ase import Atoms
    from ase.build import molecule
    from ase.data import atomic_numbers, chemical_symbols
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False


@dataclass
class UKRmolGeometry:
    """
    UKRmol+計算用の幾何学設定
    geometry.plに対応
    """
    
    # 幾何学制御
    suffix: str = ".Ne_atom"
    geometry_labels: str = "   Ne atom   "
    correct_cm: int = 1
    length_unit: int = 0  # 0 - atomic units, 1 - Angstroms
    
    # 幾何学処理制御
    start_at_geometry: int = 1
    stop_at_geometry: int = 0
    
    # 幾何学データ - 幾何学辞書のリスト
    geometries: List[Dict[str, Any]] | None = None
    
    # 直接原子指定
    atoms: List[List] | None = None
    
    # ASE分子名での自動生成
    molecule_name: Optional[str] = None
    
    # 希ガス原子の自動生成
    noble_gas: Optional[str] = None
    
    def __post_init__(self):
        if self.geometries is None:
            if self.atoms is not None:
                # 原子リストから幾何学を作成
                atom_name = self.atoms[0][0] if self.atoms else "Ne"
                self.geometries = [
                    {
                        'description': f"   {atom_name} atom  ",
                        'gnuplot_desc': f"{atom_name} atom",
                        'atoms': self.atoms
                    }
                ]
                self.suffix = f".{atom_name}_atom"
                self.geometry_labels = f"   {atom_name} atom   "
            elif self.molecule_name is not None:
                # ASE分子名から自動生成
                self._generate_from_ase_molecule()
            elif self.noble_gas is not None:
                # 希ガス原子の自動生成
                self._generate_noble_gas()
            else:
                # デフォルト: 原点にNe原子
                self.geometries = [
                    {
                        'description': "   Ne atom  ",
                        'gnuplot_desc': "Ne atom",
                        'atoms': [["Ne", 0.0, 0.0, 0.0]]
                    }
                ]

    def _generate_from_ase_molecule(self):
        """ASE分子名から幾何学構造を自動生成"""
        if not ASE_AVAILABLE:
            raise ImportError("ASE is not available. Please install ASE to use molecule_name feature.")
        
        try:
            mol = molecule(self.molecule_name)
            symbols = mol.get_chemical_symbols()
            positions = mol.get_positions()
            
            # 原子リストを作成
            atoms_list = []
            for symbol, pos in zip(symbols, positions):
                atoms_list.append([symbol, float(pos[0]), float(pos[1]), float(pos[2])])
            
            # 分子の説明を作成
            unique_symbols = list(set(symbols))
            if len(unique_symbols) == 1:
                desc = f"   {unique_symbols[0]}{len(symbols)} molecule  "
            else:
                desc = f"   {self.molecule_name} molecule  "
            
            self.geometries = [
                {
                    'description': desc,
                    'gnuplot_desc': f"{self.molecule_name} molecule",
                    'atoms': atoms_list
                }
            ]
            self.suffix = f".{self.molecule_name}_molecule"
            self.geometry_labels = f"   {self.molecule_name} molecule   "
            
        except Exception as e:
            raise ValueError(f"Failed to generate molecule {self.molecule_name}: {e}")

    def _generate_noble_gas(self):
        """希ガス原子の幾何学構造を自動生成"""
        noble_gases = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn']
        
        if self.noble_gas not in noble_gases:
            raise ValueError(f"Invalid noble gas: {self.noble_gas}. Valid options: {noble_gases}")
        
        self.geometries = [
            {
                'description': f"   {self.noble_gas} atom  ",
                'gnuplot_desc': f"{self.noble_gas} atom",
                'atoms': [[self.noble_gas, 0.0, 0.0, 0.0]]
            }
        ]
        self.suffix = f".{self.noble_gas}_atom"
        self.geometry_labels = f"   {self.noble_gas} atom   "

    @classmethod
    def from_molecule(cls, molecule_name: str, **kwargs):
        """ASE分子名から幾何学構造を作成するクラスメソッド"""
        return cls(molecule_name=molecule_name, **kwargs)
    
    @classmethod
    def from_noble_gas(cls, noble_gas: str, **kwargs):
        """希ガス原子から幾何学構造を作成するクラスメソッド"""
        return cls(noble_gas=noble_gas, **kwargs)

    def _format(self, value) -> str:
        """Python値をPerl形式にフォーマット"""
        if isinstance(value, str):
            return f'"{value}"'
        elif isinstance(value, (int, float)):
            return str(value)
        elif isinstance(value, list):
            if len(value) == 0:
                return "[]"
            formatted_items = []
            for item in value:
                if isinstance(item, dict):
                    # Format geometry dictionary
                    formatted_items.append(self._format_geometry_dict(item))
                else:
                    formatted_items.append(self._format(item))
            return f"[{', '.join(formatted_items)}]"
        elif isinstance(value, dict):
            return self._format_geometry_dict(value)
        else:
            return str(value)

    def _format_geometry_dict(self, geom_dict: Dict[str, Any]) -> str:
        """幾何学辞書をPerl用にフォーマット"""
        lines = ["\n       # start copy"]
        lines.append(f"       {{ 'description', {self._format(geom_dict['description'])}, # string to use in output files")
        lines.append(f"         'gnuplot_desc', {self._format(geom_dict['gnuplot_desc'])}, # used in gnuplot files")
        lines.append("         # specify ALL atoms (even redundant with respect to symmetry elements)")
        
        atoms_str = "         'atoms', [ "
        for i, atom in enumerate(geom_dict['atoms']):
            atom_str = f"[ {self._format(atom[0])}, {atom[1]:12.6f}, {atom[2]:12.6f}, {atom[3]:12.6f} ]"
            if i < len(geom_dict['atoms']) - 1:
                atoms_str += atom_str + ",\n                    "
            else:
                atoms_str += atom_str + " ]"
        
        lines.append(atoms_str)
        lines.append("       },")
        lines.append("       # end copy")
        
        return "\n".join(lines)

    def transform(self):
        """
        Pythonデータクラスをperlハッシュ形式に変換
        """
        perl_lines = [
            "# Geometries can be specified either manually in the array 'geometries'",
            "# by copying what is between start copy and end copy for each geometry",
            "# or automatically as shown below the hash array %geometry",
            "#",
            "# Note that 'length_unit' should be set to the actual length unit used for values in 'atoms' array,",
            "# later (in &generate_geometries) $model{'r_unit'} is set also to this unit for consistency reasons",
            "#",
            "# Note that there's a restriction on the orientation of the molecules in the scripts",
            "# Your molecular target should be oriented as follows for the corresponding point groups:",
            "# D2h: any orientation",
            "# C2v: C2 axis along Z",
            "# C2h: C2 axis along Z",
            "# D2:  any orientation",
            "# C2:  C2 axis along Z",
            "# Cs:  Molecule in the YZ plane",
            "# Ci:  any orientation",
            "# C1:  any orientation",
            "",
            "%geometry = (",
            "",
            f"  'suffix', {self._format(self.suffix)},   # string added to model directory",
            "                              # with different geometry settings",
            "",
            f"  'geometry_labels',  {self._format(self.geometry_labels)},         # labels used on the first line of output files",
            "                                          # it should correspond to numbers given in 'geometries'->'description'",
            f"  'correct_cm',   {self.correct_cm},                      # correct the center of mass to be at the origin",
            f"  'length_unit',  {self.length_unit},                      # 0 - atomic units, 1 - Angstroms",
            "",
            "  'geometries', [",
        ]
        
        # Add geometries
        for geom in self.geometries:
            perl_lines.append(self._format_geometry_dict(geom))
        
        perl_lines.extend([
            "  ],",
            "",
            "  # the following option can be used e.g. to continue an interrupted run",
            "  # all geometries are generated but codes will run only for specified geometries",
            f"  'start_at_geometry', {self.start_at_geometry},                # codes will run only for geometries with index >= than this number",
            f"  'stop_at_geometry',  {self.stop_at_geometry},               # this can be used to stop at certain geometry",
            "                                         # if zero then codes will run for all geometries",
            ");",
            "",
        ])
        
        return "\n".join(perl_lines)


# 使用例
if __name__ == "__main__":
    # He原子を指定
    he_geometry = UKRmolGeometry(atoms=[["He", 0.0, 0.0, 0.0]])
    print("He原子の設定:")
    print(he_geometry.transform())
    
    # 複数原子（水素分子）
    h2_geometry = UKRmolGeometry(
        atoms=[
            ["H", 0.0, 0.0, 0.0],
            ["H", 0.0, 0.0, 1.4]
        ],
        suffix=".H2_molecule",
        geometry_labels="   H2 molecule   "
    )
    print("\nH2分子の設定:")
    print(h2_geometry.transform())
    
    # デフォルト（Ne原子）
    default_geometry = UKRmolGeometry()
    print("\nデフォルト（Ne原子）の設定:")
    print(default_geometry.transform())
    
    # ASE分子名から自動生成
    print("\n" + "="*50)
    print("ASE自動生成の例:")
    
    # 希ガス原子（クラスメソッド使用）
    ar_geometry = UKRmolGeometry.from_noble_gas("Ar")
    print("\nAr原子の設定:")
    print(ar_geometry.transform())
    
    if ASE_AVAILABLE:
        # H2分子（ASE使用）
        h2_ase_geometry = UKRmolGeometry.from_molecule("H2")
        print("\nH2分子（ASE）の設定:")
        print(h2_ase_geometry.transform())
        
        # H2O分子（ASE使用）
        h2o_geometry = UKRmolGeometry.from_molecule("H2O")
        print("\nH2O分子（ASE）の設定:")
        print(h2o_geometry.transform())
        
        # CO2分子（ASE使用）
        try:
            co2_geometry = UKRmolGeometry.from_molecule("CO2")
            print("\nCO2分子（ASE）の設定:")
            print(co2_geometry.transform())
        except:
            print("\nCO2分子は利用できません")
    else:
        print("\nASEが利用できません")
