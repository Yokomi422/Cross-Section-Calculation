from dataclasses import dataclass
from typing import List, Dict, Any


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
                self.suffix = f"   {atom_name}    "
                self.geometry_labels = f"   {atom_name}  atom"
            else:
                # デフォルト: 原点にNe原子
                self.geometries = [
                    {
                        'description': "   Ne atom  ",
                        'gnuplot_desc': "Ne atom",
                        'atoms': [["Ne", 0.0, 0.0, 0.0]]
                    }
                ]

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
            "# Automatic geometry generation example",
            "# Neon atom geometry - single atom at origin",
            "push(@{$geometry{'geometries'}},",
            "     { 'description', \"   Ne atom  \",",
            "       'gnuplot_desc', \"Ne atom\",",
            "       'atoms', [ [ \"Ne\", 0.0, 0.0, 0.0 ] ]",
            "     }",
            ");",
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
