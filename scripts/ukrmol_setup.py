from dataclasses import dataclass
from typing import List, Dict, Any, Optional
from pathlib import Path
import os

from model import UKRmolModel
from config import UKRmolConfig
from dirs import UKRmolDirs
from geometry import UKRmolGeometry


@dataclass
class UKRmolSetup:
    """
    UKRmol+計算の統合設定クラス
    model.py, config.py, dirs.py, geometry.pyを統合して使いやすくする
    """
    
    # 基本設定
    molecule_name: str = "Ne"
    atoms: List[List] = None
    nelectrons: int = 10
    symmetry: str = "C1"
    
    # 計算設定
    basis: str = "cc-pVTZ"
    model_type: str = "CAS"
    scattering: bool = True
    photoionization: bool = False
    
    # 軌道設定
    nfrozen: int = 1
    nactive: int = 4
    nvirtual: int = 3
    
    # 幾何学設定
    length_unit: int = 0  # 0 - atomic units, 1 - Angstroms
    
    # 出力設定
    output_dir: str = "sample_output"
    suffix: str = ""
    
    # 高度な設定（オプション）
    custom_model: Optional[UKRmolModel] = None
    custom_config: Optional[UKRmolConfig] = None
    custom_dirs: Optional[UKRmolDirs] = None
    custom_geometry: Optional[UKRmolGeometry] = None
    
    def __post_init__(self):
        """初期化後の処理"""
        if self.atoms is None:
            # デフォルト原子設定
            if self.molecule_name == "He":
                self.atoms = [["He", 0.0, 0.0, 0.0]]
                self.nelectrons = 2
            elif self.molecule_name == "Ne":
                self.atoms = [["Ne", 0.0, 0.0, 0.0]]
                self.nelectrons = 10
            elif self.molecule_name == "H":
                self.atoms = [["H", 0.0, 0.0, 0.0]]
                self.nelectrons = 1
            else:
                self.atoms = [[self.molecule_name, 0.0, 0.0, 0.0]]
    
    def create_model(self) -> UKRmolModel:
        """モデル設定を作成"""
        if self.custom_model is not None:
            return self.custom_model
        
        return UKRmolModel(
            molecule=self.molecule_name,
            atoms=[atom[0] for atom in self.atoms],
            nelectrons=self.nelectrons,
            symmetry=self.symmetry,
            basis=self.basis,
            model=self.model_type,
            nfrozen=self.nfrozen,
            nactive=self.nactive,
            nvirtual=self.nvirtual,
            r_unit=self.length_unit,
            directory=self.suffix
        )
    
    def create_config(self) -> UKRmolConfig:
        """設定を作成"""
        if self.custom_config is not None:
            return self.custom_config
        
        return UKRmolConfig(
            suffix=self.suffix,
            scattering=1 if self.scattering else 0,
            photoionization=1 if self.photoionization else 0
        )
    
    def create_dirs(self) -> UKRmolDirs:
        """ディレクトリ設定を作成"""
        if self.custom_dirs is not None:
            return self.custom_dirs
        
        return UKRmolDirs(output=self.output_dir)
    
    def create_geometry(self) -> UKRmolGeometry:
        """幾何学設定を作成"""
        if self.custom_geometry is not None:
            return self.custom_geometry
        
        return UKRmolGeometry(
            atoms=self.atoms,
            length_unit=self.length_unit,
            suffix=f".{self.molecule_name}_{self.suffix}" if self.suffix else f".{self.molecule_name}"
        )
    
    def generate_all_files(self, output_path: str = ".") -> Dict[str, str]:
        """すべての設定ファイルを生成"""
        model = self.create_model()
        config = self.create_config()
        dirs = self.create_dirs()
        geometry = self.create_geometry()
        
        files = {
            "model.pl": model.transform(),
            "config.pl": config.transform(),
            "dirs.pl": dirs.transform(),
            "geometry.pl": geometry.transform()
        }
        
        # ファイル書き込み
        output_dir = Path(output_path)
        output_dir.mkdir(exist_ok=True)
        
        for filename, content in files.items():
            file_path = output_dir / filename
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(content)
        
        return files
    
    def print_summary(self):
        """設定の概要を表示"""
        print(f"分子: {self.molecule_name}")
        print(f"原子: {self.atoms}")
        print(f"電子数: {self.nelectrons}")
        print(f"対称性: {self.symmetry}")
        print(f"基底: {self.basis}")
        print(f"モデル: {self.model_type}")
        print(f"散乱計算: {self.scattering}")
        print(f"光イオン化: {self.photoionization}")
        print(f"軌道設定: frozen={self.nfrozen}, active={self.nactive}, virtual={self.nvirtual}")
        print(f"出力ディレクトリ: {self.output_dir}")


# 便利な関数
def create_atom_setup(atom_symbol: str, **kwargs) -> UKRmolSetup:
    """原子計算用の設定を作成"""
    electron_count = {
        "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
        "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
        "S": 16, "Cl": 17, "Ar": 18
    }
    
    return UKRmolSetup(
        molecule_name=atom_symbol,
        atoms=[[atom_symbol, 0.0, 0.0, 0.0]],
        nelectrons=electron_count.get(atom_symbol, 1),
        **kwargs
    )


def create_h2_setup(bond_length: float = 1.4, **kwargs) -> UKRmolSetup:
    """H2分子計算用の設定を作成"""
    return UKRmolSetup(
        molecule_name="H2",
        atoms=[
            ["H", 0.0, 0.0, 0.0],
            ["H", 0.0, 0.0, bond_length]
        ],
        nelectrons=2,
        symmetry="D2h",
        **kwargs
    )


def create_co_setup(bond_length: float = 2.13, **kwargs) -> UKRmolSetup:
    """CO分子計算用の設定を作成"""
    return UKRmolSetup(
        molecule_name="CO",
        atoms=[
            ["C", 0.0, 0.0, 0.0],
            ["O", 0.0, 0.0, bond_length]
        ],
        nelectrons=14,
        symmetry="C2v",
        **kwargs
    )


# 使用例
if __name__ == "__main__":
    # He原子の設定
    he_setup = create_atom_setup("He")
    he_setup.print_summary()
    print("\n" + "="*50 + "\n")
    
    # Ne原子の設定（カスタム）
    ne_setup = UKRmolSetup(
        molecule_name="Ne",
        basis="cc-pVDZ",
        nfrozen=1,
        nactive=8,
        nvirtual=2,
        suffix="test_run"
    )
    ne_setup.print_summary()
    print("\n" + "="*50 + "\n")
    
    # H2分子の設定
    h2_setup = create_h2_setup(bond_length=1.4)
    h2_setup.print_summary()
    print("\n" + "="*50 + "\n")
    
    # すべてのファイルを生成
    files = he_setup.generate_all_files("he_output")
    print(f"生成されたファイル: {list(files.keys())}")