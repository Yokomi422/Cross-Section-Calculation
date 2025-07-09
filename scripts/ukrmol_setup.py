from dataclasses import dataclass
from typing import  Optional
from pathlib import Path

from model import UKRmolModel
from config import UKRmolConfig
from dirs import UKRmolDirs
from geometry import UKRmolGeometry

from ase.data import atomic_numbers

@dataclass
class UKRmolSetup:
    """
    UKRmol+計算の統合設定クラス
    model.py, config.py, dirs.py, geometry.pyを統合して使いやすくする
    """
    
    # 基本設定
    molecule_name: str = "Ne"
    atoms: list[list] | None = None
    nelectrons: int = 10
    symmetry: str = "C1"
    
    # 計算設定
    basis: str = "cc-pVTZ"
    model_type: str = "CAS"
    scattering: bool = True
    photoionization: bool = False
    max_energy_eV: float = 20.0  # 衝突断面積の最大エネルギー[eV]
    
    # 軌道設定
    nfrozen: int = 1
    nactive: int = 4
    nvirtual: int = 8  # デフォルト: 2 * nactive
    nreference: int = 5  # デフォルト: nfrozen + nactive
    
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
            self.atoms = [[self.molecule_name, 0.0, 0.0, 0.0]]
        
        # ASEを使って電子数を自動計算
        if self.atoms:
            total_electrons = 0
            for atom in self.atoms:
                symbol = atom[0]
                total_electrons += atomic_numbers[symbol]
            self.nelectrons = total_electrons
        
        # 軌道設定のデフォルト値を自動計算
        # nreferenceとnvirtualは明示的に設定されていない場合のみ自動計算
        if self.nreference == 5:  # デフォルト値の場合
            self.nreference = self.nfrozen + self.nactive
        if self.nvirtual == 8:  # デフォルト値の場合
            self.nvirtual = 2 * self.nactive
    
    def create_model(self) -> UKRmolModel:
        """モデル設定を作成"""
        if self.custom_model is not None:
            return self.custom_model
        
        # nescat設定を自動計算（20eV/400の比率で計算）
        nescat_value = int(self.max_energy_eV / 20.0 * 400)
        nescat_str = f"20, {nescat_value}"
        
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
            nreference=self.nreference,
            r_unit=self.length_unit,
            directory=self.suffix,
            nescat=nescat_str
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
    
    def generate_all_files(self, output_path: str = ".") -> dict[str, str]:
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
    return UKRmolSetup(
        molecule_name=atom_symbol,
        atoms=[[atom_symbol, 0.0, 0.0, 0.0]],
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
