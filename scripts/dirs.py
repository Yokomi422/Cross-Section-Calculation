from dataclasses import dataclass


@dataclass
class UKRmolDirs:
    """
    UKRmol+計算用のディレクトリ設定
    dirs.plに対応
    """
    
    # UKRmol+実行ファイル
    bin_in: str = "/home/pj25000080/ku50001398/ukrmol+/ukrmol-in-3.2/build/bin"
    bin_out: str = "/home/pj25000080/ku50001398/ukrmol+/ukrmol-out-3.2/build/bin"
    
    # 量子化学計算プログラム
    molpro: str = "/home/app/Molpro/2024.1.0_mpipr/bin/"
    psi4: str = ""
    molcas: str = ""
    
    # リソースディレクトリ
    basis: str = "/home/pj25000080/ku50001398/ukrmol+/projects/resources/basis.sets"
    templates: str = "/home/pj25000080/ku50001398/ukrmol+/projects/resources/input.templates"
    libs: str = "../resources/lib"
    
    # 出力ディレクトリ
    output: str = "sample_output"

    def _format(self, value) -> str:
        """Python値をPerl形式にフォーマット"""
        if isinstance(value, str):
            return f'"{value}"'
        else:
            return str(value)

    def transform(self):
        """
        PythonデータクラスをPerlハッシュ形式に変換
        """
        perl_lines = [
            "# Settings related to a specific installation of UK R-matrix codes",
            "",
            "# Path to executables can be specified directly in %dirs or later if you are using several computers",
            "# (see switch($run{'computer'}) below",
            "# Use ${bs} (see dirfile.pm) in relative paths instead of \\ or / for portability",
            "# Some directories are determined later automatically",
            "# If the full path must be used then it is not necessary to use ${bs}",
            "",
            "%dirs = (",
        ]
        
        # Add directory entries
        entries = [
            ("bin_in", self.bin_in, "UKRmol+ inner region executables"),
            ("bin_out", self.bin_out, "UKRmol+ outer region executables"),
            ("molpro", self.molpro, "Molpro installation directory"),
            ("psi4", self.psi4, "Psi4 installation directory"),
            ("molcas", self.molcas, "Molcas installation directory"),
            ("basis", self.basis, "Basis sets directory"),
            ("templates", self.templates, "Input templates directory"),
            ("libs", self.libs, "Libraries directory"),
            ("output", self.output, "Output directory"),
        ]
        
        for key, value, _ in entries:
            formatted_value = self._format(value)
            perl_lines.append(f" '{key}',{' ' * (12 - len(key))}{formatted_value},")
        
        perl_lines.append(")")
        
        return "\n".join(perl_lines)
