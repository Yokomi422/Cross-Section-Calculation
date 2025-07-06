# UKRmol - UK R行列法による分子散乱計算

## 概要

UKRmolは、電子-分子および陽電子-分子散乱計算のためのUK R行列法の実装です。このプロジェクトは、分子の電子散乱断面積、光イオン化断面積、共鳴位置と幅、およびターゲット状態の特性を計算するための包括的なツールセットを提供します。

## 主な機能

- **電子散乱計算**: 弾性散乱および非弾性散乱断面積の計算
- **光イオン化**: 分子の光イオン化断面積の計算
- **共鳴解析**: 共鳴位置と幅の同定
- **複数の計算モデル**: SE、SEP、CAS、CHF、MASなど
- **並列計算**: 幾何構造と対称性に対する並列実行をサポート

## ディレクトリ構造

```
UKRmol/
├── basis.sets/          # 量子化学基底関数セット
│   ├── sto-3g/         # 最小基底
│   ├── 6-31g/          # ポプル基底
│   ├── cc-pvtz/        # 相関一貫基底
│   └── continuum/      # 連続状態基底
├── input.templates/     # 計算用テンプレートファイル
│   ├── molpro.inp      # Molpro入力テンプレート
│   ├── psi4.inp        # Psi4入力テンプレート
│   ├── scatci_integrals.inp  # 積分変換
│   ├── congen.inp      # 配置生成
│   └── scattering.*.inp # 散乱計算
├── projects/           # 分子別計算プロジェクト
│   ├── Ne/            # ネオン
│   ├── He/            # ヘリウム
│   ├── Ar/            # アルゴン
│   ├── Kr/            # クリプトン
│   ├── Xe/            # キセノン
│   └── H2O/           # 水分子
├── resources/          # ライブラリとユーティリティ
│   ├── lib/           # Perlモジュール
│   └── scripts/       # 自動化スクリプト
└── photoxrepo/        # 追加ツールとBEBモデル実装
```

## 前提条件

### 必須ソフトウェア

- **量子化学パッケージ** (以下のいずれか):
  - Molpro
  - Psi4
  - Gaussian
  - ORCA
  - Molcas
- **R行列コード**:
  - SCATCI (散乱計算)
  - CONGEN (配置生成)
  - RSOLVE (外部領域)
- **プログラミング言語**:
  - Perl 5.x
  - Python 3.x (プロット用)
- **並列計算** (オプション):
  - MPI実装 (OpenMPI, MPICH等)

## インストール

1. リポジトリのクローン:
```bash
git clone https://github.com/yourusername/UKRmol.git
cd UKRmol
```

2. 環境変数の設定:
```bash
export UKRMOL_HOME=/path/to/UKRmol
export PATH=$UKRMOL_HOME/resources/scripts:$PATH
```

3. 必要なPerlモジュールのインストール:
```bash
cpan install Config::General
cpan install Math::Complex
```

## 使用方法

### 基本的なワークフロー

1. **プロジェクトディレクトリの作成**:
```bash
cd projects/
mkdir my_molecule
cd my_molecule
```

2. **設定ファイルの作成**:
   - `geometry.pl`: 分子構造の定義
   - `model.pl`: 計算モデルのパラメータ
   - `config.pl`: 実行設定
   - `main.pl`: メイン実行スクリプト

3. **計算の実行**:
```bash
perl main.pl
```

### 設定例 (H2O分子)

`geometry.pl`:
```perl
$geometry = {
    atoms => [
        { element => 'O',  x =>  0.000,  y => 0.000,  z => 0.000 },
        { element => 'H',  x =>  0.757,  y => 0.586,  z => 0.000 },
        { element => 'H',  x => -0.757,  y => 0.586,  z => 0.000 }
    ],
    charge => 0,
    multiplicity => 1
};
```

`model.pl`:
```perl
$model = {
    type => 'SEP',              # Static Exchange + Polarization
    basis => 'cc-pvtz',         # 基底関数セット
    active_space => [8, 4],     # CAS(8,4)
    states => 5,                # 計算する状態数
    energy_range => [0, 20],    # エネルギー範囲 (eV)
};
```

## 計算モデル

- **SE (Static Exchange)**: 静的交換近似
- **SEP (Static Exchange + Polarization)**: 分極効果を含む
- **CAS (Complete Active Space)**: 完全活性空間法
- **CHF (Coupled Hartree-Fock)**: 結合ハートリー・フォック法
- **MAS (Multiple Active Spaces)**: 複数活性空間アプローチ

## 出力ファイル

計算後、以下のファイルが生成されます:

- `cross_sections.dat`: 散乱断面積データ
- `eigenphases.dat`: 固有位相シフト
- `resonances.dat`: 共鳴パラメータ
- `photoionization.dat`: 光イオン化断面積

## プロット

結果の可視化:
```bash
python plot.py --input cross_sections.dat --output plot.png
```

## トラブルシューティング

### よくある問題

1. **メモリ不足エラー**:
   - `config.pl`でメモリ割り当てを増やす
   - より小さい基底関数セットを使用

2. **収束しない場合**:
   - 初期推定を改善
   - 収束基準を緩和

3. **並列実行の問題**:
   - MPIが正しくインストールされているか確認
   - ノード間通信を確認

## 貢献

プルリクエストを歓迎します。大きな変更の場合は、まずissueを開いて変更内容を議論してください。

## ライセンス

[ライセンスタイプを指定]

## 引用

このソフトウェアを使用した場合は、以下を引用してください:
```
[適切な引用情報を追加]
```

## 連絡先

質問や問題がある場合は、[メールアドレス]までご連絡ください。