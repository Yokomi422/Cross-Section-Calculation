# BEB模型による電子衝突イオン化断面積計算

## 概要

Binary Encounter Bethe（BEB）模型は、分子の電子衝突イオン化断面積を計算するための理論的枠組みです。このリポジトリには、量子化学計算を処理してイオン化断面積を予測するPython実装（`MISC/BEB/beb.py`）が含まれています。

## 理論

BEB模型は以下を組み合わせています：
- Binary encounter理論（古典力学）
- Bethe理論（量子力学）

各分子軌道に対して2つの重要なパラメータが必要です：
- **B**: 電子結合エネルギー（イオン化エネルギー）
- **U**: 軌道電子の平均運動エネルギー

- B: 結合エネルギー = イオン化エネルギー
- U: 軌道運動エネルギー

## Gaussianでの使用方法

### ステップ1: Gaussian計算の実行

軌道運動エネルギーを出力するためにIOP(6/81=3)オプションを含むGaussian入力ファイルを作成します：

```
%mem=550Mb
%NProcShared=1
#HF/aug-cc-pVTZ Symmetry=None OPT IOP(6/81=3)

分子の説明

0 1
[座標]
```

**重要**: IOP(6/81=3)オプションが軌道運動エネルギーを出力するため、BEB模型に必要な全ての情報（軌道エネルギーと運動エネルギー）を一度の計算で取得できます。

### ステップ2: イオン化断面積の計算

```bash
cd photoxrepo/MISC/BEB
./beb.py -i gaussian_output.log --Tmax 1000 -m beb > cross_section.out
```

これにより、第一イオン化閾値から1000 eVまでの入射電子エネルギーに対するイオン化断面積が生成されます。

### 出力フォーマット

出力は以下を含みます：
- 第1列: 入射電子エネルギー [eV]
- 第2列: 全イオン化断面積 [Å²]
- 第3列以降: 個別軌道の寄与 [Å²]

## Gaussianの優位性

Gaussianは以下の理由でBEB計算に最適です：

1. **一貫性のある計算**: 構造最適化と軌道情報取得を同じ計算条件で実行
2. **運動エネルギーの直接出力**: IOP(6/81=3)により正確なUパラメータを取得
3. **シンプルなワークフロー**: 一つのプログラムで全ての必要な情報を出力
4. **高い信頼性**: 複数プログラム間でのデータ受け渡しエラーを回避

## 高度な使用法

### 特定軌道の計算

単一軌道の断面積を計算：
```bash
./beb.py -N 2 -U 15.98 -B 15.43 -T 100 -m beb
```

パラメータ：
- N: 軌道内の電子数
- U: 軌道運動エネルギー [eV]
- B: 軌道結合エネルギー [eV]
- T: 入射電子エネルギー [eV]

### 荷電種

一価陽イオンの場合、電荷フラグを追加：
```bash
./beb.py -i gaussian.log --Tmax 1000 -m beb -c 1 > cation_cross_section.out
```

## 代替手法: Talukder模型

単純な原子系では、Talukder経験模型が利用できます：
```bash
./beb.py -N 2 -B 15.0 -T 100 -n 1 -l 0 -m talukder
```

ここで、nとlは主量子数と方位量子数です。

## 制限事項

1. 現在の実装は閉殻分子のみに対応
2. 精度は量子化学計算の品質に依存

## 例：水分子の計算ワークフロー

**Gaussianでの完全なワークフロー:**

1. **Gaussian入力ファイルの作成 (water_hf_avtz.com):**
   ```
   %mem=550Mb
   %NProcShared=1
   #HF/aug-cc-pVTZ Symmetry=None OPT IOP(6/81=3)

   水分子のBEB計算用最適化

   0 1
   O    0.000000    0.000000    0.000000
   H    0.757000    0.000000    0.586000
   H   -0.757000    0.000000    0.586000
   ```

2. **Gaussian計算の実行:**
   ```bash
   g09 < water_hf_avtz.com > water_hf_avtz.log
   ```

3. **BEB計算の実行:**
   ```bash
   cd photoxrepo/MISC/BEB
   ./beb.py -i water_hf_avtz.log --Tmax 1000 -m beb > water_cross_section.out
   ```

4. **結果のプロット:**
   ```bash
   xmgrace water_cross_section.out
   ```

この方法により、たった3つのステップで構造最適化から断面積計算まで完了します。

## 参考文献

1. Kim & Irikura (2000) - Electron-impact ionization cross sections for polyatomic molecules
2. Kim & Rudd (1994) - Binary-encounter-dipole model for electron-impact ionization
3. Talukder et al. (2008) - Empirical model for electron impact ionization cross sections

## ディレクトリ構造

```
photoxrepo/
├── MISC/
│   └── BEB/
│       ├── beb.py          # メインBEB実装
│       ├── README.md       # 基本的な使用法
│       └── SAMPLES/        # 計算例
└── INPUTS/
    ├── MOLPRO/            # Molpro例
    └── G09/               # Gaussian例
```
