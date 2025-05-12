#!/bin/bash

# メッシュ生成
blockMesh
surfaceFeatures
snappyHexMesh -overwrite

## 初期条件の設定
# cp 0.orig/* 0/

# シミュレーションの実行とログの保存
scalarSimpleFoam | tee log

## 結果の可視化
#paraFoam
