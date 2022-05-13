# テーブル定義書 <!-- omit in toc -->

全体概要．

- [実体関連モデル](#実体関連モデル)
- [CRUD表](#crud表)
- [テーブルごとの定義](#テーブルごとの定義)
    - [entity_1](#entity_1)
    - [entity_2](#entity_2)
    - [entity_3](#entity_3)

## 実体関連モデル

![実体関連モデル](./figures/entity-relationship.svg)

## CRUD表

| 機能 | `entity_1` | `entity_2` | `entity_2` |
|:--|:-:|:-:|:-:|
| 機能1 | C | R | UD |

## テーブルごとの定義

### entity_1

概要．

| カラム名 | 概要 | データ型 | バイト長 | 非NULL | ユニーク | 備考 |
| --- | --- | --- | --- | --- | --- | --- |
| `entity_1_id` | サロゲートキー | `serial4` | 4 | x | x | `PRIMARY KEY` |
| `attribute_1` | 単精度浮動小数の属性 | `float4` | 4 ||||
| `note` | 備考 | `varchar(32)` | 32 ||||

#### 索引

- `entity_1_id`: b-tree
    - ユニーク制約．検索・結合用．

#### 備考

### entity_2

概要．

| カラム名 | 概要 | データ型 | バイト長 | 非NULL | ユニーク | 備考 |
| --- | --- | --- | --- | --- | --- | --- |
| `entity_2_id` | サロゲートキー | `serial4` | 4 | x | x | `PRIMARY KEY` |
| `attribute_2` | 8バイト整数の属性 | `int8` | 8 ||||
| `note` | 備考 | `varchar(32)` | 32 ||||

#### 索引

- `entity_2_id`: b-tree
    - ユニーク制約．検索・結合用．

#### 備考

### entity_3

概要．

| カラム名 | 概要 | データ型 | バイト長 | 非NULL | ユニーク | 備考 |
| --- | --- | --- | --- | --- | --- | --- |
| `entity_3_id` | サロゲートキー | `serial4` | 4 | x | x | `PRIMARY KEY` |
| `entity_1_id` | 対応する実体1 | `int4` | 4 | x || `REFERENCES entity_1(entity_1_id)` |
| `entity_2_id` | 対応する実体2 | `int4` | 4 | x || `REFERENCES entity_1(entity_1_id)` |
| `attribute_3` | 倍精度浮動小数の属性 | `float8` | 8 ||||

#### 索引

- `entity_3_id`: b-tree
    - ユニーク制約．
- `entity_1_id`: b-tree
    - 検索・結合用．
- `entity_2_id`: b-tree
    - 検索・結合用．

#### 備考
