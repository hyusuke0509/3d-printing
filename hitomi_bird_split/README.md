# hitomi_bird_egg_keycap 白い羽 分離

元データ: `hitomi_bird_egg_keycap.3mf`（Bambu Studio、白=フィラメント1 / 黒=フィラメント2 で塗り分け）

| ファイル | 内容 | 体積 |
|---|---|---|
| body.stl | 本体（卵・首・頭・黒い胴体） | 20244 mm³ |
| tail.stl | 白い尾羽 | 1008 mm³ |
| wing_R.stl | 右翼先の白い羽 (-Y側) | 317 mm³ |
| wing_L.stl | 左翼先の白い羽 (+Y/-X側、卵に接している部分) | 143 mm³ |
| hitomi_bird_split.obj/.mtl | 上記4パーツを1ファイルに（Blender用、オブジェクト分割済み・白/黒/切断面マテリアル付き） | |

- 全パーツ水密（watertight）。座標は元の位置のまま（底面 z=0, mm）。組み合わせると元形状に一致。
- 分割方法: 3MFの `paint_color` を解読し、黒く塗られていない羽の白領域を抽出 → 境界ループを膜状の面で塞いで羽パーツを閉じ、本体側にも同じ面（反転）を付与。左翼の羽は卵と白領域が連続していたため、卵の球面フィット（中心(-2.93,-0.90,-35.64) R17.58, 元座標）から0.3mm以上外側の部分を羽とした。
- 切断面は曲面。接着前提。羽パーツはサポートが必要な可能性あり。
- Blender: File > Import > Wavefront (.obj) で `hitomi_bird_split.obj`。
- `scripts/` は再現用（parse → decode → comps → lwing → split → export の順）。
