* dark_analysis.cc
  input: 生波形が入ったroot file, output: 生波形を積分した結果を詰めたroot file
  生波形をハイブリッド法(疋田修論3.3.2節)により積分する。
  depth, divisionなどのパラメータは適宜変更する必要あり。
* dark_prop.cc
  input: dark_analysis.rootの出力root file
  1 p.e.ゲインやクロストーク・アフターパルス確率、ダークカウントレートを計算する。
  クロストーク・アフターパルス確率、ダークカウントレートについては、ある程度ずれてしまうことがわかっているため注意(疋田修論3.5.1節)
* LED_analysis.cc
  input: LED信号の生波形が入ったroot file, output: LED信号波形を積分した結果を詰めたroot file
  LED信号の生波形をLED信号区間だけ積分する。
  LED信号のstart, endは適宜調節する必要あり。
* saturation.cc
  input: LED_analysis.rootの出力root file
  信号強度の異なるLED測定の結果から、saturation curveをプロットし、回復時間を計算する。
  光量変換の際に、測定で使用した10倍アンプの正確な値を使用していたり、測定時にLED光量が多いところと少ないところでアンプの増幅率を変えていたりするため、この部分は要調整。
