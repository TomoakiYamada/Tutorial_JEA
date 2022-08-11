# 日本経済学会 2021年度春季大会 チュートリアルセッション

## Tutorial_JEA.pdf
* チュートリアルセッションで使用したスライド

## Install_Julia.pdf
* 一応残してありますが**情報が古い**ので注意してください
* Juliaをインストールする方法を解説したスライド
  * チュートリアルの動画内では使用していません
* スライド内のアンダーラインはハイパーリンクになっています

## Install_Julia_2nd.pdf
* チュートリアルセッション時はJuliaProをオススメしていたのですが、更新が止まってダウンロードが出来なくなったため、VS CodeからJuliaを使う方法を簡潔に紹介しています
* それ以外の内容はほぼ変更ありません

## コードを動かすために必要なパッケージ
```
using Pkg

Pkg.add("Dierckx")
Pkg.add("IJulia")
Pkg.add("Interpolations")
Pkg.add("Optim")
Pkg.add("Plots")
```
## benchmark_spline
* 価値関数の近似にSpline interplationを使用
* スライドで使われている図はすべてここから生成

## benchmark_linear
* 価値関数の近似にlinear interpolationを使用
* 上のSpline近似と比較して途中の政策関数がギザギザしたり(fig_value_optim_li.pdf)、収束に時間がかかる(fig_diff_value_li.pdf)
* 計算時間もかかっているけど、これはSpline近似はSplineの係数を1回だけ計算しているのに対して、線形補間の場合には毎回傾きを計算しているため

## Tutorial_JEA.ipynb
* Jupyter notebook
* benchmark_splineの内容について解説を加えたもの