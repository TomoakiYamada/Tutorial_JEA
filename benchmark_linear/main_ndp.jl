## 動的計画法を使って政策関数を近似的に計算する
# チュートリアルセッション@日本経済学会

# コードが置いてあるフォルダに移動
cd("/Users/TomoakiYamada/Documents/GitHub/Tutorial_JEA/benchmark_linear")

## 使用するパッケージを呼び出す
using Interpolations # 線形補間
using Optim # 最適化(最小値)
using Printf # @printf
using Plots # 図を描く
pyplot() # 図のバックエンドとしてPyPlotを使用

@printf(" \n")
@printf("Solving dynamic programming problem \n")

## calibration

β = 0.99 # discount factor
γ = 2.0 # inverse of ies
r = 0.01 # interest rate
w = 1.0 # wage

@printf("discount factor:        %f \n", β)
@printf("relative risk aversion: %f \n", γ)
@printf("interest rate           %f \n", r)
@printf("wage                    %f \n", w)
@printf(" ")

# グリッドの範囲を指定
# 今回は特に深いことは考えないで適当に範囲指定
amax = 10.0 # グリッド(資産)の最大値
amin = 0.0 # グリッド(資産)の最大値
na = 101 # グリッドの数

## 関数を呼び出す

include("bellman_eq.jl")
include("crra.jl")
include("grid_uni.jl")

## グリッドを生成

agrid = grid_uni(amin, amax, na)

# 以下のように書くことも可能
# ただし、この場合はagridの型がVectorではなくなる点に注意
#agrid = range(amin, amax, length=na)
#typeof(agrid)

## 価値関数の初期値を設定

# 空の変数を用意
vfcn0 = zeros(na)
vfcn1 = zeros(na)
pfcn0 = zeros(na)
pfcn1 = zeros(na)

# T期のvalueがゼロならT-1期の最適貯蓄はゼロになる
# そこからT-1期前のvalueを計算する
# for i = 1:na
#     cons = (1 + r)*agrid[i] + w - 0.0 # consume total wealth
#     vfcn1[i] = crra(cons, γ)
# end

# 上の書き方でも良いけど関数のブロードキャストという方法を使えばスッキリ書く事ができる
# .をつけることでVector全体を一括で計算
# 上の書き方だと一つ一つのagridの要素を取り出して計算していた点に注意
cons = (1 + r).*agrid .+ w
vfcn1 = crra.(cons, γ)

## デバッグ用の図
# あくまでデモンストレーション用なので、普通はここまで逐一チェックをする必要は必ずしもない
# 日本語を使いたい場合、フォントを指定する必要がある

plot(agrid, vfcn1,
    xlab = "現在の資産保有量：a",
    ylab = "価値関数：V(a)",
    title = "価値関数のInitial Guess",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, 10),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_value_ini_li.pdf")

## 最適化にはOptimを使う
# 1回だけOptimizationを試してみる
# デモンストレーション用：実際はすぐにメインループに入ってOK

# Optimには様々な最適化の方法がある
# 今回はGolden Section Searchという方法で最小値を探す
# 詳細：https://julianlsolvers.github.io/Optim.jl/stable/#
for i = 1:na
    res = optimize(x->bellman_eq(x, agrid[i], agrid, vfcn1, β, γ, r, w), -10.0, 15.0, GoldenSection())
    pfcn0[i] = Optim.minimizer(res)
    vfcn0[i] = Optim.minimum(res)
end

plot(agrid, pfcn0,
    xlab = "現在の資産保有量：a",
    ylab = "次期の資産保有量：a'",
    title = "価値関数：1回だけOptimization (線形近似)",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, 10),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_value_optim_li.pdf")

plot(agrid, -1*vfcn0,
    xlab = "現在の資産保有量：a",
    ylab = "価値関数：V(a)",
    title = "価値関数：1回だけOptimization",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, 10),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_value_optim_li.pdf")

## 1000回繰り返し計算をしてみる

maxit = 1000
vdif = zeros(maxit)
pdif = zeros(maxit)

# とりあえず1000回計算してみる
for it = 1:maxit

    for i = 1:na
        # 最適化する変数はa'のみで、それ以外の変数(agrid[i]以降)はすべて最適化の際には所与
        res = optimize(x->bellman_eq(x, agrid[i], agrid, vfcn1, β, γ, r, w), -10.0, 15.0, GoldenSection())
        pfcn0[i] = Optim.minimizer(res)
        vfcn0[i] = Optim.minimum(res)
    end

    # ベルマン方程式(bellman_eq)の内部で一度-1をかけて符号を反転させている
    # これはOptimが最大値ではなく最小値を探しているため
    # そのため再度、符号を反転して元のベルマン方程式にする
    vfcn0 = -1*vfcn0

    # 計算誤差
    vdif[it] = maximum(abs.(vfcn0 - vfcn1))
    pdif[it] = maximum(abs.(pfcn0 - pfcn1))

    @printf("iteration counter:        %i \n", it)
    @printf("iteration error (value):  %f \n", vdif[it])
    @printf("iteration error (policy): %f \n", pdif[it])

    # 関数をアップデート
    vfcn1 = deepcopy(vfcn0)
    pfcn1 = deepcopy(pfcn0)

end

## デバッグ用の図

temp = range(1, maxit, length=maxit)

plot(temp, vdif,
    xlabel = "計算回数",
    ylabel = "誤差",
    title = "価値関数の繰り返し計算誤差",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, maxit),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_diff_value_li.pdf")

plot(temp, pdif,
    xlabel = "計算回数",
    ylabel = "誤差",
    title = "政策関数の繰り返し計算誤差",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, maxit),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_diff_policy_li.pdf")

## 数値計算結果

plot(agrid, pfcn1,
    xlabel = "現在の資産保有量：a",
    ylabel = "次期の資産保有量：a' ",
    title = "政策関数",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (amin, amax),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_policy_function_li.pdf")

plot(agrid, vfcn1,
    xlab = "現在の資産保有量：a",
    ylab = "価値関数：V(a)",
    title = "価値関数",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (amin, amax),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_value_function_li.pdf")
