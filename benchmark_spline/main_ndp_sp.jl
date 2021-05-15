## 動的計画法を使って政策関数を近似的に計算する
# チュートリアルセッション@日本経済学会 2021年度春季大会
# 価値関数の近似にスプライン補間(spline interpolation)を使用

# コードが置いてあるフォルダに移動
cd("/Users/TomoakiYamada/Documents/GitHub/Tutorial_JEA/benchmark_spline")

## 使用するパッケージを呼び出す
using Dierckx # スプライン補間
using Optim # 最適化(最小値を探す)
using Printf # @printfを使う
using Plots # 図を描く
pyplot() # 図のバックエンドとしてPyPlotを使用

@printf(" \n")
@printf("Solving dynamic programming problem using spline interpolation \n")

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
amin = 0.0 # グリッド(資産)の最小値
na = 101 # グリッドの数

## 関数を呼び出す

include("bellman_eq_sp.jl")
include("../benchmark_linear/crra.jl")
include("../benchmark_linear/grid_uni.jl")

## グリッドを生成

agrid = grid_uni(amin, amax, na)

## 価値関数・政策関数の初期値を設定

# 空の変数を用意
vfcn0 = zeros(na)
vfcn1 = zeros(na)
pfcn0 = zeros(na)
pfcn1 = zeros(na)

cons = (1 + r).*agrid .+ w
vfcn1 = crra.(cons, γ)

plot(agrid, vfcn1,
    xlab = "現在の資産保有量：a",
    ylab = "価値関数：V(a)",
#    title = "価値関数の初期値",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, 10),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_value_ini_sp.pdf")

## 最適化にはOptimというパッケージを使う
# 1回だけOptimizationを試してみる
# デモンストレーション用：実際はすぐにメインループに入ってOK

# Optimには様々な最適化の方法(オプション)がある
# 今回はGolden Section Searchという方法で最小値を探す
# 詳細：https://julianlsolvers.github.io/Optim.jl/stable/#

spl_coef = Spline1D(agrid, vfcn1, bc="extrapolate")

for i = 1:na
    # Golden Section Searchは両端を設定した上でその範囲内で極地(最小値)を探す方法
    # (-10,12)の範囲を指定：この範囲次第で極地が見つからなかったりするので、一般的には試行錯誤が必要な部分
    # あるいは速度を早くするために、別の最適化手法を試したり別のパッケージを探したり試行錯誤する部分
    res = optimize(x->bellman_eq_sp(x, agrid[i], spl_coef, β, γ, r, w), -10.0, 12.0, GoldenSection())
    pfcn0[i] = Optim.minimizer(res)
    vfcn0[i] = Optim.minimum(res)
end

plot(agrid, pfcn0,
    xlab = "現在の資産保有量：a",
    ylab = "次期の資産保有量：a'",
#    title = "政策関数：1回だけOptimization (スプライン近似)",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, 10),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_policy_optim_sp.pdf")

plot(agrid, -1*vfcn0,
    xlab = "現在の資産保有量：a",
    ylab = "価値関数：V(a)",
#    title = "価値関数：1回だけOptimization",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, 10),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_value_optim_sp.pdf")

## メインループ

maxit = 1000
vdif = zeros(maxit)
pdif = zeros(maxit)
ε1 = 1.0e-006
ε2 = 1.0e-006

# 計算速度を測る
start_time = time()

# 繰り返し計算誤差がε1、ε2未満になるまで繰り返し計算
for it = 1:maxit

    # スプラインの係数を予め計算しておく
    spl_coef = Spline1D(agrid, vfcn1, bc="extrapolate")

    for i = 1:na
        res = optimize(x->bellman_eq_sp(x, agrid[i], spl_coef, β, γ, r, w), -10.0, 15.0, GoldenSection())
        pfcn0[i] = Optim.minimizer(res)
        vfcn0[i] = Optim.minimum(res)
    end

    # ベルマン方程式(bellman_eq)の内部で一度-1をかけて符号を反転させている
    # これはOptimが最大値ではなく最小値を探しているため
    # そのため再度、符号を反転して元のベルマン方程式にする
    vfcn0 = -1*vfcn0

    # 計算誤差
    vdif[it] = maximum(abs.((vfcn0 .- vfcn1)./vfcn0))
    pdif[it] = maximum(abs.((pfcn0 .- pfcn1)./pfcn0))

    @printf("iteration counter:        %i \n", it)
    @printf("iteration error (value):  %f \n", vdif[it])
    @printf("iteration error (policy): %f \n", pdif[it])

    # 関数をアップデート
    vfcn1 = deepcopy(vfcn0)
    pfcn1 = deepcopy(pfcn0)

    if vdif[it] < ε1 || pdif[it] < ε2
        break
    end

end

time_elapsed = time() - start_time
@printf("")
@printf("Time elapsed: %0.3f \n", time_elapsed)

## 結果の図

#temp = range(1, maxit, length=maxit)
temp = range(1, 794, length=794) # 手抜き

plot(temp, vdif[1:794],
    xlabel = "計算回数",
    ylabel = "誤差",
#    title = "価値関数の繰り返し計算誤差",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, 800),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_diff_value_sp.pdf")

plot(temp, pdif[1:794],
    xlabel = "計算回数",
    ylabel = "誤差",
#    title = "政策関数の繰り返し計算誤差",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (0, 800),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_diff_policy_sp.pdf")

## 数値計算結果

plot(agrid, pfcn1,
    xlabel = "現在の資産保有量：a",
    ylabel = "次期の資産保有量：a' ",
#    title = "政策関数",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (amin, amax),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_policy_function_sp.pdf")

plot(agrid, vfcn1,
    xlab = "現在の資産保有量：a",
    ylab = "価値関数：V(a)",
#    title = "価値関数",
    linewidth = 4,
    color = :blue,
    legend = :none,
    framestyle = :semi,
    xlims = (amin, amax),
    titlefont = font("IPAexGothic",12),
    guidefont = font("IPAexGothic",12),
    legendfont = font("IPAexGothic",8),
    tickfont = font("IPAexGothic",8))
savefig("fig_value_function_sp.pdf")
