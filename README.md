# euler-solver-go

一维可压缩 Euler 方程的本地有限体积求解器，支持含**激波**与**稀疏波**的守恒律问题、
**自适应网格细化/合并（AMR）**、**周期/反射边界**、**自适应时间步**、**逐位可复现**
以及**快照持久化与崩溃恢复**。纯 Go 实现，无第三方依赖。

## 1. 控制方程与变量

守恒律

```
∂U/∂t + ∂F(U)/∂x = 0,
U = [ρ, ρu, E]ᵀ,
F = [ρu, ρu²+p, u(E+p)]ᵀ,
E = p/(γ−1) + ½ρu²        （理想气体，默认 γ=1.4）
```

* 守恒变量 `Conserved{ρ, ρu, E}` 与原始变量 `Primitive{ρ, u, p}` 提供**双向转换**
  （`euler/state.go`）。
* 任何 `ρ ≤ 0` 或导出 `p ≤ 0` 的状态都以 `ErrNonPhysical` **拒绝**（初始化、分裂、
  两个 RK 阶段均校验）。
* 单元存储的是守恒变量**体积平均**；质量/动量/能量积分为 `Σ U·Δx`。

## 2. 空间离散与激波捕获

* 有限体积：`dU_i/dt = −(F_{i+1/2} − F_{i−1/2})/Δx_i`。
* 界面通量可选：
  * `rusanov`（默认，局部 Lax–Friedrichs）——最稳健、正性友好；
  * `hllc`（Toro 的接触分辨 HLLC + PVRS 波速估计）。
* 高阶重构：等宽三单元 stencil 上的 **MUSCL** 斜率重构，配 **TVD 限制器**
  （`minmod` 默认，或 `vanleer`），因此跨激波/接触/稀疏波**不产生新的极值、无吉布斯振荡**。
* 粗细网格交界（hanging node）界面与反射壁面**强制一阶**（单元均值/镜像对），
  保证该界面只有一个公共 Riemann 通量，跨层级严格守恒且不振荡。
* 间断位置与峰值随网格加密收敛（见 `TestSodConvergesUnderRefinement`，密度 L1 误差
  在 128/256/512 网格上单调下降）。

## 3. 边界条件

* `periodic`：越界索引回绕；域切口两侧呈现为**同一个 Riemann 问题**，界面通量逐位相消，
  全局质量/动量/能量在整段积分中保持到机器精度（约 1e-16）。
* `reflecting`：鬼元 `(ρ,−u,p)` 精确镜像，界面质量与能量通量对称地为零，壁面不泄漏、
  不注入、不振荡（壁面通过对称冲量改变动量，这是物理结果）。

边界单元被锁定在最粗级别，保证周期切口两侧几何完全一致、反射 stencil 干净对称。

## 4. 自适应网格（AMR）

* 单元用位树编码（`id = 1` 后每位代表左/右子节点，物理根级别 `RootLevel=24`）。
  `N0` 必须是 2 的幂且 ≤ 2^24；配置中的 `MinLevel/MaxLevel` 是相对初始网格的层数
  （`MaxLevel=3` 即最多 8 倍细化）。
* **细化判据（可配置）**：尺度无关的跳变指示子
  `score = (2|fc−fl|+2|fc−fr|)/(|fl|+2|fc|+|fr|+ε(1+|fc|))`，
  可针对密度（`density-gradient`）或压力（`pressure-gradient`）。
  光滑区随网格 `O(Δx)` 减小，间断处为 O(1)，故固定阈值即能稳定识别波区。
  `score > RefineAt` 分裂；同胞对两侧都 `score ≤ CoarsenAt` 时合并。
* **守恒插值**：
  * 分裂直接在**守恒变量**上做受限斜率 `U_L = U − ½dU`、`U_R = U + ½dU`，
    因分量线性，`U_L Δx_L + U_R Δx_R = U_parent Δx_parent` **精确成立**；
    任一子单元非正即回退分段常数（仍守恒）。
  * 合并用精确体积加权平均。
* 每次重分后迭代做 **2:1 尺寸平衡闭包**（只提级，需要时取消合并或追加分裂），
  保证相邻单元宽度差不超过 2 倍；网格始终无重叠、无缝隙、恰好覆盖全域。
* `RegridEvery` 控制重分频率。

## 5. 时间推进与自适应步长

* 时间积分用保强稳定性的二阶 **SSP-RK2（Heun）**；任一 RK 阶段产生非物理状态即整体回退
  并返回错误，不会污染结果。
* 每步 CFL 自适应：`dt = CFL · min_i Δx_i/(|u_i|+c_i)`，并在接近终点时缩短到恰好落在
  `TEnd`。任何配置下 `dt(|u|+c)/Δx ≤ CFL`（默认 CFL=0.8）。
* 长时间周期运行（多个声波穿越时间）保持有限、严格正、不发散（`TestLongTimeStability`）。

## 6. 可复现性

* 全程确定、串行、固定求和次序（网格恒按 `x0` 排序）；无随机数、无时间戳进入状态。
* 相同输入与配置下，两次运行的终态与快照文件**逐位一致**
  （`TestDeterministicOutput`）。
* 快照用 Go `encoding/json`（float64 最短往返十进制），保存/加载不丢任何一位。

## 7. 快照持久化与崩溃恢复

* 快照含 `format_version`、完整配置、`time`、`step` 以及全部叶单元（id/level/几何/U）。
* 写入采用**临时文件 + 原子 rename**，崩溃不会留下半截检查点。
* 从最近快照 `Resume` 后继续到 `TEnd`，其轨迹与连续运行**逐位相同**
  （`TestResumeMatchesContinuous`）。物理/数值配置全部取自快照；瞬态的 `MaxSteps`
  在新进程中清零，输出路径由本次命令行覆盖。

## 8. 构建与运行

```bash
go build ./...
go vet ./...
go test ./...

# 经典 Sod 激波管，反射边界，AMR
go run ./cmd/euler -ic sod -n0 128 -tend 0.2 -amr -maxlevel 3 \
  -boundary reflecting -flux rusanov -limiter minmod -cfl 0.8 \
  -snapshot snapshot.json -out solution.csv

# 周期边界，关闭 AMR，换用 HLLC + vanLeer
go run ./cmd/euler -boundary periodic -amr=false -flux hllc -limiter vanleer

# 从检查点恢复（继续到配置中的 TEnd）
go run ./cmd/euler -resume snapshot.json -out solution.csv
```

### 配置项（CLI 标志 ↔ `euler.Config`）

| 标志 | 字段 | 含义 |
|---|---|---|
| `-ic` | `ICName` | 命名初值（内置 `sod`） |
| `-xmin/-xmax` | `Xmin/Xmax` | 计算域 |
| `-n0` | `N0` | 初始单元数（2 的幂，≤ 2^24） |
| `-gamma` | `Gamma` | 比热比 γ |
| `-boundary` | `Boundary` | `periodic` / `reflecting` |
| `-flux` | `Flux` | `rusanov` / `hllc` |
| `-limiter` | `Limiter` | `minmod` / `vanleer` |
| `-cfl` | `CFL` | CFL 数 (0,1] |
| `-tend` | `TEnd` | 终止物理时间 |
| `-maxsteps` | `MaxSteps` | 步数硬上限（0 不限；恢复时清零） |
| `-amr` | `Refine.Enabled` | 启用自适应 |
| `-minlevel/-maxlevel` | `Refine.Min/MaxLevel` | 相对细化层级 |
| `-refine-at/-coarsen-at` | `Refine.RefineAt/CoarsenAt` | 分裂/合并阈值 |
| `-refine-criterion` | `Refine.Criterion` | `density-gradient`/`pressure-gradient` |
| `-regrid-every` | `RegridEvery` | 每 N 步重分（0 关闭） |
| `-snapshot-every` | `SnapshotEvery` | 每 N 步检查点（0 仅终点） |
| `-snapshot` | `SnapshotPath` | 检查点路径 |
| `-resume` | — | 从检查点恢复 |
| `-out` | `OutputPath` | 终态 CSV（`x,dx,level,ρ,u,p,ρu,E`） |

作为库使用时，可用 `euler.RegisterProblem` 注册自定义 `InitialCondition`（按 x 给原始
变量；初值用两点 Gauss 积分成守恒均值），再调用 `NewSimulation/Advance/Run/SaveSnapshot/Resume`。

## 9. 测试覆盖

`go test ./...` 共 21 个用例，覆盖需求 1–8：

* **变量转换/正性**：双向往返、负 ρ/负 p 拒绝、推进中负状态拒绝。
* **激波管**：对 Toro 精确 Riemann 解的 L1 误差、网格加密单调收敛、无超界极值。
* **边界**：周期三量机器守恒、反射壁质量/能量守恒与镜像、鬼元。
* **AMR**：分裂/合并单元级守恒、2:1 平衡、铺砌完整、波区确实加密且远场保持粗、
  周期 AMR 三量机器守恒、AMR 解精度。
* **时间步/稳定**：CFL 逐单元不违规、长时间正性稳定、dt 随细化缩小且恒正。
* **崩溃恢复/可复现**：快照逐位往返、中断后恢复轨迹逐位相同、同配置输出与快照文件逐位一致。
