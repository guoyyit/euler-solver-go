# euler1d — 一维可压缩流动求解器

一个用纯 Go（无第三方依赖）实现的本地一维可压缩欧拉方程数值求解器，面向含
**激波与稀疏波**的黎曼问题。它在机器精度内守恒质量/动量/能量，拒绝负密度与负
压力，支持周期性与反射边界、激波捕获 MUSCL 重构、守恒自适应网格细化（AMR）、
自适应时间步、确定性逐位复现，以及基于快照的崩溃恢复。

- 模块路径：`github.com/guoyyit/euler-solver-go`
- 代码包：`pkg/euler`（库）、`cmd/euler1d`（命令行）
- Go 版本：见 `go.mod`

---

## 1. 控制方程与变量

一维理想气体（量热完全气体）欧拉方程：

```
∂t U + ∂x F(U) = 0,

        ⎡ ρ   ⎤            ⎡    ρu     ⎤
  U  =  ⎢ ρu  ⎥,   F(U)  = ⎢ ρu² + p   ⎥
        ⎣ E   ⎦            ⎣ u(E + p)  ⎦
```

- 守恒变量 `Conserved{R, M, E}`：质量、动量、单位体积总能量。
- 原始变量 `Primitive{Rho, U, P}`：密度、速度、压力。
- 状态方程：`p = (γ−1)(E − ½ρu²)`，声速 `c = √(γp/ρ)`。

**双向转换**（`state.go`）：

- `EOS.PrimToCons` / `EOS.ConsToPrim`
- `ConsToPrim` 会对 `ρ ≤ 0` 返回 `ErrNegativeDensity`，对 `p ≤ 0` 返回
  `ErrNegativePressure`。任何 RK 阶段一旦得到非物理状态都会被拒绝，时间步会以
  半步重试，确保不会静默地产生负状态。

---

## 2. 空间离散与激波捕获

- **数值通量**（`flux.go`）：
  - `FluxRusanov`：局部 Lax–Friedrichs/Rusanov，单调、最稳健；
  - `FluxHLL`：Harten–Lax–van Leer（Davis 波速估计），激波更锐利。
- **MUSCL 分段线性重构**（`reconstruct.go`），在原始变量上逐分量做斜率限制：
  - `LimitNone`（无限制，二阶但非 TVD，仅供对照）
  - `LimitMinmod`（minmod，最稳健）
  - `LimitMC`（monotonized-centered，更锐利，仍 TVD）

  minmod/MC 都是 TVD 限制器：外推面值必落在三点模板范围内，配合物理状态即可
  保持正性，跨激波/稀疏波不产生超调振荡。间断位置与峰值随网格加密收敛
  （见 `L1` 自收敛测试，约一阶收敛率，符合 TVD 激波捕获方法的理论预期）。

时间格式为 **SSP-RK2（Heun，强稳定保序二阶）**：

```
U¹ = U⁰ + dt R(U⁰)
Uⁿ⁺¹ = U⁰ + ½dt(R(U⁰) + R(U¹))
```

---

## 3. 边界条件（`boundary.go`）

- `BoundaryPeriodic`：虚拟层从对侧内部精确复制，总守恒量不变。
- `BoundaryReflecting`：固壁镜像虚拟层——`ρ,p` 镜像、`u` 取反（守恒量中
  `R,E` 镜像、`M` 取反）。壁面面通量的质量分量恒为零，无穿透；边界只做镜像
  不做外推，因此不引入非物理振荡。

---

## 4. 自适应网格细化（AMR）

两级（或多级，`Ratio` 默认 2）块结构化细化，关键性质如下：

- **细化判据**（`criterion.go`）：基于密度/压力/速度的一阶归一化跳跃
  `JumpCriterion{RefineThreshold, CoarsenThreshold, Eps}`，带迟滞（细化/粗化阈
  值可分别配置）以抑制抖动。判定依据完全由可序列化参数决定。
- **缓冲区**：`MergeBufferCells` 把标记向外扩张若干格，保证波不会立刻跑出细网
  格；间距很近的组件会被桥接。
- **守恒插值**（`interp.go`）：
  - 粗→细：对守恒量做分段常数注入（细子格守恒密度=父格），逐位满足
    `Σ_k U_child·dx/r = U_parent·dx`；
  - 细→粗：体积加权平均，权重之和恰为 1。
- **粗细耦合与回流（reflux）**：SSP-RK2 的两个阶段在时间上一致地耦合——细网格
  阶段 0 的虚拟层取自时间 n 的粗网格，阶段 1 取自粗网格 n+1 预测值。记录两个 RK
  阶段在粗细界面上的通量差并进行 Berger–Oliger 式回流校正，使复合更新严格
  守恒（测试中达到 ~1e-15）。
- **唯一覆盖统计**：总质量/动量/能量对每个空间点只计一次——活动细格覆盖其粗
  父格，间隙/未覆盖区域使用粗格。
- 周期域上的细化区间被强制保持在内部（远离周期缝），使每个粗细界面都是普通
  内部面，单次回流即可精确闭合全局平衡。

---

## 5. 自适应时间步

- 每个时间步由当前状态确定 `dt = 0.9 · CFL · dx / max_i(|u_i|+c_i)`，对所有
  层级取最细间距的最小值；0.9 的安全系数保证激波发展、波速增大时仍满足 CFL。
- 最后一步缩短以精确落在目标时刻。
- `dt` 是当前状态的确定性函数：同输入/同配置下，相同状态必然产生相同 `dt`，
  因而整个步进序列可复现。
- 若某阶段出现非物理状态（极端舍入情形），以半步重试（最多 20 次）；正常配置
  下不会触发。

---

## 6. 逐位可复现性

- 数值核只用 IEEE-754 `float64`，不使用 wall-clock、map 迭代序、goroutine 调度
  或任何随机源。
- 相同输入与配置两次运行的所有单元格逐位一致（见确定性测试）。
- 允许以快照恢复历史结论。

---

## 7. 快照与崩溃恢复（`snapshot.go`）

- 版本化、小端、与架构无关的二进制格式（魔数 + `snapshotVersion`），只含标量
  与定长数组，无指针/map，故同状态产生同字节。
- 持久化内容：完整配置、当前时间与步号、每个层级的几何/间距/活动区间以及全部
  守恒状态。
- **检查点只在完整时间步结束后写入**（`StepOne` 之后）。从该快照恢复并继续，
  与从不间断的连续运行逐位一致（因为 `dt` 由恢复后的状态确定性地重算）。
- API：`solver.Snapshot()` / `euler.Restore(data)`、
  `euler.WriteSnapshot(w,s)` / `euler.ReadSnapshot(r)`，以及每步回调
  `Solver.CheckpointFunc`。

---

## 8. 命令行用法

```bash
go build -o euler1d ./cmd/euler1d

# Sod 激波管，200 个基格，2 级 AMR，跑到 t=0.2
./euler1d -problem sod -n 200 -levels 2 -t 0.2 -out prof.txt

# 周期性写自然步检查点（每 50 步），之后从某个检查点恢复到更晚时刻
./euler1d -problem sod -n 200 -levels 2 -t 0.24 \
    -snapshot cp -checkpoint-every 50
./euler1d -resume cp.150 -t 0.24 -out prof_resume.txt
```

| 参数 | 默认 | 说明 |
| --- | --- | --- |
| `-problem` | `sod` | 内置问题：`sod` / `lax` / `2raref` |
| `-n` | `200` | 基础网格格数 |
| `-t` | `0.2` | 终止时刻 |
| `-cfl` | `0.8` | CFL 数 `(0,1]` |
| `-flux` | `hll` | `rusanov` / `hll` |
| `-limiter` | `minmod` | `none` / `minmod` / `mc` |
| `-bc` | 问题默认 | `periodic` / `reflect` |
| `-levels` | `1` | AMR 层级数（1=均匀网格） |
| `-refine` | `0.05` | 细化阈值 |
| `-coarsen` | `0.005` | 粗化阈值 |
| `-snapshot` | — | 结束时写快照的路径 |
| `-checkpoint-every` | `0` | 每 N 个自然步写 `<snapshot>.<step>` |
| `-resume` | — | 从快照恢复 |
| `-out` | — | 输出基础网格上的 `x ρ u p` 剖面 |

内置问题（`problems.go`）：Sod 激波管、Lax 激波管、双稀疏波（近真空，强正性
考验，使用反射边界）。

---

## 9. 库的最小用法

```go
cfg := euler.Config{
    X0: 0, X1: 1, NBase: 200,
    Boundary: euler.BoundaryPeriodic,
    Flux:     euler.FluxHLL,
    Limiter:  euler.LimitMinmod,
    CFL:      0.8, Gamma: 1.4,
    AMR: func() euler.AMRConfig {
        c := euler.DefaultAMRConfig()
        c.MaxLevels = 2
        return c
    }(),
}
p := euler.GetProblem(euler.ProblemSod)
s, _ := euler.NewSolver(cfg, p.Init)
if err := s.Advance(0.2); err != nil { /* ... */ }
mass, momentum, energy := s.TotalConserved()
data, _ := s.Snapshot()
restored, _ := euler.Restore(data)
```

---

## 10. 自动化测试

```bash
go test ./...
go test ./... -race     # 竞态检测
```

覆盖的场景：

- **变量转换/正性**：守恒↔原始往返、负密度/负压力拒绝、声速。
- **激波管**：Sod/Lax 的机器精度守恒、密度 L1 自收敛、峰值与间断位置收敛、
  双稀疏波近真空正性。
- **边界**：周期/反射虚拟层正确性、壁面零质量通量、反射长时间稳定。
- **通量/重构**：等状态通量一致性、超声迎风、线性分布精确重构、间断处 minmod
  无超调。
- **AMR**：守恒注入/限制恒等式、Sod 全程机器精度守恒、AMR 相对纯粗网格降低
  误差、AMR 反射壁质量/能量守恒与正性。
- **崩溃恢复**：自然步检查点→恢复→继续与连续运行逐位一致；快照字节稳定、文件
  I/O 往返一致、总量保持。

---

## 11. 范围与限制

- 显式二阶 MUSCL + SSP-RK2；在激波处为一阶（TVD 限制器），光滑区为二阶。
- 当前 AMR 默认演示两级；`Ratio`、层级数、阈值、缓冲区、判据 Eps 均可配置。
- 一维、理想量热完全气体。
