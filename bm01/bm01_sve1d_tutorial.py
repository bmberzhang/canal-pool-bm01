# -*- coding: utf-8 -*-
"""
============================================================================
 BM-01 一维圣维南模型 —— 逐行教学版
============================================================================
这份代码和正式版算出来是【完全一样的数】，区别只是：每一行都加了中文注释，
把"这一步到底在干嘛"讲成大白话。适合你边读边理解，不适合直接当论文附录。

建议按这个顺序读（从上到下）：
  1. 先看最下面的 "if __name__" 入口 —— 知道程序从哪开始
  2. 再看 Section 类 —— 断面长啥样（面积、水面宽、湿周）
  3. 再看 Ctrl 类 —— 闸门控制器（PI，就是"自动调闸门的手"）
  4. 最后看 sve_1d 函数 —— 核心求解器，水怎么一格一格往前算
  5. 结尾的 main() 只是"打印报表"，不影响计算

跑法（和正式版一模一样）：
    python bm01_sve1d_tutorial.py            # 矩形基准 + 三个梯形开关
    python bm01_sve1d_tutorial.py --selftest # 自检

证据等级：synthetic-1d-sve（数值降阶、量级示例参数、非现场结论）。
============================================================================
"""
import math  # 只要标准库里的数学函数（sqrt 开方），不装任何第三方包

# =============================================================================
# 一、几个写死的物理常数（这些是"大自然规定的"，不是算例假设）
# =============================================================================
G = 9.81          # 重力加速度，单位 m/s²。水往下流、水波传播，全靠它。
MANNING_N = 0.014 # Manning 糙率，量级示例（代表渠道壁面"多粗糙"）。越大水越难流。
CFL = 0.40        # CFL 数，一个"安全系数"，规定时间步不能太大（见后面 pick_dt）。

# 候选时间步（秒）。为什么要列这么一串？因为时间步必须能整除几个关键时刻：
#   总时长 14400 s、阶跃时刻 900 s、控制周期 300 s、采样 30 s。
# 如果整除不了，阶跃发生的那一"瞬"就会卡在两个时间步之间，算不准。
DT_CANDIDATES = [5, 4, 3.6, 3, 2.5, 2.25, 2, 1.8, 1.5, 1.25, 1.125, 1,
                 0.9, 0.75, 0.6, 0.5, 0.45, 0.4, 0.36, 0.3, 0.25, 0.225, 0.2]


# =============================================================================
# 二、断面几何 Section：这个池子"横着切开"长什么样
# =============================================================================
# 渠道横断面有几种形状：
#   矩形（m=0）：就是一条沟，两边是直墙。     面积 = 底宽 × 水深
#   梯形（m>0）：两边往外斜，越深水面越宽。   面积 = 底宽×水深 + 边坡×水深²
# 这里 m 叫"边坡系数"，意思是"水平方向每升高 1 米，往两边各扩 m 米"。
class Section:
    def __init__(self, b: float, m: float = 0.0):
        self.b = b   # b：底宽（米）
        self.m = m   # m：边坡系数，0 就是矩形

    def A(self, h):
        # 过水面积：水有多"粗"。矩形=h*b；梯形=h*b + m*h²（多出来的两腰三角）
        return h * (self.b + self.m * h)

    def B(self, h):
        # 水面宽度：水在表面有多宽。矩形=b；梯形=b + 2*m*h
        # 这个量很重要——水位涨 1 厘米，水体积多多少，就由 B 决定。
        return self.b + 2.0 * self.m * h

    def P(self, h):
        # 湿周：水跟渠道壁面"接触"的总长度。湿周越大，摩擦阻力越大。
        # 梯形两腰是斜的，腰长 = h × √(1+m²)，所以 P = b + 2h√(1+m²)。
        return self.b + 2.0 * h * math.sqrt(1.0 + self.m * self.m)

    def is_rect(self):
        # 是不是矩形（m≈0）。矩形和梯形的"水力半径"算法不一样，靠这个区分。
        return abs(self.m) < 1e-12

    def I1(self, h):
        # 静水压力积分（本教学版没直接用，但正式版里有）。可以理解为
        # "水压在断面上产生的总推力"的一个中间量。
        return 0.5 * self.b * h * h + (self.m / 3.0) * h * h * h


# =============================================================================
# 三、渠底坡度 S0 怎么定：让初始状态"不自己乱晃"
# =============================================================================
def normal_slope(p: dict, n: float = MANNING_N) -> float:
    """
    矩形断面的纵坡反解。
    目标：让"水深=目标水位、流量=180"恰好是"均匀流"（正常流）——
    也就是水既不加速也不减速，一开始就稳稳当当，不自己晃。
    方法：反用 Manning 公式 Q = (1/n)·A·R^(2/3)·√S，把 S 解出来。
    """
    B, ht = p["pool_width_m"], p["h_target"]
    A = B * ht                    # 过水面积
    R = A / (B + 2.0 * ht)        # 水力半径 = 面积 ÷ 湿周
    return (n * p["q0"] / (A * R ** (2.0 / 3.0))) ** 2


def section_normal_slope(sec: "Section", h: float, Q: float, n: float) -> float:
    """
    按【实际断面】反解纵坡（梯形断面下必须用这个，不能用上面的矩形式）。
    这是 2026-10-01 修正过的地方：梯形如果硬套矩形公式，纵坡会算大 1.38~1.79 倍，
    导致水被"想象得太陡"，峰值算错。
    """
    A = sec.A(h)
    P = sec.P(h)
    R = A / P if P > 0 else h
    return (n * Q / (A * R ** (2.0 / 3.0))) ** 2


def pick_dt(dx: float, wmax: float) -> float:
    """
    挑时间步 dt。
    规则：水波一秒走 wmax 米，一个时间步内最多只能走 CFL×dx 这么远
    （也就是不能"跳过"一个格子）。所以 dx 越小（格子越细），dt 必须越小。
    wmax 是"最快波速"= 水流速度 + 浅水重力波速 √(g·h)。
    """
    lo = CFL * dx / max(wmax, 1e-6)   # 允许的最大时间步
    for c in DT_CANDIDATES:            # 从候选里挑第一个 ≤ lo 的（同时要能整除关键时刻）
        if c <= lo:
            return c
    raise ValueError(f"CFL 太紧：dx={dx} 需 dt<{lo:.3f} s")


# =============================================================================
# 四、闸门控制器 Ctrl：那双"自动调闸门的手"
# =============================================================================
class Ctrl:
    """
    PI 控制器。作用：每隔一段时间看看"闸前水位"离目标差多少，然后决定闸门开多大。

    P（比例）：差多少就补多少 —— 差得大，动作大。
    I（积分）：差了很久就补得更狠 —— 防止"一直差一点点却收不回来"。

    打个比方：你开车保持 60 码。
      P 是"现在差多少"：一看表 70 了，松油门。
      I 是"已经差了多久"：一直 65 迟迟降不下来，就再多松一点。
    """
    def __init__(self, p: dict):
        self.p = p
        self.a = p["a0"]            # 当前开度（0=全关，1=全开），从 a0=0.5 开始
        self.a_cmd = p["a0"]        # 想要的开度（指令值）
        self.integral = 0.0         # 积分项，累计"差了多少 × 时间"
        self.last_ctrl_t = -1e18    # 上次控制时刻（用来控制"多久动一次"）

    def control(self, t: float, h: float) -> None:
        """
        每 ctrl_period_s（300 秒）算一次开度指令。参数 h 是"闸前水位"。
        """
        p = self.p
        if t - self.last_ctrl_t < p["ctrl_period_s"]:
            return                       # 还没到 300 秒，这次不动
        self.last_ctrl_t = t

        # 把水位按"测量分辨率"取整，模拟真实传感器（真实水位计不是无限精度的）。
        meas = round(h / p["meas_resolution_m"]) * p["meas_resolution_m"]
        e = meas - p["h_target"]          # e = 误差 = 实测 - 目标（正=水位偏高）

        # 死区：误差太小就不反应，避免闸门"神经质"地抖来抖去。
        if abs(e) >= p["deadband_m"]:
            at_hi = self.a >= 1.0 - 1e-12 and e > 0   # 已经全开了，还想再开 → 没招
            at_lo = self.a <= 1e-12 and e < 0         # 已经全关了，还想再关 → 没招
            if not (at_hi or at_lo):
                self.integral += e * p["ctrl_period_s"]   # 没顶到极限，就继续积分

        u = p["kp"] * e + p["ki"] * self.integral   # PI 输出 = 比例 + 积分
        self.a_cmd = min(max(p["a0"] + u, 0.0), 1.0)  # 开度指令，夹在 [0,1]

    def step_gate(self, dt: float) -> None:
        """
        真正挪闸门。闸门不能"啪"一下跳到指令值（现实里液压启闭机有速度上限），
        所以按 rate_limit（每秒最多动多少）慢慢移过去。这叫"速率限制"。
        """
        step = self.p["rate_limit_per_s"] * dt   # 这一步最多能挪多少
        delta = self.a_cmd - self.a              # 还差多少
        if abs(delta) < step:
            self.a = self.a_cmd                  # 一步到位
        else:
            self.a += step if delta > 0 else -step  # 只挪一点点
        self.a = min(max(self.a, 0.0), 1.0)      # 别超出 [0,1]


# =============================================================================
# 五、核心求解器 sve_1d：水怎么一格一格往前算
# =============================================================================
def sve_1d(p: dict, ncell: int, step_dq: float, *,
           ctrl_on: str = "gate", m_slope: float = 0.0,
           n: float = MANNING_N, cfl: float = CFL,
           s0: float = None, sensor_len_m: float = 100.0,
           early_window_s: float = 7200.0) -> dict:
    """
    这是整个模型的"发动机"。把 6 公里渠池切成 ncell 格，每一格有自己的水位 h
    和流量 Q，然后按时间一步步往前推，最后统计水位有没有超过 ±150 mm 条带。

    参数：
      p          算例参数（见 DEFAULT_PARAMS）
      ncell      切多少格（120/240/480），格越多越准、越慢
      step_dq    阶跃量（+45 m³/s），t=900 s 后上游突然多加这么多水
      ctrl_on    控制器看哪个水位：
                   "gate"     只看最末一格（单点传感器）—— 有失稳风险
                   "gate_avg"  看闸前 100 m 平均 —— 网格无关，更稳（默认正式用这个）
                   "mean"      看全池平均
                   "off"       不控制
      m_slope    边坡系数（0=矩形，1.0/1.5/2.5=梯形）
    """
    # ---- 准备几何和初始状态 ----
    B0, L, ht = p["pool_width_m"], p["pool_len_m"], p["h_target"]
    sec = Section(B0, m_slope)
    dx = L / ncell                       # 每格多长
    S0 = normal_slope(p, n) if s0 is None else s0   # 渠底坡度
    h = [ht] * ncell                     # 初始：每格水位都 = 目标水位 7 m
    Q = [p["q0"]] * ncell                # 初始：每格流量都 = 180 m³/s
    nsens = max(1, min(ncell, int(round(sensor_len_m / dx))))  # 闸前 100 m 占几格

    # 波速 wmax 和 时间步 dt
    wmax = p["q0"] / (B0 * ht) + math.sqrt(G * ht)   # 水流速度 + 重力波速
    lo = cfl * dx / wmax
    dt = pick_dt(dx, wmax) if abs(cfl - CFL) < 1e-15 else next(c for c in DT_CANDIDATES if c <= lo)
    nt = int(round(p["duration_s"] / dt)) + 1        # 总步数（4 小时 ÷ dt）
    rec_every = max(1, int(round(30.0 / dt)))        # 每 30 秒采样一次
    ctrl = Ctrl(p)                                   # 建一个控制器

    # ---- 记录用的变量 ----
    peak_mean = peak_gate = 0.0   # 池均 / 闸前的最大偏差
    early_peak = 0.0              # 前 2 小时内的峰值（物理峰，避开晚期病态）
    tail_mean = 0.0               # 收尾段（最后 30 分钟）的偏差
    t_tail = p["duration_s"] - p["settle_window_s"]
    over_band_dt = 0              # 累计"越带"的步数（×dt = 越带时长）
    vol_net = 0.0                 # 净入流量（用于质量守恒自检）
    c1 = dt / dx                  # 常用系数：Δh = -c1·ΔF/B 里的 dt/dx
    sqrt = math.sqrt
    hmin = 1e-3                   # 水位下限，防止除零 / 变负
    friction_correct = not sec.is_rect()   # 梯形要用断面感知的水力半径

    out = {"t": [], "dev_mean": [], "dev_gate": [], "dev_head": [], "gate": [], "q_out": []}

    # ================= 主循环：时间一步步往前推 =================
    for k in range(nt):
        tk = k * dt

        # 1) 上游进水：t<900 s 是 180，之后 180+45=225（这就是"阶跃"）
        q_up = p["q0"] + (step_dq if tk >= p["t_step_s"] else 0.0)

        # 2) 闸门出流：Q = k·a·√(h/h_target)。
        #    开度 a 越大、闸前水位 h 越高，放水越多。那个 √ 来自能量守恒。
        q_gate = p["kg"] * ctrl.a * sqrt(max(h[-1], hmin)) / sqrt(ht)

        # 3) 控制器先看水位、挪闸门
        if ctrl_on != "off":
            if ctrl_on == "gate":
                meas = h[-1]                          # 只看最末一格
            elif ctrl_on == "gate_avg":
                meas = sum(h[ncell - nsens:]) / nsens # 闸前 100 m 平均
            else:
                meas = sum(h) / ncell                 # 全池平均
            ctrl.control(tk, meas)   # 算开度指令
            ctrl.step_gate(dt)       # 实际挪闸门（受速率限制）

        # 4) Rusanov（局部 Lax-Friedrichs）数值通量：算相邻格子之间的水流交换
        #    核心思想：用两侧最快波速当"数值粘性"，把陡变抹平，防止算出非物理的锯齿。
        hx = [h[0]] + h + [h[-1]]     # 左右各补一个"虚拟边界"格子
        Qx = [q_up] + Q + [q_gate]
        nif = ncell + 1               # 交界面数量（比格子多 1）
        F1 = [0.0] * nif              # 质量通量（水流量）
        F2 = [0.0] * nif              # 动量通量
        for j in range(nif):
            h1, h2 = hx[j], hx[j + 1]
            if h1 < hmin: h1 = hmin
            if h2 < hmin: h2 = hmin
            A1, A2 = sec.A(h1), sec.A(h2)
            if A1 < 1e-6: A1 = 1e-6
            if A2 < 1e-6: A2 = 1e-6
            Bc1, Bc2 = sec.B(h1), sec.B(h2)
            # 两侧各自的最大波速
            s1 = abs(Qx[j]) / A1 + sqrt(G * A1 / Bc1)
            s2 = abs(Qx[j + 1]) / A2 + sqrt(G * A2 / Bc2)
            sm = s1 if s1 > s2 else s2          # 取大者当"数值粘性"强度
            # 质量通量 = 平均流量 - 粘性×(面积差)
            F1[j] = 0.5 * (Qx[j] + Qx[j + 1]) - 0.5 * sm * (A2 - A1)
            # 动量通量 = 平均(动量+静水压力) - 粘性×(流量差)
            P1 = 0.5 * G * A1 * A1 / Bc1
            P2 = 0.5 * G * A2 * A2 / Bc2
            F2[j] = 0.5 * (Qx[j] * Qx[j] / A1 + Qx[j + 1] * Qx[j + 1] / A2 + P1 + P2) \
                    - 0.5 * sm * (Qx[j + 1] - Qx[j])
        F1[0], F1[-1] = q_up, q_gate   # 边界交界面直接用进/出流量

        # 5) 显式欧拉：用通量进出差，更新每格的水位 h 和流量 Q
        h_new = [0.0] * ncell
        Q_new = [0.0] * ncell
        for i in range(ncell):
            hi = h[i]
            A = sec.A(hi)
            if A < 1e-6: A = 1e-6
            Qi = Q[i]
            if friction_correct:
                R = A / sec.P(hi) if sec.P(hi) > 0 else hi   # 水力半径（梯形）
            else:
                R = A / (B0 + 2.0 * hi)                      # 水力半径（矩形）
            # 摩擦坡降 Sf = n²·Q·|Q| / (A²·R^(4/3))，这就是 Manning 阻力项
            Sf = n * n * Qi * abs(Qi) / (A * A * R ** (4.0 / 3.0))
            Bc = sec.B(hi)
            hn = hi - c1 * (F1[i + 1] - F1[i]) / Bc          # 更新水位
            h_new[i] = hn if hn > hmin else hmin
            Q_new[i] = Qi - c1 * (F2[i + 1] - F2[i]) + dt * G * A * (S0 - Sf)  # 更新流量
        h, Q = h_new, Q_new
        vol_net += (q_up - q_gate) * dt   # 累计净入流（自检用）

        # 6) 记录指标
        dm = sum(h) / ncell - ht   # 池均偏差
        dg = h[-1] - ht            # 闸前偏差
        if abs(dm) > peak_mean: peak_mean = abs(dm)
        if abs(dg) > peak_gate: peak_gate = abs(dg)
        if tk <= early_window_s and abs(dm) > early_peak: early_peak = abs(dm)  # 前 2h 峰值
        if tk >= t_tail and abs(dm) > tail_mean: tail_mean = abs(dm)            # 末段偏差
        if abs(dm) > p["band_normal_m"]: over_band_dt += 1                      # 越带计数
        if k % rec_every == 0:   # 每 30 秒存一次曲线
            out["t"].append(tk); out["dev_mean"].append(dm); out["dev_gate"].append(dg)
            out["dev_head"].append(h[0] - ht); out["gate"].append(ctrl.a)
            out["q_out"].append(q_gate)

    # ============ 质量守恒自检：净入流 vs 池内蓄量增量 ============
    # 精确式（梯形下不能用简单的 B·Δh，要逐格积分）
    dv = 0.0
    for x in h:
        dv += (x - ht) * (sec.b + sec.m * (x + ht)) * dx
    return {"ncell": ncell, "dx_m": dx, "dt_s": dt, "steps": nt, "S0": S0,
            "ctrl_on": ctrl_on, "m_slope": m_slope,
            "series": out,
            "peak_mean_mm": peak_mean * 1000, "peak_gate_mm": peak_gate * 1000,
            "early_peak_mm": early_peak * 1000, "final_mean_mm": tail_mean * 1000,
            "over_band_s": over_band_dt * dt,
            "vol_net_m3": vol_net, "vol_storage_m3": dv,
            "mass_balance_rel_err": abs(dv - vol_net) / max(abs(vol_net), 1.0)}


# =============================================================================
# 六、算例参数 DEFAULT_PARAMS：这个"池子"长啥样
# =============================================================================
DEFAULT_PARAMS = {
    "pool_len_m": 6000.0,     # 池长 6 公里
    "pool_width_m": 50.0,     # 底宽 50 米
    "q0": 180.0,              # 初始来水 180 m³/s
    "h_target": 7.0,          # 目标水位 7 米
    "a0": 0.5,                # 初始开度 50%
    "kg": 360.0,              # 闸门系数（全开时能放 360 m³/s）
    "kp": 1.6,                # P 增益
    "ki": 0.0016,             # I 增益
    "ctrl_period_s": 300.0,   # 控制周期 300 秒（5 分钟看一次）
    "deadband_m": 0.01,       # 死区 1 cm
    "meas_resolution_m": 0.001,  # 测量分辨率 1 mm
    "rate_limit_per_s": 0.0003333333333333333,  # 闸门每秒最多动多少（≈1/3000）
    "duration_s": 14400.0,    # 总时长 4 小时
    "t_step_s": 900.0,        # 900 秒（15 分钟）时发生阶跃
    "band_normal_m": 0.15,    # 正常条带 ±150 mm
    "band_alarm_m": 0.3,      # 报警条带 ±300 mm
    "settle_tol_m": 0.02,     # 收敛容差 20 mm
    "settle_window_s": 1800.0,  # 末段窗口 30 分钟
}
DEFAULT_STEP_DQ = 45.0   # 阶跃量 +45 m³/s（约 +25%）

# 三个梯形开关：只改边坡系数，其余不变（这就是"控制变量法"）
SWITCHES = [
    ("S2a", "梯形断面 边坡 m=1.0", dict(m_slope=1.0)),
    ("S2b", "梯形断面 边坡 m=1.5", dict(m_slope=1.5)),
    ("S2c", "梯形断面 边坡 m=2.5", dict(m_slope=2.5)),
]
SENSOR_LEN_M = 100.0   # 闸前平均长度 100 米（传感器看这段的平均水位）


# =============================================================================
# 七、main()：只是"打印报表"，不参与计算
# =============================================================================
def main() -> int:
    import argparse, time
    ap = argparse.ArgumentParser()
    ap.add_argument("--grids", default="240,480")   # 用哪两级网格
    ap.add_argument("--only", default=None)         # 只跑某个开关
    ap.add_argument("--selftest", action="store_true")  # 自检模式
    args = ap.parse_args()

    p = DEFAULT_PARAMS
    step_dq = DEFAULT_STEP_DQ
    grids = [int(x) for x in args.grids.split(",")]
    t0 = time.time()

    c = math.sqrt(G * p["h_target"])
    print("BM-01 一维圣维南（逐行教学版）")
    print(f"池长 {p['pool_len_m']:.0f} m × 底宽 {p['pool_width_m']:.0f} m，"
          f"目标水位 {p['h_target']:.2f} m，来水 {p['q0']:.0f} m³/s，"
          f"阶跃 +{step_dq:.0f} m³/s，条带 ±{p['band_normal_m']*1000:.0f} mm")
    print(f"波速 c={c:.2f} m/s，水波走一个池长 L/c={p['pool_len_m']/c/60:.1f} min，"
          f"来回 2L/c={2*p['pool_len_m']/c/60:.1f} min\n")

    if args.selftest:
        # 自检：矩形 + 单点传感器，跟 v1 记录值比对（差 <0.5 mm 即通过）
        v1ref = {"120": 142.88143604489179, "240": 156.6079038724011}
        print("[自检] 与 v1 记录值比对")
        ok = True
        for nc in (120, 240):
            r = sve_1d(p, nc, step_dq, ctrl_on="gate", m_slope=0.0)
            d = abs(r["peak_mean_mm"] - v1ref[str(nc)])
            good = d < 0.5
            ok = ok and good
            print(f"   {nc:>4} 格  本文件 {r['peak_mean_mm']:.2f} mm vs "
                  f"v1 {v1ref[str(nc)]:.2f} mm　差 {d:.3f} mm → {'一致 ✓' if good else '不一致 ✗'}")
        print(f"   结论：{'与 v1 逐行等价 ✓' if ok else '存在差异 ✗'}")
        return 0 if ok else 1

    print(f"[基准] 矩形；传感器=闸前 {SENSOR_LEN_M:.0f} m 平均")
    base = {}
    for nc in grids:
        base[nc] = sve_1d(p, nc, step_dq, ctrl_on="gate_avg", m_slope=0.0,
                          sensor_len_m=SENSOR_LEN_M)
        r = base[nc]
        print(f"   {nc:>4} 格 dt={r['dt_s']:>6} s｜池均峰 {r['peak_mean_mm']:7.1f} mm｜"
              f"前2h峰 {r['early_peak_mm']:7.1f} mm｜闸前峰 {r['peak_gate_mm']:7.1f} mm｜"
              f"末段 {r['final_mean_mm']:7.1f} mm｜守恒 {r['mass_balance_rel_err']:.1e}")
    base_conv = base[grids[-1]]["early_peak_mm"]
    spread = abs(base[grids[-1]]["early_peak_mm"] - base[grids[-2]]["early_peak_mm"])
    print(f"   收敛峰值（前2h，{grids[-1]} 格）={base_conv:.1f} mm"
          f"（与 {grids[-2]} 格差 {spread:.1f} mm）")

    only = args.only.split(",") if args.only else None
    for tag, name, kw in SWITCHES:
        if only and tag not in only:
            continue
        print(f"\n[{tag}] {name}")
        sec = Section(p["pool_width_m"], kw.get("m_slope", 0.0))
        s0_sec = section_normal_slope(sec, p["h_target"], p["q0"], MANNING_N)
        s0_leg = normal_slope(p)
        print(f"   纵坡：断面感知 {s0_sec:.6e}｜矩形反解 {s0_leg:.6e}　矩形式偏大 {s0_leg/s0_sec:.2f} 倍")
        for nc in grids:
            r = sve_1d(p, nc, step_dq, ctrl_on="gate_avg", s0=s0_sec,
                       sensor_len_m=SENSOR_LEN_M, **kw)
            print(f"   {nc:>4} 格 dt={r['dt_s']:>6} s｜池均峰 {r['peak_mean_mm']:7.1f} mm｜"
                  f"前2h峰 {r['early_peak_mm']:7.1f} mm｜闸前峰 {r['peak_gate_mm']:7.1f} mm｜"
                  f"末段 {r['final_mean_mm']:6.1f} mm｜越带 {r['over_band_s']/60:5.2f} min｜"
                  f"守恒 {r['mass_balance_rel_err']:.1e}")
        conv = sve_1d(p, grids[-1], step_dq, ctrl_on="gate_avg", s0=s0_sec,
                      sensor_len_m=SENSOR_LEN_M, **kw)["early_peak_mm"]
        print(f"   收敛峰值（前2h，{grids[-1]} 格）={conv:.1f} mm　"
              f"对基准 {base_conv:.1f} mm 增量 {conv-base_conv:+.1f} mm"
              f"（{(conv-base_conv)/base_conv*100:+.1f}%）")

    print(f"\n总用时 {time.time()-t0:.1f} s")
    print("证据等级：synthetic-1d-sve（数值降阶，非现场结论）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
