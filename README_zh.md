# 这个包里是什么（中文说明，给自己看的）

对应投稿稿件：*Feedback quantity can flip the absorbability verdict for upstream steps in canal pools*
（投 Journal of Hydroinformatics，IWA）。

## 内容

| 目录 | 内容 |
|---|---|
| `bm01/` | 合成渠池算例（反例主体）：1D 圣维南求解器、0D 蓄量模型、控制器、SWMM 复核、`raw/` 里 99 个存档结果、`figures/` |
| `bm01r/` | 真实渠段（唐河—放水河，25.6 km）：文献参数化重跑，含质量守恒与 0D/1D 核对 |
| `README.md` | 英文主说明，含"稿件哪一节 ↔ 哪个脚本 ↔ 哪个结果文件"对照表 |
| `RUN.md` | 每一条复现命令 + 应当跑出来的数字（别人据此自查） |
| `PROVENANCE.md` | 打包时的环境版本、复跑一致性验证记录、随机性说明 |
| `checksums_sha256.txt` | 全部文件 SHA-256，防篡改 |
| `LICENSE` / `LICENSE-DATA` | 代码 MIT；数据与图 CC BY 4.0 |
| `CITATION.cff` / `.zenodo.json` | 引用信息与 Zenodo 元数据，上传时自动读取 |

## 打包时做过的验证（可以写进 response to reviewers）

- `bm01_sve1d_v2.py` 在包内重跑 → 与存档 `raw/bm01_sve1d_v2.json` **逐字节相同，只有 `runtime_s` 不同**
  （129.56 s → 140.59 s）。所有峰值、判定码、网格档、守恒指标完全一致。
- `bm01_band_sensitivity.py` 重跑 → 输出与存档**逐字节相同**。
- `bm01_ext_swmm.py` 默认档（ncell=20）重跑 → 正常结束，池均峰值 71.55 mm，与存档一致。

## 发布的两条路（都要你本人点一下，我这边没有你的账号授权）

### 路线 A：GitHub 公开库（3 分钟，立刻有 URL）

1. 浏览器打开 https://github.com/new
2. Repository name 填 `canal-pool-bm01`，选 **Public**，**不要**勾 Add a README / .gitignore / license
3. 点 Create repository，把新页面上的 SSH 地址（`git@github.com:bmberzhang/canal-pool-bm01.git`）发我
4. 我这边 `git init` + 提交 + 推送，然后校对网页

### 路线 B：GitHub 公开库 → Zenodo 拿 DOI（推荐，10 分钟）

在路线 A 的基础上：

1. 打开 https://zenodo.org，用 **GitHub 账号**登录
2. 右上头像 → Settings → **GitHub** → 找到 `canal-pool-bm01`，把开关拨到 **ON**
3. 回 GitHub，在仓库页右侧 **Releases → Create a new release**，Tag 填 `v1.0.0`，Publish release
4. 等 1–2 分钟，回 Zenodo → 头像 → **My dashboard**，会出现一条带 **DOI** 的记录（形如 `10.5281/zenodo.XXXXXXX`）
5. 把 DOI 发我，我把它和 URL 写进稿件的 Data and code availability，重出 docx + PDF

> 也可以不走 GitHub：直接 https://zenodo.org/uploads/new 拖整个文件夹上传，`.zenodo.json` 会被自动读取，
> 然后 Publish。缺点是没有 git 版本历史。

## 稿件那句话我会改成（等 URL/DOI 定下来再落笔）

> The 1D solver, the 0D model, the SWMM driver and input deck, all raw JSON/CSV results and the figure
> scripts are archived at **<URL>** (release v1.0.0; DOI **<DOI>**) under the MIT licence (code) and
> CC BY 4.0 (data). The main-case parameters are listed in Supplementary Table S20, and the reproduction
> commands, with the expected numbers for each run, are given in Supplementary Material note S11.

注意：现稿的 Highlights 第 5 条写了 "open scripts and data is provided"，摘要也写了
"open scripts and data"，**所以不能退回 "available on reasonable request"** —— 那句话和"open"自相矛盾，
审稿人一眼能看出。必须真发出去。
