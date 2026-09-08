# Agent 使用指南:rnaseq-downstream

本文件是 AI agent(Codex、Kimi Code、DeepSeek harness 或任何能执行 shell 的
agent)操作本工具包的入口文档。执行任何命令前请完整阅读。本文说明如何驱动
工具包、每种统计方法做什么、以及现有证据支持和不支持哪些声明。

英文版见 `AGENTS.md`;两版内容应保持一致,如有出入以英文版为准并提请维护者修正。

## 这个工具包是什么

一个证据导向的 bulk RNA-seq 下游分析工具包。它刻意只实现少量分析路径,
并拒绝不支持的请求。只有通过相应闸门的路径才列为 evidence-gated;DESeq2
目前已实现但未过闸。每个 run bundle 都绑定带哈希的输入与锁定运行时;
evidence-gated 声明还会绑定已通过的机器可读 benchmark 证据。

**成熟度:research preview(研究预览)。** 本工具包的任何输出本身都不构成
"某个研究或结论已可发表"的声明。benchmark 证据可以证明同引擎数值忠实;
只有存在且通过适用的模拟闸门时,才能证明被测场景下的校准。它不认证你的研究。

## 不可协商的操作纪律

1. **永不绕过 hard fail。** 拒绝(`DESIGN_RANK_DEFICIENT`、
   `CONTRAST_NOT_ESTIMABLE`、`INPUT_EVIDENCE_REQUIRED` 等)是工具在按设计
   工作。不要修改输入去骗过闸门;修复底层问题,或把拒绝如实报告给用户。
2. **永不修改已发布的输出目录。** 证据包原子发布、永不覆盖。要改分析,
   发起一次新的 run。
3. **stdout 是一份 JSON 文档。** 解析它;不要抓 stderr(那里只有日志)。
4. **一次 run 只用一个后端。** 永不把两个引擎的结果合并成一个结论
   ("取交集当真值"被禁止)。证据门控的 edgeR QL 路径是主路径;DESeq2
   路径是研究预览级的备选。
5. **仅离线。** 基因集、注释、参考文件全部是本地冻结文件。不要用在线
   资源替代。
6. **如实报告证据等级。** 向用户总结结果时,区分 `evidence_gated` 路径与
   `implemented_ungated` 路径。

## 五个命令

入口:`rnaseq-downstream <command>` 或 `python -m rnaseq_downstream <command>`。

| 命令 | 用途 |
| --- | --- |
| `capabilities` | 报告已实现的命令、证据门控与未过闸路径、锁定运行时身份。每个 session 先调用它。 |
| `inspect` | 只读审计声明的输入:provenance、SHA-256、Salmon 元数据与 inferential replicate 检测。不做分析。 |
| `validate` | 验证输入语义,产出输入验证包(`run` 的前置条件)。 |
| `run` | 对已验证输入包执行分析请求,原子发布证据包。 |
| `summarize` | 独立验证 run 产物(哈希、schema、状态网格、BH 重算)并汇总结论。绝不盲信后端报告值。 |

响应信封(固定):`schema_version`、`command`、`status`
(`success`/`error`/`partial`)、`data`、`warnings`、`errors`、`artifacts`。
退出码:`0` 成功 · `2` 请求无效 · `3` 科学验证失败 · `4` 后端失败 ·
`5` 显式部分完成 · `70` 内部错误。

稳定错误码包括:`INVALID_REQUEST`、`INPUT_READ_FAILED`、
`INPUT_EVIDENCE_REQUIRED`、`INPUT_INTEGRITY_FAILED`、`COUNT_VALUES_INVALID`、
`SAMPLE_SET_MISMATCH`、`GENE_ID_INVALID`、`ASSAY_PROTOCOL_REQUIRED`、
`SALMON_OFFSET_REQUIRED`、`DESIGN_RANK_DEFICIENT`、`COVARIATE_CONFOUNDED`、
`CONTRAST_NOT_ESTIMABLE`、`BACKEND_FAILED`、`PARTIAL_RUN`、
`FEATURE_NOT_IMPLEMENTED`、`INTERNAL_ERROR`。把它们当作机器可读契约:
匹配 `code` 字段,不要匹配 message 文本。

## 标准工作流

```bash
rnaseq-downstream capabilities
rnaseq-downstream inspect   --request examples/input-requests/salmon-full-length.request.json
rnaseq-downstream validate  --request <input-request.json>
rnaseq-downstream run       --request <analysis-request.json>
rnaseq-downstream summarize <run-bundle-dir>
```

请求模板在 `examples/input-requests/` 和 `examples/analysis-requests/`。
精确 schema 见 `docs/contracts/input-request-v1.md`、
`docs/contracts/analysis-request-v1.md`、`docs/contracts/cli-v1.md`。

## 输入语义(三条已验证路线)

输入类型必须**显式声明**;工具包绝不从文件内容或 dtype 猜测。

| 类型 | 含义 | 关键规则 |
| --- | --- | --- |
| `featurecounts_integer` | featureCounts 整数 counts | 严格整数文法;必须带来源证据(manifest 或 per-sample 文件);裸 merged 矩阵一律拒绝。 |
| `salmon_quant_dirs_full_length` | 每样本 Salmon `quant.sf` 目录,全长转录组协议 | tximport 且强制 transcript-length offset(`countsFromAbundance="no"`)。inferential replicates:0 个允许且不传播不确定性;恰好 1 个拒绝;≥2 个一致则由 edgeR 作为 inferential overdispersion 传播。 |
| `salmon_quant_dirs_three_prime` | Salmon,3′ 标签协议 | 要求显式声明 `assay_protocol: three_prime`。用 `countsFromAbundance="no"` 的原始 counts,**不加**长度 offset(长度校正会给 3′ 数据引入偏差)。inferential replicates:0 个允许且不传播不确定性;恰好 1 个拒绝;≥2 个一致则允许并记录,但不用于长度校正或 inferential-overdispersion 调整。 |

对 DESeq2 而言,两条 Salmon 路线都记录允许的 inferential replicates,但不消费
它们。全长路线使用 `DESeqDataSetFromTximport`,其内部转换调用 R `round()`;
三撇号路线在 `DESeqDataSetFromMatrix` 前显式调用 R `round()`。两条路线均记录
取整前后矩阵哈希、变化单元数、最大绝对变化和每样本总计数变化。R 采用
IEC 60559 ties-to-even 规则,不允许截断。

内部主键一律为稳定 gene ID;gene symbol 仅用于展示。样本与 metadata 必须
精确匹配——永不静默取交集、填补或取整;上述 DESeq2 Salmon 转换是显式且
带审计的。

**信任边界:** featureCounts combined-matrix 的 manifest 是自证证据。闸门能
验证其内部一致性与摘要绑定,但无法从密码学上证明矩阵确实由 featureCounts
产生。

## 统计方法详解

### 1. 差异表达:edgeR v4 准似然(证据门控路径 `edger_ql_p0_v1`)

管线固定并记入 provenance:

```
filterByExpr(design) → normLibSizes(method="TMM")
  → glmQLFit(abundance.trend=TRUE, robust=TRUE,
             winsor.tail.p=c(0.05, 0.1), legacy=FALSE,
             top.proportion=NULL, keep.unit.mat=FALSE)
  → glmQLFTest(contrast=weights)  (lfc_threshold = 0)
  → glmTreat(contrast=weights, lfc=threshold, null="interval")
      (lfc_threshold > 0)
```

- **`filterByExpr`** 移除表达量低到无法检验的基因;它使用设计矩阵,
  在设置表达过滤条件时考虑实验的有效重复结构。被过滤的基因以
  `filtered` 状态报告,绝不静默丢弃。
- **TMM 归一化**(`normLibSizes`)估计每样本的组成因子,使计数比较不被
  少数高表达基因或文库大小差异混淆。
- **准似然 GLM**(`glmQLFit`):RNA-seq 计数相对 Poisson 过度离散;NB
  离散度逐基因估计,而 *QL 离散度*描述相对 NB 均值—方差关系的逐基因
  剩余变异,并向丰度趋势做 moderation。`robust=TRUE` 限制离散度异常基因
  的影响。这是选择稳健 QL 路线的设计依据;现有证据并未隔离证明稳健性
  就是其与 DESeq2 Wald 校准差异的因果原因。
- **`glmQLFTest`** 是对 contrast 的准似然 F 检验——默认显著性检验。
- **`glmTreat(null="interval")`** 是正式效应阈值检验:原假设为
  |log2FC| ≤ 阈值。它*不是*在 p 值之上做事后 |logFC| 筛选;当科学问题是
  "变化至少 X 倍"时用它。注意 `glmTreat` 不报告 F 统计量
  (`statistic_status: not_reported`)。
- **Design lint(拟合前,fail-closed):** 设计矩阵 QR 秩、残余自由度为正、
  完全混杂检测(`limma::nonEstimable` 别名)、contrast 可估计性(contrast
  向量在设计行空间内)。违规即以 `DESIGN_RANK_DEFICIENT` 或
  `CONTRAST_NOT_ESTIMABLE` 中止,不发布任何产物。
- **支持的设计:** 仅加性项(带显式参照水平的因子、连续协变量)。
  配对/批次设计用加性 block 表达(如 `subject + condition`)。交互项、
  splines、随机效应不在契约内。
- **Contrast** 是设计矩阵列上的任意权重向量,因此 `(A+B)/2 − C` 这类组合
  也可表达。FDR 为每个 contrast 内的 Benjamini–Hochberg,并由
  `summarize` **独立重算**。

证据:同引擎 airway oracle(归档运行中系数/logFC/PValue/F/FDR 的观测数值
差异为 0;闸门容差 rtol=1e-6)+ compcodeR 负二项模拟 gate(FDP/TPR 限值)。
范围:后端内核;输入路线由契约 + 锁定集成测试覆盖。

### 2. 基因集检验:limma fry / mroast / camera(路径 `edger_ql_p0_v1_limma_gene_sets_v1`)

与 DE 分析在同一拟合、同一设计、同一 contrast 上运行。

- **`fry`**(主)与 **`mroast`**(佐证)是 **self-contained** 检验:
  原假设为"该 set 内没有任何基因差异表达"。回答"这条通路到底有没有变",
  通过设计矩阵天然支持配对/block 设计。`fry` 是 `roast` 的快速旋转近似;
  `mroast` 对全部 set 做旋转检验,种子策略有记录。
- **`camera`** 是 **competitive** 检验:原假设为"该 set 内基因不比补集
  更差异",逐 set 估计基因间相关(`inter.gene.cor=NA`,VIF 有记录)。
  回答"这个 set 相对背景是否富集"。**camera 目前只有同引擎 oracle 证据,
  没有模拟 FDR/TPR gate,标注为 `supplementary`。**
- 两类检验的原假设不同。它们的结果、BH 族、报告章节严格分开;
  绝不要混排成一个显著性列表。
- 基因集来自冻结的本地 GMT 文件,带版本 + SHA-256 和逐 set 的 ID 映射审计
  (mapped/unmapped/ambiguous 计数、映射率、tested universe 内计数)。
  低于最小 tested 基因数的 set 以 `not_tested` 报告并附明确原因。
- **已知证据边界(已披露):** 模拟 fixture 逐基因独立生成,因此 set 内
  相关结构下的校准**不在** gate 覆盖范围内。

### 3. DESeq2 路径(`deseq2_p1_v1_gate_pending`)—— 研究预览,未过闸

由 DESeq2 1.52.0 + apeglm 实现:Wald contrast、LRT omnibus 模式(reduced
design 必须是正确嵌套的加性真子集、恰好一个 reporting contrast、
`lfc_threshold=0`、输出标注 omnibus)、经
`results(lfcThreshold=θ, altHypothesis="greaterAbs")` 的正式阈值检验,
以及 apeglm LFC 收缩(contrast 必须恰好是一个非截距设计系数上的 `+1`)。

**使用前必读:**

- 同引擎 airway oracle:数值完全一致(Wald 与 LRT)。
- compcodeR 校准 gate(6v6,NB 模拟):**未通过**——mean FDP 0.11821,
  预设限值 0.065(相同 seed 上 native edgeR QL 为 0.04834)。
- 机理诊断(已归档,假设生成级):假阳性集中于高离散 null 基因
  (最高真离散度五分位贡献约 38.9% 的 FP);FP 基因的 final/true 离散度
  比值偏低。independent filtering 与拟合退化不能解释该结果,对齐 tested
  gene family 后差距仍然存在。直接运行 DESeq2 复现了工具包结果;这排除了
  本次诊断中可观测的 wrapper/路由数值偏差,而非所有可能的实现错误。
- 解读:归档的同引擎 oracle 中,封装与官方 DESeq2 数值忠实一致;但该被测
  场景中的 DESeq2 结果反保守。可用于探索性工作与文献对照;正式结论以
  edgeR QL 路径为准。完整证据:[探索报告](tests/simulation/deseq2-compcoder-exploratory-report.json)、
  [方法审计](tests/simulation/deseq2-compcoder-method.md),以及机理诊断的
  [人读版](tests/simulation/deseq2-compcoder-mechanism-diagnostic.md)和
  [机器版](tests/simulation/deseq2-compcoder-mechanism-diagnostic.json)。

### 4. QC 数学(展示层)

- 公开 edgeR 展示 PCA 使用 R 后端导出的、经过 `filterByExpr` 与 TMM 后的
  logCPM 矩阵。它选取 top 高变基因、只居中、**不再缩放**。导出的 logCPM
  与 PCA 都仅用于展示,永不进入推断。
- 另有经过测试的 legacy 数学工具可拟合联合模型
  `Y = X_biology·β + X_nuisance·γ + ε`,只移除 nuisance 分量(等价于
  `limma::removeBatchEffect(..., design=X_biology)`)并拒绝混杂校正;但 adjusted
  PCA 不属于当前公开的 v1.1/v1.2 展示契约。
- volcano 与 MA 图从已发布结果表渲染;PCA 从单独导出的 display-only logCPM
  矩阵渲染。展示层不重算推断统计量(每个 display manifest 中
  `statistical_recalculation: false`)。

### 5. 结果语义

每个 基因 × contrast 行恰好携带一个状态:`filtered`、`not_tested`、
`not_estimable`、`failed`、`tested`。数值字段只在显式状态配合下为空;
不存在含糊的 NA。`summarize` 会在每个 family 内重算 BH、验证全部产物
摘要与五/六文件清单,验证不通过的 bundle 一律拒绝。

## 环境与 provenance

锁定运行时:双层锁(conda-lock 管 Python/R 4.6.1/工具链 + renv.lock 管
Bioconductor 3.23 包,源码归档以 SHA-256 钉死)。环境演进遵循
`environment/README.md` 的快照政策:旧环境文件保存在 append-only 快照中,
每次变更重跑全部已发布 live evidence gate 并出兼容性报告,断言历史数值
产物逐字节不变。CI 在每次 push 时从零重跑这些 live gate。未通过的 DESeq2
校准研究仅作归档;CI 阻塞式检查其披露与完整性,但不重跑或放宽该闸门。

## 目前不提供的功能

交互作用(difference-in-differences)、随机效应/纵向模型(dream)、
voomLmFit 样本加权、排名型探索富集(fgsea)、TF 活性(decoupleR)、
转录本层面分析(DTE/DTU)、WGCNA、meta 分析。不要在工具包之外临时拼凑
这些分析;通过维护者提出需求,以便走认证流程加入。

## 出错了怎么办

把稳定错误码和失败阶段原样报告给用户。不要改动输入重试,不要自己取整
counts,不要放松阈值。一次带清晰错误码的拒绝,就是契约的一次成功执行。
