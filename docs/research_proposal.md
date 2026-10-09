# 缺血性卒中诱导的星形胶质细胞 APA 事件及其对 RNA 突触周分布的影响

版本：2026-09-21

## 摘要

公共数据分析提示，缺血性卒中后皮层星形胶质细胞出现时间相关的 3′UTR 使用变化。GSE238125 中的 Atp2a2 事件在第 3 天和第 7 天表现为远端 poly(A) 位点使用比例降低，并具有外周突起定位基因注释。现有自动化结果在负链事件的远端丢失片段定义上存在方向问题，涉及 QKI 重叠、motif、保守性和候选分级。对 Atp2a2 的人工链方向校正显示，真正的远端区间仍与另一个 QKI 结合峰重叠；这一结果需通过全量重建确认。现有 APA 结果来自标准短读长 RNA-seq 的算法推断，每个条件只有两个生物学样本，尚未确定真实剪切位点，也没有长、短 RNA 异构体的空间测量。

本研究拟检验：缺血性卒中是否诱导一个可重复的星形胶质细胞 APA 事件，并改变相应长、短 mRNA 异构体在胞体和突触周星形胶质细胞突起中的相对分布。研究以 Atp2a2 为暂定主候选，以完成链方向审计后的事件作为备选池。研究依次完成链方向校正、候选读段核查、3′RACE 位点确认、独立动物中的异构体比例测定、原位双探针成像、3′UTR 顺式片段和内源 PAS 扰动。后续机制和功能实验以每一阶段的预设结果为前提。主要空间终点为卒中状态与亚细胞区室对长异构体比例的交互效应。该设计将公共数据中的 APA 预测转化为可检验的事件级、异构体级和空间级证据。

## 1. 研究背景

选择不同的剪切和加尾位点可产生具有不同 3′端的 mRNA 异构体。3′UTR 长度变化能够改变 RNA 上的蛋白结合位点、miRNA 结合位点和其他顺式调控序列，从而影响 RNA 稳定性、翻译和亚细胞分布。标准短读长 RNA-seq 可以用于筛选 APA 事件，但对 poly(A) 位点的识别和使用比例估计存在局限；3′端测序和长读长测序在位点定义方面更可靠。[Shah et al., 2021](https://doi.org/10.1186/s13059-021-02502-z)

星形胶质细胞的细突起接近突触并参与局部离子和递质稳态。已有研究在皮层和海马星形胶质细胞的外周或突触周突起中检测到核糖体结合 RNA 和局部翻译，并发现部分 3′UTR 序列能够影响 RNA 定位或局部翻译。[Sakers et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28439016/)、[Mazaré et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32846133/)。神经元研究进一步说明，同一基因的不同 3′UTR 异构体可以在胞体和神经突起中呈现不同分布。[Tushev et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29656876/)

目前缺少的证据位于两个研究层级之间。一类研究测量卒中后的基因表达或推断 APA；另一类研究描述正常条件下的星形胶质细胞 PAPome，或利用少数 3′UTR 序列研究 reporter 的定位。二者尚未确定一个内源、卒中诱导的具体 APA 事件是否改变相应 RNA 异构体在 PAP 中的分布。

## 2. 前期数据与现有材料

### 2.1 GSE238125 的 APA 探索结果

GSE238125 包含 FACS 分选的小鼠皮层星形胶质细胞 bulk RNA-seq，条件为 sham 以及卒中后第 1、3、7、21、60 天，每个条件两个生物学样本。公开样本信息显示，卒中模型采用永久性远端大脑中动脉结扎并联合 7 分钟双侧颈总动脉阻断。[GSE238125 样本记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7658907)

现有 DaPars2 流程量化了 9,693 个事件，并按探索性统计规则得到 977 个候选事件。五个卒中时间点共用同一组 sham 样本，因此跨时间点结果反映时间轨迹和稳定性，不能作为相互独立的实验重复。

### 2.2 主候选与备选

Atp2a2 event 6953 位于 mm10 `chr5:122,453,512–122,457,153`，基因位于负链，预测近端位置为 `chr5:122,456,498`。逐样本 PDUI 为：

| 条件 | 重复 1 | 重复 2 |
|---|---:|---:|
| sham | 0.19 | 0.20 |
| day 3 | 0.08 | 0.06 |
| day 7 | 0.09 | 0.09 |
| day 21 | 0.11 | 0.13 |
| day 60 | 0.12 | 0.12 |

该事件在 day 3 和 day 7 达到现有探索性筛选标准，效应量分别为 ΔPDUI −0.125 和 −0.105。Atp2a2 属于现有 PAP 参考集合。Atp2a2 位于负链；现有 `lost_segments.bed` 将预测位点至高坐标端标记为丢失段，该方向与负链转录方向不符。按链方向人工校正后，暂定远端丢失段为 `chr5:122,453,512–122,456,498`，该区间与 QKI peak `chr5:122,454,200–122,454,250` 重叠。原报告使用的 `chr5:122,456,550–122,456,650` 位于共有或保留侧。Ndrg2 event 2072 和 Agpat3 event 216 的 APA 时间轨迹可保留为备选信息，其原有 QKI 丢失段重叠经同样校正后未得到支持。

链方向问题不改变 DaPars2 输出中的逐样本 PDUI 值。它会改变远端片段序列、QKI/其他 RBP peak 重叠、motif、保守性和事件分级。上述注释在全量链方向重建完成前均标记为待复核。

### 2.3 方法一致性边界

第二条计算路线使用 Salmon 转录本定量推算长异构体比例。在可比较的 27 个重点结果中，13 个与 DaPars2 同向，14 个反向；这 27 个结果中不包含 Atp2a2。QAPA 构建没有得到可用结果。因此 Atp2a2 目前只有 DaPars2 事件预测、关联注释和待核查的原始比对读段证据。

### 2.4 当前可用文件

`ATP2A2_IGV_CHECK_PACKAGE_20260920(1).zip` 已包含 sham×2、day 3×2、day 7×2 的 indexed BAM slices，覆盖 Atp2a2、Ndrg2 和 Agpat3 全基因及两侧区间；同时包含样本映射、GENCODE vM25 DaPars2 3′UTR BED、QKI peak、QC 文件和主要结果表。这些文件足以完成三个区域的读段核查。最终交付包中的事件表、PDUI 矩阵、注释矩阵和脚本用于重建链方向相关结果。完整 FASTQ 和全量 BAM 作为可追溯资产保留，当前阶段不安排全量重跑。

## 3. 科学问题与研究假设

### 3.1 科学问题

缺血性卒中是否诱导一个可验证的星形胶质细胞 APA 事件，并改变所得长、短 mRNA 异构体在胞体和 PAP 中的相对分布？

### 3.2 中心假设

卒中促进一个经验证的候选基因近端 poly(A) 位点使用，使保留远端 3′UTR 的长异构体比例下降；若远端 3′UTR 含有影响突起分布的顺式序列，长异构体下降将伴随该 RNA 在 PAP 中的相对分配改变。Atp2a2 是当前优先检验对象，最终主靶由链方向重建、BAM 核查和 3′RACE 共同确定。

该假设包含两个需要分别检验的效应：

1. 卒中改变长、短异构体的全细胞比例。
2. 卒中改变异构体在 PAP 与胞体之间的相对分配。

第二个效应以条件和区室的交互作用为判定依据，并同步控制总 RNA 丰度和星形胶质细胞形态。

## 4. 研究目标

### 目标一：确认卒中相关 APA 事件的读段基础、真实剪切位点和独立样本复现性

重建全部候选事件的链方向相关远端片段、RBP 重叠、motif、保守性和候选分级。完成 Atp2a2、Ndrg2 和 Agpat3 的 BAM/IGV 核查，形成主靶和备靶排序。随后在独立 sham 与卒中动物中用 3′RACE 确认近端和远端 RNA 3′端，并测定每只动物的异构体比例。

### 目标二：测定长、短异构体在胞体和 PAP 中的分布

在脑组织中使用共同序列探针和长异构体远端特异探针，结合星形胶质细胞膜、突触和血管标记，测量候选 RNA 在胞体、突触邻近突起和血管周终足区域的分布。主要比较为卒中状态与亚细胞区室对长异构体比例的交互效应。

### 目标三：检验远端 3′UTR 片段及内源 PAS 对 RNA 分布的作用

在前两个目标通过后，使用长/短 3′UTR reporter、候选元件删除和回补构建体测试顺式序列作用；随后采用 PAS 阻断或定点编辑改变内源异构体比例，并检测 PAP 分布。QKI 和 PABPN1 扰动属于条件性机制分支。

## 5. 研究设计与方法

### 5.1 目标一：候选事件确认

#### 5.1.1 链方向审计与注释重建

以 DaPars2 3′UTR 区间、预测 proximal PAS、基因链方向和事件 ID 为输入，为每个事件输出共有区和远端特异区。坐标统一为 BED 的 0-based、half-open 体系，并同时保存原始 DaPars2 坐标及其来源字段。区间规则为：

- 正链远端特异区：`[proximal PAS, UTR end)`；
- 负链远端特异区：`[UTR start, proximal PAS)`；
- 共有区为同一 3′UTR 中靠近终止密码子的一侧。

在正链和负链人工构造事件上验证规则后，重新计算 QKI peak、其他 RBP peak、motif 和保守性重叠，并重建 `p3_p4_event_annotation_matrix.tsv` 与候选卡。旧的 `lost_segments.bed`、QKI overlap、motif、保守性和 class 字段保留用于审计，不继续用于选靶。Atp2a2 人工校正结果作为单元测试：远端特异区应覆盖 `chr5:122,454,200–122,454,250` 的 QKI peak，不应将 `chr5:122,456,550–122,456,650` 计为丢失段。

#### 5.1.2 BAM 与 IGV 核查

核查对象为 sham×2、day 3×2、day 7×2 的 Atp2a2、Ndrg2 和 Agpat3 BAM slices。所有样本使用相同的显示范围、覆盖尺度和过滤条件。检查内容包括：

- 预测近端位置两侧的覆盖变化；
- 两只生物学重复的方向一致性；
- 末端外显子结构、剪接连接和重叠转录本；
- 低比对质量、重复区域和异常局部堆积；
- DaPars2 预测位置与 GENCODE、PolyASite 或其他 PAS 记录的关系；
- 链方向校正后的共有区和远端特异区的样本内覆盖比。

每个候选形成一张事件审查表和一组固定尺度截图。Atp2a2 的 day 7 为主比较；day 3 作为辅助证据，因为该时间点原文库存在较短插入片段和较高重复率。

#### 5.1.3 公共 PAP 与 reporter 参考

GSE143531 用于检查候选基因在健康小鼠海马 PAP polysome 与全星形胶质细胞 polysome 中是否可测及其方向。该数据的 GEO 条目包含三个生物学样本及其技术重复；技术重复在生物学样本内合并或纳入重复结构，不能直接作为独立 `n`。[GSE143531](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143531)

GSE330741 用于复核 SN-MPRA 的样本、fraction、barcode 和 tile 对应关系，并提取阳性对照片段和分析方法。该文库主要覆盖 Glt1/Slc1a2、Sparc、Glt1b 和 Hsbp1，不能提供 Atp2a2 的直接功能结果。[GSE330741](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330741)、[Koester et al., 2026](https://pubmed.ncbi.nlm.nih.gov/42094343/)

#### 5.1.4 独立动物和取材

湿实验尽量匹配 GSE238125 的卒中模型、取材皮层区域和卒中后时间。day 7 为预设主时间点；day 3 在 day 7 事件复现后作为时间扩展。动物按手术批次分层随机分配至 sham 和卒中组，样本编号在 3′RACE、比例测定和成像分析阶段保持盲法。独立实验单位为小鼠。

动物使用成年小鼠，性别构成、年龄范围、手术模型、取材坐标、术后纳入排除规则和死亡处理在实验开始前登记。优先复现 GSE238125 的永久性远端 MCA 结扎联合短时双侧颈总动脉阻断，并采集同侧皮层的预设缺血相关区域。若平台采用其他模型，模型差异写入方案，公共数据只用于候选发现。

事件复现先导阶段计划每组 4 只小鼠，用于确认检测成功率并估计每只动物 long/total 的方差。正式样本量依据 day 7 的 long/total 组间差和目标二的 `Δspatial` 计算，双侧 α=0.05，效能 80%，另按预登记的手术死亡率和样本失败率增加动物数。正式分析以小鼠为单位；不得以切片、视野或细胞数替代动物数。

#### 5.1.5 3′RACE 与位点测序

从独立动物的目标皮层星形胶质细胞 RNA 制备锚定 oligo(dT) cDNA。使用位于预测近端位点上游的基因特异引物与通用接头引物扩增 3′端，设置 no-RT 和 no-template 对照。产物采用 Sanger 测序或靶向 amplicon sequencing，确认基因特异序列与非模板 poly(A) 尾的接合处。[Ma and Hunt, 2015](https://pubmed.ncbi.nlm.nih.gov/25487210/)

候选位点附近的基因组富 A 序列单独检查。仅有 oligo(dT) 内部结合且缺少真实 RNA–poly(A) 接合证据的产物不计为 PAS。

#### 5.1.6 异构体比例测定

建立两个基础检测：

- 共同区域检测：测量长、短异构体总量；
- 远端特异检测：测量长异构体。

优先使用数字 PCR 或经扩增效率校准的 qPCR，并用已知拷贝数标准品评估线性范围。若 3′RACE 结果允许设计近端 cleavage/poly(A) junction 特异检测，则增加短异构体直接测量；否则将短异构体量标记为校准后的推算值。主要分子终点为每只动物的长异构体比例及其 sham–stroke 差值。

#### 5.1.7 目标一的进入标准

主靶进入空间实验需同时满足：

1. BAM 核查未发现足以解释事件的明显比对或注释伪影；
2. 3′RACE 确认至少两个具有 poly(A) 接合证据的 RNA 3′端；
3. 独立动物中长异构体比例的卒中效应与发现数据同向；
4. 结果不由单只动物或单一手术批次驱动；
5. 预先冻结的效应量和不确定性标准达到进入条件。

Atp2a2 未通过时，从链方向重建后的备选池中按事件质量、效应重复性、PAP 可测性和探针可设计性选择备靶。Ndrg2 和 Agpat3 只有在新的候选排序中达到门槛后进入湿实验。三个优先候选均未通过时，定位机制实验停止在事件确认阶段。

### 5.2 目标二：异构体空间分布

#### 5.2.1 探针与检测逻辑

采用两套 RNA 探针：

- common probe：位于近端 PAS 上游，检测全部候选 RNA；
- distal probe：位于经 3′RACE 确认的远端特异区，仅检测长异构体。

已知长、短构建体或合成标准 RNA 用于估计两套探针的检出率、交叉反应和双通道共定位效率。长异构体由 distal-positive 信号直接确定；common-only 信号只有在校准后才用于估计短异构体。若双探针共定位效率不足，改用 cleavage-junction 识别方案、BaseScope 或 padlock-based 3′端检测。

#### 5.2.2 组织标记和区室定义

星形胶质细胞使用能显示细突起的膜标记或遗传 reporter。GFAP 免疫信号用于细胞状态和较粗突起标记，不单独承担 PAP 边界分割。突触邻近区域由预设的突触前和突触后标记组合定义；血管以 CD31 或 lectin 标记，AQP4 用于辅助识别血管周终足。PAP、胞体和血管周区域的距离规则、图像分辨率和分割参数在读取组别前冻结。

常规共聚焦不能直接证明纳米尺度突触接触。研究将图像中的区域表述为“突触邻近星形胶质细胞突起”。一部分样本可用 Airyscan、结构光、扩增显微或等效高分辨率平台验证空间定义。

#### 5.2.3 主要和次要终点

主要终点：

\[
\Delta_{spatial}=\left(\frac{long}{total}\right)_{PAP}-\left(\frac{long}{total}\right)_{soma}
\]

比较 sham 与卒中组的 `Δspatial`，等价统计检验为 condition × compartment 交互效应。

次要终点包括：

- PAP 和胞体内的 long/total；
- 单位突起长度或体积中的总 RNA 和长异构体密度；
- 每个细胞的突起长度、分支、体积和突触邻近面积；
- 血管周终足中的相同指标；
- 候选基因总 RNA 丰度。

这些指标用于区分异构体选择性分配、总 RNA 下降和形态丢失。

#### 5.2.4 统计分析

小鼠为独立实验单位。每只小鼠可采集多个切片、视野和细胞，这些观测属于嵌套技术或生物学子样本。主要分析采用每只小鼠的预设汇总值，或采用包含小鼠随机截距的混合效应模型。模型包含 condition、compartment 及其交互项；若纳入多个时间点，则增加 time 及预设交互，并限定检验家族。报告效应量、置信区间、确切样本数和多重比较校正。探针批次、手术批次和成像批次在设计中平衡，并在必要时作为批次项处理。

### 5.3 目标三：顺式元件和内源 PAS 因果实验

#### 5.3.1 Reporter 构建

所有 reporter 使用相同启动子、编码区、载体骨架和 poly(A) 设计，仅改变候选 3′UTR。核心构建包括：

1. 经 3′RACE 确认的长 3′UTR；
2. 经 3′RACE 确认的短 3′UTR；
3. 长 3′UTR 删除候选定位元件；
4. 短 3′UTR 回补同一候选元件。

候选元件依据链方向校正后的 RBP peak、序列保守性、motif 和重叠 tile reporter 预筛结果确定。GSE330741 的 Glt1 或 Sparc 活性片段可作为方法阳性对照。RNA 定位采用 reporter RNA 的 smFISH 或 barcode 检测；总 RNA 和 reporter 蛋白分别测量。定位效应在控制总 RNA 丰度后单独报告。

初筛可在原代星形胶质细胞中完成，组织背景确认使用器官型切片或体内 astrocyte-specific 表达。体外培养结果不单独承担 PAP 机制结论。

#### 5.3.2 内源 PAS 操纵

对通过 reporter 实验的候选，采用 proximal PAS 阻断寡核苷酸或经验证的定点编辑策略改变内源长、短异构体比例。每次操纵先用 3′RACE 和比例测定确认 PAS 使用改变，再测定 RNA 的 PAP/胞体分布。加入 scrambled 对照和不影响 PAS 的序列对照。

因果证据要求：操纵改变内源异构体比例，并以相应方向改变 `Δspatial`；恢复 PAS 使用或回补候选元件后，空间表型随之恢复。

#### 5.3.3 QKI 和 PABPN1 条件分支

QKI 扰动仅在下列条件满足后启动：真实候选异构体存在，远端片段影响定位，候选区域存在可重复的 QKI 结合证据。检测内容包括 QKI 蛋白、候选 RNA 结合、异构体比例和空间分布。

PABPN1 扰动用于检验 PAS 选择，不作为预设上游结论。PABPN1 蛋白量、核定位、异构体比例和 rescue 同步测量。PABPN1 缺失促进近端位点使用的既往研究仅提供机制先验。[Jenal et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22502866/)

#### 5.3.4 条件性功能终点

若 Atp2a2 通过事件、空间和因果三层验证，则增加局部 SERCA2 蛋白分布和星形胶质细胞突起 Ca²⁺ 处理实验。候选方案包括经验证的 SERCA2 免疫检测，以及适用于星形胶质细胞突起的 cytosolic 或 ER Ca²⁺ reporter。功能实验比较恢复长异构体与维持短异构体条件下的局部 Ca²⁺ 动态，并保留总蛋白、细胞活性和形态作为控制。该部分的具体 reporter、刺激方式和成像平台在设备确认后冻结。

## 6. 结果解释框架

| 结果组合 | 解释范围 | 后续处理 |
|---|---|---|
| 链方向重建改变 QKI、motif、保守性或候选 class | 旧的序列层和候选分级不能继续使用 | 采用重建结果重新排序，保留 PDUI 层 |
| BAM 可疑或 3′RACE 不支持两个位点 | 公共数据中的目标事件未获得验证 | 换备靶；三个候选均失败则停止机制路线 |
| 位点和比例复现，空间交互不成立 | 存在卒中相关 APA，未发现异构体选择性 PAP 重分布 | 报告 APA 结果，停止定位机制实验 |
| 空间分布变化伴随总 RNA 和形态同步变化 | 结果可能由供给量或结构变化解释 | 增加归一化和形态匹配，限制定位主张 |
| Reporter 改变总 RNA或蛋白，但不改变空间分配 | 远端片段影响稳定性或翻译 | 转向 RNA 稳定性或翻译机制 |
| Reporter 与内源 PAS 操纵均改变空间分配 | 支持 APA 事件参与异构体空间分布调控 | 进入 QKI/PABPN1 和功能分支 |
| Atp2a2 通过空间与因果验证并改变局部 Ca²⁺ 终点 | 支持 APA—RNA 分布—局部功能链 | 形成完整机制论文主线 |

## 7. 项目贡献

本研究将分析单位从“某基因出现在 PAP RNA 集合中”收敛到“一个经测序确认的卒中相关 APA 事件及其长、短 RNA 异构体”。研究以同一候选的位点、比例、空间分布和序列扰动构成连续证据。公共 PAPome 和 SN-MPRA 数据承担候选筛选与实验设计作用，内源卒中样本承担事件和空间结论。

## 8. 可行性与主要风险

### 8.1 已具备条件

- GSE238125 的 APA 结果表、候选卡、样本 QC 和分析记录；
- sham、day 3、day 7 六个样本的候选区域 BAM/BAI；
- Atp2a2、Ndrg2、Agpat3 的预测位点和逐样本效应；
- GRCm38/mm10、GENCODE vM25 和样本—SRR 对应记录；
- 可直接采用的 3′RACE、数字 PCR、smFISH/RNAscope 和 reporter 方法路线。

### 8.2 关键风险

1. **链方向相关注释需要重建。**先冻结坐标体系和单元测试，再批量重跑区间、RBP、motif、保守性和候选 class；重建前不确定湿实验靶标。
2. **Atp2a2 长异构体比例较低。**采用数字 PCR 或靶向 3′端计数，提高少数异构体的定量稳定性。
3. **短异构体缺少独有内部序列。**主要检测 long/total，短异构体采用经校准的差值或 cleavage-junction 特异检测。
4. **脑切片中 PAP 边界难以分割。**使用膜标记、突触与血管多重标记，并在高分辨率平台验证区域定义。
5. **候选片段同时影响稳定性或翻译。**同步测定总 RNA、RNA 半衰期和蛋白，分别报告定位、稳定性与翻译效应。
6. **公共数据样本量较小。**公共结果承担候选发现；主要结论来自独立、按主要终点设计的动物实验。

## 9. 工作包与时间顺序

| 工作包 | 内容 | 进入条件 | 主要交付 |
|---|---|---|---|
| WP0A | 链方向与坐标重建 | 现有结果表、3′UTR BED、脚本齐全 | 新 lost/retained BED、RBP/motif/保守性、候选重排 |
| WP0B | BAM 和事件核查 | WP0A 完成、IGV ZIP 完整 | 候选审查表、IGV 图、主靶决定 |
| WP1 | 3′RACE 和异构体比例 | 主靶读段通过 | 真实 PAS、独立动物效应 |
| WP2 | 脑组织异构体空间成像 | WP1 通过 | condition × compartment 结果 |
| WP3 | Reporter 和内源 PAS 操纵 | WP2 支持空间变化 | 顺式元件及因果结果 |
| WP4 | QKI/PABPN1 和局部功能 | WP3 支持因果 | 调控和生物学后果 |

具体周期按动物模型排期、探针定制和成像平台确定。WP0–WP1 完成前不启动大规模 reporter、AAV 或功能成像。

## 10. 预期论文结构

- Figure 1：公共数据中的时间相关 APA 图谱、链方向审计、候选重排和主靶读段核查；
- Figure 2：3′RACE 确认的近端/远端位点和独立动物异构体比例；
- Figure 3：长异构体在胞体、PAP 和终足中的原位分布及空间交互；
- Figure 4：长/短 3′UTR、元件删除和回补 reporter；
- Figure 5：内源 PAS 操纵及条件性局部蛋白/Ca²⁺ 结果；
- Supplementary：候选备选、公共 PAP 数据、GSE330741 manifest、QC、阴性和停止结果。

## 11. 伦理、随机化与预登记

动物实验在机构动物伦理批准后启动，并按 ARRIVE 2.0 记录动物来源、年龄、性别、饲养条件、手术、镇痛、排除、死亡和 humane endpoint。随机序列由不参与结果判读的人员生成；样本编号在分子检测、成像和主要统计分析完成前保持盲态。主要终点、样本量规则、排除标准、批次处理、Gate 进入条件和分析脚本版本在正式实验前冻结。

## 12. 数据与代码管理

长期保留以下资产：原始 FASTQ 的位置与 MD5、全量 BAM 的位置、参考 FASTA/GTF 的版本与 checksum、候选 BAM slices、所有脚本、参数、软件版本、统计表和图像原始文件。当前 IGV 核查包不要求内嵌完整 FASTA；文件中必须保留 GRCm38/mm10、GENCODE vM25 的来源和校验信息。湿实验数据按“一只小鼠一个主目录”组织，技术重复、切片、视野和细胞均记录其层级关系。

## 参考文献与数据资源

1. Shah A, et al. Benchmarking sequencing methods and tools that facilitate the study of alternative polyadenylation. *Genome Biology*. 2021;22:291. [DOI: 10.1186/s13059-021-02502-z](https://doi.org/10.1186/s13059-021-02502-z)
2. Sakers K, et al. Astrocytes locally translate transcripts in their peripheral processes. *PNAS*. 2017;114:E3830–E3838. [PMID: 28439016](https://pubmed.ncbi.nlm.nih.gov/28439016/)
3. Mazaré N, et al. Local Translation in Perisynaptic Astrocytic Processes Is Specific and Changes after Fear Conditioning. *Cell Reports*. 2020;32:108076. [PMID: 32846133](https://pubmed.ncbi.nlm.nih.gov/32846133/)
4. Tushev G, et al. Alternative 3′ UTRs Modify the Localization, Regulatory Potential, Stability, and Plasticity of mRNAs in Neuronal Compartments. *Neuron*. 2018;98:495–511.e6. [PMID: 29656876](https://pubmed.ncbi.nlm.nih.gov/29656876/)
5. Jenal M, et al. The poly(A)-binding protein nuclear 1 suppresses alternative cleavage and polyadenylation sites. *Cell*. 2012;149:538–553. [PMID: 22502866](https://pubmed.ncbi.nlm.nih.gov/22502866/)
6. Ma L, Hunt AG. A 3′ RACE protocol to confirm polyadenylation sites. *Methods in Molecular Biology*. 2015;1255:135–144. [PMID: 25487210](https://pubmed.ncbi.nlm.nih.gov/25487210/)
7. Koester SK, et al. In Vivo Massively Parallel Reporter Assay Reveals Sequence Determinants of mRNA Localization in Astrocytes. *bioRxiv*. 2026. [PMID: 42094343](https://pubmed.ncbi.nlm.nih.gov/42094343/)
8. [GSE238125: Dynamic astrocytic gene expression profiles after ischemic stroke](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE238125)
9. [GSE143531: Local translation in perisynaptic astrocytic processes](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143531)
10. [GSE330741: In Vivo Massively Parallel Reporter Assay](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330741)
