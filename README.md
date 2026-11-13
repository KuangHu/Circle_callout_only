# Circle_callout_only / circle-stream

在 LRC 上持续执行 **下载 → circularization candidate calling → 证据归档 → 校验 → 清理 → 下一批**。
支持 Nanopore 和 Illumina paired-end，**不以 insertion call 为入口**。

这是 **v0.1 集群试跑版**：控制层可以本地测试，生物信息学软件、环境安装、真实数据速度和 LRC 节点权限必须在集群验收。不是已经在 LRC 跑过的 production pipeline。

## 第一版范围

| 输入 | Caller | 前提 |
|---|---|---|
| Nanopore 普通 WGS | CReSIL `trim` + `identify_wgls` | 对应 genome/reference FASTA |
| Illumina 普通 WGS，paired-end，所有 reads ≥75 bp | Circle_finder bwa-mem-samblaster | 对应 genome/reference FASTA |

- 目前每个 manifest entry = 一个 ONT FASTQ 或一个有序 R1/R2 pair。不同 runs 分开登记；需要合并同一样本的多个 runs 时，先另外准备合并文件和校验值。
- 下载使用带大小和 MD5 的公开 HTTPS FASTQ URL，推荐 ENA。SRA accession 本身不能代替下载链接。v0.1 不包含 SRA → FASTQ 转换。
- 必须显式提供未压缩 reference。v0.1 **不自动组装、不选择参考、不发现参考外序列**。含 IS 的 reference 可以来自独立样本组装，不需要 insertion catalog。
- RCA/enriched/ATAC、PacBio、单端 Illumina、低于 75 bp 的普通 WGS 不混入本配置。
- 已做 adapter trimming 的 Illumina：默认 `caller.short_read_policy="reject"` 拒绝含 <75 bp reads 的样本；设为 `"drop_pair"` 时任一 mate <75 bp 则整对丢弃，其余 reads 保持原长度，`qc.json` 记录总/丢弃/保留 pairs 和保留碱基。FASTQ 截断、配对 ID、质量长度等错误仍然报错。过滤会去掉部分短片段信号，结果需注明该策略。
- 数据源：ENA HTTPS 在 LRC 实测每连接约 0.7–2.5 MB/s 且整体被限速；NCBI SRA 的 AWS 公开桶（`https://sra-pub-run-odp.s3.amazonaws.com/sra/<run>/<run>`）实测单连接 12–18 MB/s、8 连接 41–55 MB/s。manifest 中可用 `<run>.sra` 单文件条目（`size_bytes` + `s3_etag`，由 `scripts/sra_manifest.py` 通过 HEAD 生成）；下载后校验 S3 multipart ETag，preflight 再跑 `vdb-validate`，qc 用 `fasterq-dump --split-3`（临时文件在节点本地 /tmp）转 FASTQ，Illumina 无配对的 reads 丢弃并记录。需要在 control env 安装 `sra-tools`。以比对形式提交的 run（如 BAM 上传）依赖外部参考序列（refseq）；fasterq-dump 找不到时会联网下载，在计算节点上可能无限挂起。因此下载作业用 `align-info` 列出依赖，把每个 refseq 从 NCBI 下载到 `.sra` 同目录（`vdb-validate` 校验），并记入 download.json；preflight 确认全部依赖已在本地，qc 在断网（`unshare -rn`，节点支持时）和 `sra_conversion_timeout_minutes`（默认 60）限制下转换，超时即失败而不是占满作业时限。
- `download_mode`：`"dtn"`（默认：DTN 上 `download --watch` + 登录节点 `submit --watch`）或 `"slurm"`（只需在登录节点用 tmux 保持一个 `submit --watch`；它自动提交 1 核下载作业和计算作业，**不再需要 DTN 的 `download --watch`**）。slurm 模式下：
  - 下载：最多 `fetch.max_jobs` 个 1 核作业并发；每个作业逐个下载，样本按空闲 fetch 槽位均分，且不超过 `fetch.max_samples` 个、`fetch.batch_gb` GB。实际供给速度需按批次实测，不能只用 `max_jobs × 单流速度` 推算。
  - 暂存：`scratch_dir` 只放下载的源文件，预算按源文件字节计入 `scratch_budget_gb`。Lustre 元数据延迟高时（实测每次操作 4–30 s），请把暂存目录放在 VAST 组目录。
  - 计算：在节点本地 `/tmp` 完成，归档并校验后删除本地目录（失败时也删除）。`auto_cleanup=true` 时作业会顺带释放暂存的源文件；为 `false` 时保留，需显式执行 `cleanup`。
- `reference_store`（可选，强烈建议设置）：init 时把每个参考基因组按 SHA-256 命名、只读地永久保存一份。归档只记录路径和 SHA-256，不再内嵌 FASTA；路径和哈希只能确认身份，不能恢复丢失的文件，所以参考库必须永久保留，并与状态目录、日志、归档、暂存一起计入空间预算。
- 检出的对象是 **circularization-compatible candidates**。tandem duplication 等结构可能有同样 junction；caller 输出不是物理闭环证明。

## 架构

```text
DTN: lrc-xfer              共享 scratch             Slurm compute
download --watch   ->   raw / prepared / work  ->   Snakemake sample job
                              ^                    qc → reference → mapping
                              |                    → call → evidence → archive
                              |                                      |
                              +---- 归档验证后清理、释放 reservation --+

持久 state: spec、状态、job ID、日志、每条 rule 的测量
持久 archive: 每个样本的校验归档包
```

下载器、提交器是两个独立进程，通过共享文件系统协调。**不要从 compute job SSH 到 DTN**。
DTN 每次只下载一个样本；计算可并行，默认最多两个作业。字节预算和 active-sample 上限阻止下载无限堆积。
这里不假设 LRC 有 `xfer` Slurm partition，也不假设 DTN 能运行 `sbatch`。

## 安装

```bash
git clone https://github.com/KuangHu/Circle_callout_only.git
cd Circle_callout_only
bash scripts/install-tools.sh
conda activate circle-control
python -m unittest discover -s tests -v
python scripts/check_dag.py
```

在管理员允许联网安装软件的节点执行安装。安装脚本创建三个独立 conda environments；已存在的环境不会自动覆盖。使用其他环境名时同步修改 config 的 command prefixes。

- ONT 主线（v3，`scripts/ont_junctions.py`）：全基因组独立扫描 circularization junction，不依赖 CReSIL、insertion calling 或 Illumina 位点。同一 read 相邻两段比对在同一 contig、同一链上“回跳”（到达 X 后从 Y<X 继续）即一次 junction 事件，候选单元 [Y, X)；不要求 read 覆盖整个环，也不要求绕圈两次。比对：`minimap2 -c -x map-ont --secondary=yes -N 5 -p 0.9`，次级比对只用于标记重复/多重比对歧义。
  - 每个事件记录两侧锚定长度、mapq、是否唯一（mapq≥20 且无 ≥0.9×得分的次级比对）、伸出单元外的长度（同一比对从 junction 连到外侧序列）、两段重叠长度；全部写入 `ont_junction_events.tsv`（原始坐标、归属）。长 read 穿过 tandem 时 minimap2 常把多出的一份写成 CIGAR 插入：唯一比对上长度≈X−Y、位置在 [Y, X] 内、并伸出 ≥100 bp 外侧的插入也计为 tandem 分子（`insertion_tandem_molecules`）。
  - 聚类按滚动中位端点（容差 50 bp）；再用参考序列把微同源等价断点左对齐后合并，`raw_clusters` 保留原始簇及坐标。独立分子：同一 read 的所有事件先合并为一个分子（多次绕环的 read 只算一次）；再把参考上起止位置（该 read 比对链的两端，所有圈共用）相差 ≤10 bp 的 read（如 duplex 正反链、重复记录）合并；阈值都用 `independent_molecules`，且它永远不超过 `junction_reads`。限制：比读长短的小复制子/小环上，read 两端未比对上的残圈使比对链两端常落在单元边界（0, L），不同 read 会被当作重复合并，所以这里 `independent_molecules` 是保守下限（1 TB 试跑 200 样本：≥3 条不同 read 但 <3 个分子的 junction，internal 34/273，whole_contig 107/531）；同时报告 `junction_reads`。
  - 标签（证据类别，不是拓扑证明）：`tandem_supported`（≥2 个独立分子把 junction 连到唯一外侧序列；不排除同一样本同时有环）、`junction_supported`（≥1 个分子两侧锚定均唯一且 ≥300 bp）、`ambiguous`（只有重复/多重比对或短锚定）。`circle_compatible` 只是一列：没有任何分子连到外侧、且按本文库读长若为 tandem 应看到 ≥3 个时为 yes——“兼容”，不是 confirmed circle。单条越界 read 不判 tandem；同一 read 多次跨 junction（`multipass_molecules`）只记录，不作拓扑证据。
  - 位置：`whole_contig`（contig 首尾闭合，±150 bp，吸收组装器在环状 contig 两端留下的重叠；完整质粒/染色体，单独统计作技术对照）、`contig_edge`、`internal`；均标记，不删除。主结果 `support_ge3`（≥3 个独立分子），同时在 calls.json 报告 ≥1/2/3/5/10 的计数。
  - CReSIL linkage 并行运行、结果分开保存；`scripts/cresil_trace.py` 在 CReSIL 自己的 trim.txt 上复现 linkage 的区域步骤，记录每个候选在哪一步丢失：trim（mapq≥30）→ 区域长度 ≤ -mwindow（100 kb）→ 单条 read 覆盖 ≥90% 区域 → 深度/最小长度 → 进入最终图 → 成环。模拟对照与 DRR393831：3–5 kb 质粒成环保留；~60 kb 质粒在“单 read 覆盖 90%”一步丢失；>100 kb 复制子和染色体在区域长度一步丢失；三拷贝 IS 内的环在 trim（mapq）丢失。
  - 存档：候选表、事件表、支持 read ID、支持 read 的 PAF（含次级比对）和 BAM、CReSIL 图→read ID 表（`graph_reads.tsv.gz`，所有图的全部 read ID）、`cresil_trace.tsv`、`ont_library.json`（N50、比对率、覆盖度）和 `contig_depth.tsv`（同一培养核对：某质粒 contig 无 read 即提示非同一分离株/DNA 批次）。全库 PAF 不存档。
  - 支持 read 序列（v3.3）：junction counter 的全部 read 都保留；CReSIL 成环图每个最多保留 `caller.cresil_graph_read_cap`（默认 50）条 read 的序列（优先 junction read，其余按 read ID 哈希固定选取；`cresil/graph_read_sequences.json` 记录每图 read 数与保留数）。1 TB 试跑中最大的归档（1–2.7 GB，超过输入）99% 是高拷贝质粒成环图的 read。ONT 不再另存 `support_reads1.fastq.gz`：`support.bam` 每条 read 恰有一条 primary（或 unmapped）记录、含完整序列和质量值，`samtools fastq -F 0x900 support.bam` 可原样恢复；删除前核对 primary 记录数等于支持 read 数。
  - CReSIL 失败（如深文库中 worker 超内存，1 TB 试跑 1,440 例中 3 例）不再丢弃整个样本：junction counter 结果照常归档，`calls.json` 记录 `cresil.status = failed` 和错误尾部；`compare_platforms.py` 单独列出这些样本。
  - 对照：`scripts/make_ont_controls.py`（带逐 read 真值：区分“测序证据不足”和“有证据但漏检”）+ `scripts/eval_ont_controls.py`。跨平台：`scripts/compare_platforms.py`（ONT 为主，Illumina 仅交叉支持；允许断点偏移和微同源等价）。参数在 `caller.ont_junction`（可选）。
- Circle_finder 固定到 `3eb333db2ea6277dde36cbf640be9afeb710c717`；原始 MIT 源码和许可证在 `third_party/`。上游脚本不修改；适配层（`circlestream/rules.py`）只把 5 行 `grep -w -Ff ids table` 换成 `scripts/wordgrep.py`，按 GNU grep `-w -F`（C locale）语义逐字节等价地筛选，但避开 GNU grep 在大量互为前缀的 read ID 下的慢路径（单步曾达 2 小时）。空模式或含 tab 的模式交回 grep 本身处理。设置 `CS_WORDGREP_VERIFY=1` 会同时运行 GNU grep 并在任何差异时失败（仅用于验证）。
- CReSIL 固定到 `78356affa6adca6a41ea297c3d7d93f5ecb5ae6a`，沿用其环境文件，安装时 `pip --no-deps`。
- conda 求解后的完整环境清单写入 `vendor/*.explicit.txt`。安装需要网络，compute jobs 不联网安装软件。

## 配置

```bash
cp config.example.json config.local.json
cp samples.example.json samples.local.json
# 编辑这两个 JSON；示例中的 CHANGE_ME / REPLACE 故意不能运行。
python stream.py init --config config.local.json --manifest samples.local.json
```

三个根目录必须互不嵌套：

- `state_dir`：持久、小型控制和日志目录，所有相关节点可见。
- `scratch_dir`：只用于本流水线的工作目录。
- `archive_dir`：持久结果目录；不得放在会被清理的 scratch 子目录中。路径独立不等于介质已备份，请选真正持久的位置。

`slurm.account/partition/qos` 按自己的 LRC association 填写。路径不要有空格或 shell 特殊字符，因为上游 CReSIL 内部仍会构建 shell 命令。

LRC 初测注意：`lr_debug` 上限为 3 小时，示例现在用 2 小时。`lr6` 是 exclusive 节点；申请 8 核不等于按 8 核消耗资源。请以 `scontrol show partition lr6`、`sinfo -p lr6 -N -o '%N %c %m'` 和自己的 association 为准，选择有权限的共享分区，或按整节点预算做测试。增加线程数不保证 caller 等比例加速。项目账号的实际收费/配额方式需按账号核对，不能一概视为现金账单。

作业可显式加载环境：在 `slurm` 中设置 `"setup_commands": ["source /etc/profile", "module load miniforge3"]`（模块名先在本站确认）。这些是按顺序执行的 Bash 命令，失败会停止作业；配置会随 spec 冻结。默认空列表仍依赖提交时的环境。`init` 须在 `circle-control` 中执行，使保存的 Python 路径指向安装了 Snakemake 的环境。

容量单位为十进制 GB。每个样本的 `peak_gb` 应包含下载、解压后的 FASTQ/FASTA、所有比对副本、排序临时文件、caller 和 evidence。Circle_finder 原版保留多份 SAM/BAM，初测不要只按压缩 FASTQ 的两倍估算。`scripts/sra_manifest.py` 按 SRA 对象大小估算：Illumina 30×，Nanopore 30×（1 TB 试跑 1,039 个完成样本的节点本地峰值/对象大小：中位 6.2×、p99 18×、最大 20.6×；20× 时 1,440 例中 10 例在 call 步超出，最高约 23×）。slurm 模式下 `peak_gb` 只用于单个 job 的节点本地限额，不占共享 staging 预算。

`scratch_budget_gb` 是本流水线的预算，不会修改文件系统 quota。执行 rule 时每 5 秒采样一次目录大小，超过样本 reservation 会终止该进程组并保留现场；**它不是绝对的硬磁盘限额**，瞬时峰值可能漏采。archive 可用空间另行检查。共享文件系统其他用户造成的空间变化仍可能导致任务失败。

配置、workflow 和参考内容会被记录。变更线程数、阈值、caller 模式或软件版本时，用独立的 state/scratch/archive 根目录做新的 benchmark，避免旧 `.ok` 被误复用。不要在正在运行的 checkout 中 `git pull`。

## 集群上先测哪些 rule

完整验收顺序和预期在 [docs/CLUSTER_BENCHMARK.md](docs/CLUSTER_BENCHMARK.md)。建议先 synthetic controls，再各平台的小/中/大真实样本。

在 **DTN**：

```bash
conda activate circle-control
python stream.py download --config config.local.json --sample ont_example
```

在允许提交 Slurm 的节点：

```bash
conda activate circle-control
# 首先只跑到 QC；submit 自身不执行计算。
python stream.py submit --config config.local.json --sample ont_example --through qc
python stream.py status --config config.local.json
```

等待任务结束后，从已完成的 rules 接着跑：

```bash
python stream.py recover --config config.local.json --sample ont_example
python stream.py submit --config config.local.json --sample ont_example --through reference
# 也可以 through mapping / call / evidence / archive。
```

`recover` 会通过 `sacct` 确认上一个 job 已到终态；作业仍在运行或 accounting 不可用时不会重提。

只看计划、不执行任何 rule：

```bash
snakemake --snakefile workflow/Snakefile \
  --directory /YOUR_SCRATCH/circle-stream/ont_example/work \
  --cores 8 --dry-run --printshellcmds \
  --config spec=/YOUR_STATE/circle-stream/ont_example/spec.json -- archive.ok
```

## 持续运行

试跑配置默认 `auto_cleanup=false`，不会自动删除原始下载。验收通过后，创建独立 production config/run roots，设置 `auto_cleanup=true` 再初始化。

DTN 上启动下载端：

```bash
python stream.py download --config production.json --watch
```

在允许轻量 job orchestration 的节点启动提交/清理端：

```bash
python stream.py submit --config production.json --watch
```

这两个命令都在前台运行；请使用 **LRC 允许的** session/service 托管方式。`tmux` 只能解决会话断线，不会在节点重启后自动启动。v0.1 不安装 cron/systemd，也不承诺站点允许常驻进程。可由站点认可的服务定期执行单次命令；跨主机目录锁防止重复领取。

待处理清单有限，处理完所有样本后退出。可以用同一配置再次 `init` 添加新的 sample IDs，已登记样本必须保持不变。下载校验最多自动重试配置次数；计算失败保留现场，需显式恢复，防止无限消耗 SU。

每次提交前，提交器用 `sacct` 查已知 job ID，并用 `squeue` 确认它不再运行/排队。被杀而残留 RUNNING 的样本会进入 COMPUTE_FAILED、释放计算名额；不会删除其 worker lock，也不会自动重跑。查不到终态、查询失败、或 SUBMITTING 没有 job ID 时仍保留名额，避免重复计算。

一个归档清理失败会标为 CLEANUP_FAILED，其他样本继续。DOWNLOAD_FAILED / COMPUTE_FAILED / CLEANUP_FAILED 不占 active-sample 名额，但**默认仍保留完整磁盘 reservation**。1 TB 试跑中 32 个失败样本占住 43 GB 原始数据，staging 预算耗尽后准入停止。大规模运行（slurm 下载模式）可设 `"release_failed_raw": true`：COMPUTE_FAILED 样本的 Slurm job 离开队列后，删除其自有 staging 目录（ownership/symlink 校验同 cleanup），状态保留 `raw_released: true` 和原错误；`recover` 将其改回 NEW，由 fetch job 重新下载（ETag/MD5 校验）后重跑。持久 archive 满时仍须处理。空闲轮询会打印状态统计；只剩 CLEANED 和需人工处理的样本时，`submit --watch` 打印汇总后退出，不再无限轮询。

大规模运行建议把 `scratch_dir`（只放原始下载）放在 `/global/scratch/users/$USER`（Lustre，容量大，非备份），archive/state 留在组目录。

上述恢复依赖两个控制进程继续运行。DTN/提交器进程被杀、机器重启或全局锁残留仍需恢复；不要在未确认旧进程停止时自动删除锁。

## 结果与清理

```bash
python stream.py report --config config.local.json > benchmark.tsv
python stream.py cleanup --config config.local.json --sample ont_example
```

`cleanup` 要求 ARCHIVED、归档整体 SHA256 正确、每个内部文件校验正确、scratch ownership 匹配、无 worker lock、无 symlink。只删除该样本本流水线创建的 scratch 目录。`stage-local` 的源文件和参考原件不会被删除。

归档保存：

- Circle_finder 的所有 split/discordant read IDs 对应的完整原始 FASTQ（含 mates、quality）、这些模板在完整 mapping BAM 中的所有 records、split/disc 文本和 microDNA-JT。全库 concordant 文本及 read-ID 频次表不进入归档。
- ONT：junction counter 全部支持 read，以及 CReSIL 成环图（每图至多 `cresil_graph_read_cap` 条）的 read，重新对全 reference 比对为 `support.bam`（CIGAR/SA/secondary；primary 记录含完整序列与质量值，可恢复 FASTQ）；所有图的 read ID 表 `graph_reads.tsv.gz`；eccDNA_final、subGraphs.summary、graph 的 reference_regions.fa、仅支持 reads 的原始 trim 行（trim.support.txt）。重新比对的 BAM **不是原始 trim 内部 mappy alignment**。全库 trim、final_reads.fa 和 coverage/tmp 不进入归档。
- 完整 reference FASTA/FAI、代码快照、参数、MD5、版本、日志、已完成的 rule measurements。
- 零候选也保存有明确状态的结果。空 trim/不匹配参考/未识别软件错误不自动当作阴性。

不会保留完整原始文库或所有背景 BAM。未来若要更换参考/算法全面重分析，需要按 accession/URL 重新下载。当前 evidence 不完成 tandem/linear alternative 的自动排除；需要后续复核。CReSIL graph-derived sequence 不能冒称新 assembly consensus。

`archive` 自身的测量、最终状态和 Slurm stdout 留在持久 state 目录；归档里是打包前的状态快照。保守证据策略可能仍然占较多空间，benchmark 必须测结果保留量。
每个 evidence 目录的 `caller-retention.json` 列出复制的 caller 文件和字节数、未原样复制的字节数、保留的 trim 行数；`evidence.json` 记录总保留量。支持 read 集合很大时，归档仍可能很大。

## 恢复

- 下载中断：确认原下载器已停止。查看 `.download-lock/owner.json`；仅当它确实是 stale lock 时手动移走，然后在 DTN 执行 `recover --sample ID` 续传。
- `SUBMITTING` + 没有 job ID：`sbatch` 可能已被接受。按 `status.json` 中的 job_name 查询 `squeue/sacct`，补齐正确 `job.json` 后再恢复；**不要直接改回 NEW**。
- 超时/OOM：确认 Slurm 终态。处理残留 `.worker-lock` 和 Snakemake lock 前，必须确认没有进程继续使用该目录。失败的 caller 输出会被重命名保存，不混入重跑。
- 已生成归档但状态未更新：先执行 `python stream.py verify --config ... --sample ID`，仅在验证和 job 终态确认后收敛到 ARCHIVED。
- CLEANUP_FAILED：检查 `status.json` 的原因，修复问题后执行 `verify --sample ID` 重新校验并回到 ARCHIVED；下一次 submit 轮询会再尝试清理。归档不通过校验时不会释放 scratch。
- 软件/参考版本改变：新 run roots，避免复用旧结果。需要改变 memory/时间等资源时也使用新配置；可用 `stage-local` 复用保留的下载副本。

## 验证与限制

本地 tests 检查控制流、校验、配对、恢复、归档和删除保护；它们不是 ONT/Illumina caller 的准确性验证。GitHub Actions 另外在 Linux 检查 shell 和 Snakemake DAG。

速度报告区分 compressed input bytes 和 QC 后 sequencing bases。GNU time 的 MaxRSS 是进程级统计，不等于所有并发子进程内存之和；集群最终预算用 Slurm accounting/cgroup 数据。下载时间如果是恢复或本地 copy 会被标记，不能当成完整外网下载吞吐。

## 来源

- [LRC data transfer node](https://scienceit-docs.lbl.gov/hpc/data-transfer-node/)
- [LRC Slurm](https://scienceit-docs.lbl.gov/hpc/running/slurm-overview/)
- [Circle_finder](https://github.com/pk7zuva/Circle_finder)
- [CReSIL](https://github.com/visanuwan/cresil)
- [2024 方法比较](https://www.nature.com/articles/s41467-024-53496-8)（富集模拟数据的性能不能直接外推普通微生物 WGS）

第三方代码授权与适配说明见 [THIRD_PARTY.md](THIRD_PARTY.md)。
