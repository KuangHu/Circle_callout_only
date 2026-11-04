# LRC rule 验收与效率实验

目的：先确认流程正确并测到真实成本，再开持续扫描。不要从 caller 论文的运行时间推算 LRC 吞吐。

## 0. 节点与环境

在允许的节点运行 `bash scripts/preflight-cluster.sh`。确认：

1. DTN 的实际 hostname 与 config.download_host 一致，能通过 HTTPS 下载测试文件。
2. 提交节点有 `sbatch/squeue/sacct`；自己的 account/partition/QoS 可用。
3. DTN、提交节点、compute nodes 可见相同 state/scratch/archive/repo 和 conda 路径。
4. archive/state 位于持久存储；磁盘和 inode 配额有余量。
5. DTN 长期下载进程和提交器的托管方式符合 LRC 规定。v0.1 不自动部署后台服务。

## 1. 不下载真实数据的功能对照

```bash
conda activate circle-control
python scripts/make_controls.py --out /YOUR_SCRATCH/circle-functional-controls
# 编辑 config.local.json，先保留 auto_cleanup=false。
python stream.py init --config config.local.json \
  --manifest /YOUR_SCRATCH/circle-functional-controls/samples.json

python stream.py stage-local --config config.local.json --sample nanopore_circle \
  --inputs /YOUR_SCRATCH/circle-functional-controls/nanopore_circle/reads1.fastq.gz
python stream.py stage-local --config config.local.json --sample illumina_circle \
  --inputs /YOUR_SCRATCH/circle-functional-controls/illumina_circle/reads1.fastq.gz \
           /YOUR_SCRATCH/circle-functional-controls/illumina_circle/reads2.fastq.gz
```

类似方式分别 stage `nanopore_linear`、`nanopore_tandem`、`illumina_linear`、`illumina_tandem`。
同一时间只能 stage 不超过 max_active_samples 个；每个 controls entry 预留 4 GB。不要把 controls 中 `example.invalid` URL 用于 download，它们只支持 stage-local。

对照是 30 kb 随机参考、2 kb 元件，error-free reads，ONT-like 长度 1000 bp，Illumina PE150。

| 对照 | 预期解释 |
|---|---|
| circle | 应能测试跨 12000→10000 的 junction；先检查 mapping，再看 caller 是否达到阈值 |
| linear | 应无真实 circle；如果产生候选，检查其结构和阈值 |
| tandem | **线性 tandem 可能产生与 circle 同样的候选**。它测试报告是否保持歧义，不要求 caller 自动区分 |

这不是测序误差/真实菌群的 accuracy benchmark，也不能据此宣称 precision/recall。低丰度、真实 ONT indels、多拷贝 IS、混合菌株和参考缺失需要独立后续测试。

## 2. 逐条 rule

```bash
python stream.py submit --config config.local.json --sample nanopore_circle --through qc
# 等 Slurm COMPLETED，status 显示 PAUSED 后：
python stream.py recover --config config.local.json --sample nanopore_circle
python stream.py submit --config config.local.json --sample nanopore_circle --through reference
```

继续使用 `--through mapping`、`call`、`evidence`、`archive`。每次前一次 job 必须终态，先 recover 再 submit。Illumina 同样做一遍。正常继续会复用已完成的 `.ok`；改变参数必须开新 run。

| Rule | 检查对象 | 合格条件 | 测量 |
|---|---|---|---|
| download（DTN，非 Snakemake rule） | 单个真实 HTTPS 文件，含中断恢复 | MD5/字节数一致；未完成没有 READY；错 MD5 保留失败状态 | 首次完整下载 seconds/压缩 GB；续传单列 |
| preflight | FASTQ 校验、依赖、固定 CReSIL core source、路径 | 参考不变、所有软件可执行；缺失依赖必须失败 | 校验扫描成本 |
| qc | 完整 FASTQ、质量长度、重复 ID、R1/R2 配对 | 配对 ID 顺序一致；无截断；ONT 保留质量同时生成 FASTA | wall/CPU、碱基数、展开倍率 |
| reference | 本样本参考副本、FAI、BWA/MMI | genome ID 唯一；checksum 一致；索引成功 | 首次 index 成本，不算成纯 calling |
| mapping | Illumina BWA+SAMBLASTER/SAM→BAM；ONT CReSIL trim | BAM quickcheck 或 trim schema 正确 | wall、内存、临时 SAM/BAM 空间 |
| call | 固定版本 Circle_finder/CReSIL WGLS | schema 正确；CReSIL `ctc=False` 不计作圆形候选；零结果有明确记录 | wall/CPU、raw/去重候选数、状态 |
| evidence | read IDs、完整 FASTQ/quality/mates、alignment、reference | 每个 requested ID 都找回；BAM 可读/可索引 | 归档前保留量、二次比对成本 |
| archive | tar.gz + 内部 manifest + SHA256 | 每个成员读回校验；样本 spec 匹配 | wall、压缩比例、归档写入量 |
| cleanup（控制层） | ownership/worker lock/校验 | 并非 ARCHIVED 或归档损坏时拒绝；只删自身 scratch | 清理前后磁盘；原始 source/reference 仍在 |

CReSIL WGLS `depth` 模式会在 call 内部再次比对并计算 coverage，因此 `call` 的时间不是纯 graph 算法时间；详细 subprocess 日志可用于分解。

## 3. 故障测试（只在单独测试目录）

- 下载中途停止：保留 .part，确认没有活跃 worker 后恢复，完成后 MD5 相同。
- 将 manifest MD5 设置错误：不得进入 READY，更不得提交计算。
- 任务失败/超时：不得自动删除；未确认终态不得重复提交。
- 同时运行第二个 download/submit 实例：锁应阻止重复领取。
- 故意降低测试样本 `peak_gb`（仍高于压缩输入大小）：rule 达到预算时失败保留现场；不能以“少用空间”为由删除未归档结果。
- 在另一套测试归档中修改一个字节：cleanup 必须拒绝。
- CReSIL 输出 `no_candidates.json` 时：确认不是零 mapping、依赖错误或 traceback。wrapper 只接受固定版本三个明确无候选出口，其余失败不会被算作阴性。
- 检查 scratch 超限或 Slurm kill 后所有子进程都已停止，才能人工处理 stale locks。

## 4. 真实数据效率

Nanopore 和 Illumina 各取小/中/大数据集各两个。优先微生物普通 WGS，分别记录 isolate/metagenome、read length、参考来源与可比对比例。不把人类/RCA/ATAC 和微生物 WGS 合并统计。

先 8 核；每个平台选一个完整代表样本，在独立 run roots 用 4/8/16 核比较。不使用仅含已知 junction 的 read 子集测全流程性能。不同线程下候选差异必须检查；不只追求速度。

每次原始下载独立记录，复用 input/已缓存文件时标注，不把缓存命中当外网带宽。v0.1 每个样本独立建参考索引，避免共享索引的并发污染；若索引是主要开销，再做按 reference hash 共享缓存。

```bash
python stream.py report --config config.local.json > benchmark.tsv
sacct -X -P -j JOB_ID --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,ReqMem,ExitCode
sacct -P -j JOB_ID --format=JobID,State,Elapsed,MaxRSS,MaxDiskRead,MaxDiskWrite
```

`state/ID/qc.json` 有 reads/bases；`metrics/*.json` 有各次规则执行的 wall/cpu/rss/采样磁盘峰值，失败尝试也保留。

计算时请区分：

- 输入 GB：下载字节；Gb bases：测序碱基，不可混用。
- allocated core-hours = 分配 CPU 数 × 作业 elapsed 小时；实际 CPU time 单独报告。实际 LRC SU 按站点计费规则，不把这两者直接当账单。
- 单样本完成时间包含排队；纯计算效率剔除排队，但 end-to-end 吞吐必须计入。
- MaxRSS 不是整个并行进程树总内存；优先 Slurm/cgroup 峰值校准 memory request。
- `peak_gb` 要按实际最大峰值加余量，不能按最终目录大小估算。
- 候选计数必须去重；速度之外检查独立 reads、重复比对和 tandem 歧义。当前 pipeline 不输出已完成生物学验证的“高置信真环数”。

最终产出两个平台的：每 Tb bases 的 core-hours、每 TB 下载数据的网络时间、每样本临时空间与结果空间、成功率、失败原因，以及建议下载/计算并发。

## 5. 扩大到持续循环的门槛

1. 两平台完整链条都通过，阳性对照可解释，tandem 对照保留歧义。
2. 真实样本大小、运行时间和保留证据空间有实测。
3. 所有清理/恢复故障测试通过，archive 真正位于持久存储。
4. LRC 批准的两个 controller 托管方式已落实。
5. 先开小队列，不一次提交全部 accession；根据 reservation、失败积压和归档增长逐步放大。
