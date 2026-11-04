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

- Circle_finder 固定到 `3eb333db2ea6277dde36cbf640be9afeb710c717`；原始 MIT 源码和许可证在 `third_party/`。
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

容量单位为十进制 GB。每个样本的 `peak_gb` 应包含下载、解压后的 FASTQ/FASTA、所有比对副本、排序临时文件、caller 和 evidence。Circle_finder 原版保留多份 SAM/BAM，初测不要只按压缩 FASTQ 的两倍估算。

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

一个归档清理失败会标为 CLEANUP_FAILED，其他样本继续。DOWNLOAD_FAILED / COMPUTE_FAILED / CLEANUP_FAILED 不占 active-sample 名额，但**仍保留完整磁盘 reservation**。因此成功样本能持续循环；失败数据积压或持久 archive 满时必须处理，不能保证无人维护地无限运行。空闲轮询会打印状态统计，便于发现等待/积压。

上述恢复依赖两个控制进程继续运行。DTN/提交器进程被杀、机器重启或全局锁残留仍需恢复；不要在未确认旧进程停止时自动删除锁。

## 结果与清理

```bash
python stream.py report --config config.local.json > benchmark.tsv
python stream.py cleanup --config config.local.json --sample ont_example
```

`cleanup` 要求 ARCHIVED、归档整体 SHA256 正确、每个内部文件校验正确、scratch ownership 匹配、无 worker lock、无 symlink。只删除该样本本流水线创建的 scratch 目录。`stage-local` 的源文件和参考原件不会被删除。

归档保存：

- Circle_finder 的所有 split/discordant read IDs 对应的完整原始 FASTQ（含 mates、quality）、这些模板在完整 mapping BAM 中的所有 records、split/disc 文本和 microDNA-JT。全库 concordant 文本及 read-ID 频次表不进入归档。
- CReSIL 所有报告 graph（含 acyclic）的支持 read IDs 对应的完整原始 FASTQ；重新对全 reference 比对，保存 CIGAR/SA/secondary records；eccDNA_final、subGraphs.summary、graph 的 reference_regions.fa / final_reads.fa、仅支持 reads 的原始 trim 行（trim.support.txt）。重新比对的 BAM **不是原始 trim 内部 mappy alignment**。全库 trim 和 coverage/tmp 不进入归档。
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
