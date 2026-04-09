# Git 开源发布说明

本文档说明如何将本项目从本地目录整理为可公开仓库，并发布到 GitHub 或 Gitee。

## 1. 本仓库已经做好的开源准备

- 已补充 `.gitignore`，避免提交运行日志、缓存、原始下载数据和结果图表。
- 已补充 `LICENSE`，默认使用 MIT License。
- 项目代码、配置模板、环境定义和使用文档可以直接公开。

## 2. 发布前建议检查

公开前先确认以下内容不应进入仓库：

- `data/raw/cached/` 下的数据库导出文件
- `data/raw/geo_datasets/` 下的 GEO 原始数据
- `data/processed/` 下的中间结果
- `results/` 下的图表、表格、对接结果
- `logs/` 下的运行日志
- 任何账号、令牌、API key、本地绝对路径

如果你后续在 `config.yaml` 中写入了真实密钥，不要提交。

## 3. 本地 Git 初始化

在项目根目录执行：

```bash
git init -b main
git add .
git commit -m "chore: open source initial release"
```

如果 Git 版本不支持 `-b main`：

```bash
git init
git branch -M main
git add .
git commit -m "chore: open source initial release"
```

## 4. 发布到 GitHub

1. 在 GitHub 新建一个空仓库。
2. 仓库建议命名：
   `network-pharmacology-pipeline`
3. 不要勾选自动生成 README、`.gitignore`、LICENSE。
4. 在本地绑定远程并推送：

```bash
git remote add origin https://github.com/<your-name>/<repo-name>.git
git push -u origin main
```

如果使用 SSH：

```bash
git remote add origin git@github.com:<your-name>/<repo-name>.git
git push -u origin main
```

## 5. 发布到 Gitee

1. 在 Gitee 新建一个空仓库。
2. 复制仓库地址。
3. 在本地执行：

```bash
git remote add origin https://gitee.com/<your-name>/<repo-name>.git
git push -u origin main
```

如果仓库已经绑定了 GitHub 远程，也可以改成：

```bash
git remote rename origin github
git remote add origin https://gitee.com/<your-name>/<repo-name>.git
git push -u origin main
```

## 6. 首次公开后的建议动作

- 在仓库首页补充主题标签，例如 `bioinformatics`、`snakemake`、`network-pharmacology`
- 在 Release 中注明示例药物和示例疾病是模板案例，不代表最终研究结论
- 在 README 首页说明哪些步骤需要手工下载数据
- 如果准备长期维护，建议再补充 `CONTRIBUTING.md`

## 7. 推荐仓库可见性策略

- 如果只是课程、作业、方法展示，可直接公开
- 如果论文尚未投稿，建议先私有仓库内部整理，投稿前再公开
- 如果包含尚未发表的真实研究结果，建议只公开代码模板，不公开真实原始数据和结果
