# ViRdRp Mock Live Demo

**Interactive Protein Structure & Similarity Visualization Platform**

> [!IMPORTANT]
> **Public Demonstration**
>
> This public version removes real research data and showcases only the website visualization framework. Research data, sequence information, protein identifiers, downloadable structures, search functions, and other sensitive dataset content have been removed or anonymized.
>
> Displayed entries such as `Sample-001` and `XXXX` are placeholders and do not represent the underlying research dataset.

> [!IMPORTANT]
> **公开演示版本**
>
> 本公开版本已去除真实研究数据，只展示网站可视化框架。研究数据、序列信息、蛋白质标识符、可下载结构、搜索功能及其他敏感数据集内容均已删除或匿名化。
>
> 页面中显示的 `Sample-001`、`XXXX` 等内容均为占位符，不代表底层真实研究数据集。

## Website

[https://richard-he-happy.github.io/ViRdRp/](https://richard-he-happy.github.io/ViRdRp/)

## Demo Features

- Original-scale MCL network with zoom, pan, reset, and navigation using anonymized node identifiers
- Original-scale structure and sequence NJ trees with anonymized leaf labels
- Original-scale structure and sequence heatmaps with anonymized axes
- A single `Sample-001` Mol* structure page backed by the demonstration structure
- Hash-based single-page navigation with no application server
- Disabled Search, Download, upload, and analysis actions

## Run Locally

This is a fully static site. No Node.js installation, package installation, or custom backend is required.

```text
python -m http.server 4173
```

Open `http://localhost:4173/#/mcl`.

The same files can be hosted by GitHub Pages, Cloudflare Pages, Nginx, or any ordinary static HTTP server. Navigation uses a hash router so refreshes work on static hosting.

## Public-Data Policy

The files under `data/` preserve the original visualization data shapes while replacing identifiers, taxonomy, sequences, scores, and metadata with anonymized placeholders or deterministic mock values. The `Sample-001` structure page uses an anonymized demonstration filename; all textual metadata and download actions remain disabled.

The Mol* bundle is distributed under its upstream Apache 2.0 terms. See the header of `vendor/pdbe-molstar/pdbe-molstar-plugin.js` for attribution.

---

# 中文版本

**交互式蛋白质结构与相似性可视化平台**

> [!IMPORTANT]
> **公开演示版本**
>
> 本公开版本已去除真实研究数据，只展示网站可视化框架。研究数据、序列信息、蛋白质标识符、可下载结构、搜索功能及其他敏感数据集内容均已删除或匿名化。
>
> 页面中显示的 `Sample-001`、`XXXX` 等内容均为占位符，不代表底层真实研究数据集。

## 网站地址

[https://richard-he-happy.github.io/ViRdRp/](https://richard-he-happy.github.io/ViRdRp/)

## 演示功能

- 使用匿名节点标识符展示原始规模的 MCL 网络，并支持缩放、平移、重置和导航
- 使用匿名叶节点标签展示原始规模的结构与序列 NJ 树
- 使用匿名坐标轴展示原始规模的结构与序列热图
- 提供一个使用演示结构的 `Sample-001` Mol* 结构详情页
- 使用哈希路由实现无需应用服务器的单页导航
- 禁用搜索、下载、上传和分析操作

## 本地运行

这是一个完全静态的网站，无需安装 Node.js、软件包或自定义后端。

```text
python -m http.server 4173
```

打开 `http://localhost:4173/#/mcl`。

相同文件也可以部署到 GitHub Pages、Cloudflare Pages、Nginx 或其他普通静态 HTTP 服务器。网站使用哈希路由，因此在静态托管环境中刷新页面也能正常工作。

## 公开数据策略

`data/` 目录中的文件保留原始可视化数据的结构，同时将标识符、分类信息、序列、分数和元数据替换为匿名占位符或确定性的模拟值。`Sample-001` 结构页面使用匿名化的演示文件名，所有文本元数据和下载操作均保持禁用。

Mol* 组件沿用其上游 Apache 2.0 条款发布，署名信息见 `vendor/pdbe-molstar/pdbe-molstar-plugin.js` 文件头。
