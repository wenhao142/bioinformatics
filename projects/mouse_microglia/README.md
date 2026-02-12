# Mouse microglia (female, AD models) workflow

Goal: keep小鼠流程與現有人類分析分開，後續下載、處理與結果全部放在 `projects/mouse_microglia/`。

## 資料/來源（待下載）
- GSE148405：5xFAD vs WT，FACS microglia snRNA-seq，含性別/年齡（30 週）。用以 female 子集 AD vs WT DE。
- （可選）其他含 sex 標註的 microglia sc/snRNA-seq；若有 OVX/高齡模型可再補。
- 人類 AD GWAS 匯總：Kunkle 2019 (ieu-b-2) 或 Bellenguez 2022；用於 ESR1/ESR2 區域與 DE 基因交集、簡易 fine-mapping。

## 建議步驟
1) 下載 expression + metadata 到 `data/`  
2) 預處理 & 女性 microglia 子集 → DE：輸出 `results/de_microglia.tsv`、上/下調清單與簡易圖。  
3) 雌激素簽名富集（MSigDB HALLMARK_ESTROGEN_RESPONSE）。  
4) GWAS 區域擷取（ESR1 chr6:146,809,107-157,129,619; ESR2 chr14:59,226,707-69,338,613, GRCh37/38 擇一），與 DE 基因做交集；視需要做 clump/fine-mapping。

## 目前狀態
- 只建立了目錄結構，尚未下載任何檔案，避免和其他分析混在一起。

## 執行環境
- 可重用主目錄的 conda 環境 `/Users/yuewenhao/miniforge3/envs/ad-finemap`（含 Python、snakemake、samtools、bcftools、matplotlib 等）。若要 R/Seurat，另行安裝。

## 待確認
- 網路是否能順利下載 GEO 補充檔與 MSigDB GMT。若不通需換網/代理。
- OpenGWAS/EBI 是否可下載 AD GWAS；若需要 token 或鏡像，請先提供。 

## 已完成的流程與結果（2026-01-27）
- 下載：`data/GSE148405_counts.csv.gz`、`data/GSE148405_series_matrix.txt.gz`
- 套件：額外安裝 `scanpy`, `anndata`, `gseapy`
- 差異分析：5xFAD vs WT 微膠細胞（30 週，性別未標；全體視為「中年」代理）
  - 腳本：`scripts/run_gse148405_proxy.py`
  - 產物：`results/gse148405_proxy/DE_5xFAD_vs_WT.tsv`、`up_genes.txt`、`down_genes.txt`、`volcano.png`
  - 結果：上調 71 基因、下調 406 基因；ESR1 表達無顯著變化，ESR2 未檢出
- 雌激素簽名 GSEA（Hallmark）
  - 腳本：`scripts/run_estrogen_prerank.py`
  - 產物：`results/gse148405_proxy/gsea_estrogen/` 下的報告、圖、`estrogen_only.tsv`
  - 結果：Estrogen Response Early NES = -1.04 (FDR 0.77)，Late NES = -0.78 (FDR 1.0) → 未達顯著
- KEGG 富集（上調基因 ORA，Enrichr KEGG_2019_Mouse）
  - 腳本：`scripts/run_kegg_enrich.py`
  - 產物：`results/gse148405_proxy/kegg/kegg_enrichr_raw.tsv`（無 FDR<0.05；名目最小 p 為 Long-term potentiation）

## 可重用的模組
- `scripts/lib/enrich_utils.py`：包含 `load_gene_list`, `run_enrichr`, `run_prerank` 等 helper，可在其他專案直接匯入。
- `scripts/run_go_enrich.py`、`scripts/run_kegg_enrich.py` 已改為使用上述模組，支援 `--fdr`、`--topn` 等參數方便跨專案重跑。

## Snakemake 小流程（自動化 DE + GSEA/GO/KEGG）
- 檔案：`Snakefile`（此專案根目錄）
- 目標：DE、Hallmark Estrogen GSEA、GO BP (up/down)、KEGG (up)
- 主要指令：
  ```bash
  cd projects/mouse_microglia
  snakemake --cores 4
  ```
- 依賴：環境 `/Users/yuewenhao/miniforge3/envs/ad-finemap`；輸入 `data/GSE148405_counts.csv.gz`, `data/GSE148405_series_matrix.txt.gz`
- 產出：`results/gse148405_proxy/` 下的 DE、GSEA、GO、KEGG 檔案（與先前手動一致）

## 可改進／後續建議
- 取得細胞層級性別／年齡／絕經模型的元資料，再做 pre/post 子集 DE
- 放寬 DE 門檻或使用 prerank GSEA for KEGG，增加檢出率；可用上調+下調合併的 ORA
- 做 TF regulon 活性（DoRothEA/decoupler）檢查 ESR1/ESR2 或共調控因子
- ssGSEA/AUCell 在單細胞層級算雌激素簽名，避免只看偽 bulk
- 交叉人類 AD GWAS（ESR1/ESR2 區域）與小鼠 DE 同源基因；需要可用的 GWAS 檔案或 token
- 若有真實雌激素干預模型（OVX、高齡雌性、E2 補充），以該資料取代目前「中年代理」
