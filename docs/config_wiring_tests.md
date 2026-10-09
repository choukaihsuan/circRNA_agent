# 設定參數接線測試（config → 命令列）

## 要擋的 bug

「config 選項存在、UI 看起來改了，實際沒傳到腳本」——靜默接線缺口。專案已發生三次
（`components_arg`、`sig_cap`／`fc_cap`、`tools.ciriquant.anchor`），每次都要逐行比對才發現。

## 做法

`tests/test_config_wiring.py` + `tests/wiring.py`：

1. 在暫存目錄建立一個假專案（空的 genome 佔位檔、4 個樣本的 metadata，不需要任何真實資料）。
2. 用 **Snakemake dry-run**（`snakemake -n -p --forceall`）取得「真正會被執行的命令字串」。**不會執行任何指令。**
3. 對 `consensus.*`、`de.*`、`tools.*` 底下每個純量鍵（以及 rules 以 `.get(key, default)` 讀取但
   config.yaml 沒宣告的鍵，例如 `consensus.adaptive`、`consensus.adaptive_ratio`），設成**明顯偏離預設**的值
   （整數 +7、浮點 ×1.7+0.0123、字串加 `_X`、布林取反），重新 dry-run，斷言：
   - 新值以獨立 token 出現在某個被改變的命令上；**或**
   - 該鍵出現在某個 `script:` rule 的 `params` 中（見 4）；
   - 否則測試失敗，訊息為 `SILENT WIRING GAP`。
4. `script:` rule（`analysis.R`、`isoform_switching.R`、`generate_report.py`）的 params 不會出現在 dry-run 的命令裡，
   所以另做靜態檢查：rule 傳的每個 param 都必須被腳本讀取、腳本讀的每個 param 都必須在 rule 定義
   （後者就是「mock 腳本少一個 key → NULL」那一類 bug）。
5. 選擇器類的鍵（`consensus.tools`、`consensus.adaptive`、`de.tumor_label`／`normal_label`）改成斷言 job 圖／旗標的變化。
6. 發現的新鍵會自動納入；`de.biomarker.*`、`tools.ciriquant.*` 一旦出現在 `config.yaml`，就自動被測，不需改測試。

已知落差登記在 `KNOWN_GAPS`（strict xfail：有人把它接好時測試會翻成失敗，逼人把條目刪掉，清單不會腐爛）。

### 為什麼用 dry-run 而不是更輕的做法

命令字串由 `.smk` 的模組層級 Python（`adaptive_flag`、依 `USE_CIRIQUANT` 條件的 input、lambda、`expand`、
`shell.prefix`）與 Snakemake 自己的 wildcard／params 解析共同組成。用 regex 重寫一份等於測的是重寫版，而不是管線本身；
靜默接線缺口恰好就藏在這層黏合碼裡。dry-run 約 1–2 秒、唯讀、只在暫存目錄操作，符合「不在伺服器、不碰已發表結果」的限制。

### 驗證測試本身有效（突變測試）

在本機把以下四種破壞各做一次，測試都正確失敗，之後已還原：
`--min-bsj` 寫死、刪掉 `--adaptive-ratio`、script rule 少傳 `heatmap_top_n`、script rule 多傳一個腳本不讀的 param。

## 圖檔來源測試

`tests/test_figure_provenance.py`：`audit/figs/` 底下每個 `.png`，`audit/` 內（figs 以外）必須有已 commit 的腳本
（`.py/.R/.sh/.ipynb`）在內容中出現該檔名或其主檔名。目前 repo 沒有 `audit/figs/`，該測試 skip；
檢查器本身以合成目錄樹測試（孤兒會被抓、不相干的腳本不算）。

## 執行

```bash
pip install -r requirements-test.txt
pytest -q                    # 全部（約 1 分鐘，大部分是 25 次 dry-run）
REQUIRE_SNAKEMAKE=1 pytest tests/test_config_wiring.py   # 沒裝 snakemake 時失敗而不是 skip
```

CI：`.github/workflows/tests.yml`（repo 先前**沒有**任何 CI，這是新建的；沒有「既有 CI」可接）。

## 限制與誠實說明

- **repo HEAD 沒有** `rank_biomarkers_v2`、`de.biomarker.*`、`tools.ciriquant.*`，也沒有 `components_arg`、
  `sig_cap`、`fc_cap`、`anchor`（全 repo 與全 git 歷史只在 CLAUDE.md 新增的規則裡出現）。這些很可能在尚未推上 repo 的版本中。
  所以任務所指定的三組鍵，我**無法**直接斷言；改為做成會自動涵蓋它們的通用機制，並以
  `test_expected_sections_present`（xfail）把「缺這些」明確記錄下來。把那份版本推上 repo 後，這些鍵會被自動測到。
- 測試以 Snakemake 9.x 渲染；server 為 7.24（Python 3.7）。命令字串渲染在兩個版本間一致，但若日後依賴 9.x 才有的行為，需在 server 版本上再確認。
- `script:` rule 的檢查是靜態的（名稱對應），不驗證 R／Python 腳本內部是否真的用到該值。
