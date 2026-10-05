# study_title 修正紀錄（2026-10-05）

## 問題

GSE121842、GSE171011、GSE136569 的 `study_title` 結尾被附加了「breast cancer and adjacent normal tissues」。
這段字出現在報告 Samples 區塊的頂部。

## 來源：未能從紀錄中確認

- repo 與 server 的程式（`scripts/web_ui.py`、`generate_report.py`、`preview_runner*.py`）都沒有產生這段字。
  web_ui 只寫入表單的單一字串（`save_project_snapshot` 用 `yaml.safe_dump`），不會產生 YAML 的縮排延續行。
- GEO esummary 的標題（`_fetch_geo_title` 同一公式）也不含這段字：GSE171011、GSE136569 的 GEO 標題是乾淨的，
  GSE121842 的 GEO 標題也不含 breast cancer。
- 受影響的快照在 YAML 中是多行純文字（第二行縮排的 `  breast cancer and adjacent normal tissues`），
  這是手動編輯的特徵，不是程式寫出的格式。
- 這段字與 GSE113230、GSE58135 的 breast cancer 標題用字相同，最可能是在編輯時從別處貼上。

結論：無法指出是哪一次操作。沒有證據指向任何程式或批次腳本。

## 修正

| 資料集 | 修正前 | 修正後（GEO esummary，空白正規化） |
|---|---|---|
| GSE121842 | Expression profiles of circular RNAs in human colorectal cancer based on RNA deep sequencing (Wang et al. 2019) breast cancer and adjacent normal tissues | Next Generation Sequencing Analysis of Colorectal Cancer Tissue Transcriptomes |
| GSE171011 | Next Generation Sequencing of Papillary Thyroid Cancer tissue sample and adjacent normal tissue breast cancer and adjacent normal tissues | Next Generation Sequencing of Papillary Thyroid Cancer tissue sample and adjacent normal tissue |
| GSE136569 | circNFIB suppresses lymphatic metastasis of pancreatic cancer breast cancer and adjacent normal tissues | circNFIB suppresses lymphatic metastasis of pancreatic cancer |

- GSE121842 的舊值是一個完整的論文標題，只有結尾的尾巴有問題。依指示改為 GEO 標題。若你偏好保留原論文標題、只去掉尾巴，請告訴我。
- 快照以 `yaml.safe_dump` 改寫，已確認除 `study_title` 外所有鍵值逐一相同。原檔備份在 server 的
  `~/title_fix_backup_20261005/`。

## 報告

- 只重新產生三份 `report.html`（`generate_report`，不跑 rank）。
- sha256 前後比對：三個結果目錄共 687 個檔案，只有這三份 `report.html` 變動。
- 三份報告中舊尾巴已不存在，新標題已出現。

## SRP156355 與其他

- SRP156355：快照與報告都已有標題「RNA-seq of early-stage invasive ductal carcinoma (IDC) of breast, paired
  tumor and adjacent normal」。它不是「未設定」，因此未修改，也未重新產生報告。GEO gds 查不到此 SRP（`None`），
  本機 RunInfo 也沒有 study 標題欄位，所以無從核對這個字串是否為完整的 SRA study 標題。
- 「未設定」這個字樣來自附錄 A2 產生器對空白標題的顯示。實際上空白標題的是 **GSE221107 與 GSE229705**
  （快照的 study_title 為空），而不是 SRP156355。這兩組的報告目前沒有標題；本次依指示未處理，請決定是否要補上。
