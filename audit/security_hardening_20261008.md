# Web UI 資安強化報告（OWASP ASVS L1）— `scripts/web_ui.py`

日期：2026-10-08／09（雲端 session）。範圍：只改 repo 內程式碼與設定檔，**未部署、未連任何伺服器**。
測試：`pytest tests/test_web_security.py`（含後續新增測試，全套 149 passed）；同一批測試在**修改前**的 `web_ui.py`（搭配新的 `security.py`）上有 38 個失敗
（證明測試確實擋得住惡意輸入）。

## 0. 需要先知道的三件事

1. **repo 版本 ≠ 任務描述的版本。** `PROJECT_ID_RE`、`_validate_project_id` 在 repo 的 `web_ui.py` **不存在**，
   `PIPELINE_ALLOWED_EMAILS` 白名單也**沒有實作**（CLAUDE.md 有寫、程式碼沒有）。
   推測 server 上跑的是比 repo 新的版本。**這個 PR 以 repo HEAD 為基準；合併前請與 server 上的 `web_ui.py` 做 diff。**
   若 server 版較新，請把它放進 repo，我再 rebase 本 PR。
2. **沒有 shell 注入。** 所有 `subprocess` 呼叫本來就是 list 形式、無 `shell=True`、無 `os.system`
   （新增 AST 測試把這件事變成回歸保護）。實際存在的是**參數注入**（`gse_id` 進入 `agent.py --gse <值>`）與路徑穿越，已修。
3. **任務 1 結論：零憑證命中**，所以本 PR 不含任何輪換項目。

## 1. 逐項：改前會怎樣 → 改後會怎樣

| # | 問題 | 改前 | 改後 | 測試 |
|---|------|------|------|------|
| 1 | 命令／參數注入、accession 驗證 | `/run_gse`、`/run_manual`、`/run_local`、`/update`、`/api/detect_labels` 對 `gse_id`／`project_id`／`srr_id` 無驗證；值如 `--cores`、`GSE1;id` 會進入 agent.py 命令列、metadata 目錄名、config 檔名 | 白名單：`^(GSE\|SRP\|PRJNA)\d+\Z`（run_gse、detect_labels）、`^SRR\d+\Z`（每個 SRR，表單與 CSV）、`^[A-Z][A-Z0-9_-]{0,31}\Z`（手動／本地的自訂專案 id，不得以 `-` 開頭）；不合格 → 400 固定文字，**不回顯輸入**，且不入佇列 | `test_run_gse_rejects_*`、`test_run_manual_rejects_*`、`test_no_shell_string_subprocess_calls`（AST：任何 `shell=True`／字串 argv 即失敗） |
| 2 | 路徑穿越 | `project_id` 直接拼進 `metadata/{id}`、`config/projects/{id}.yaml`；`run_local` 的樣本 `name` 可寫出 raw_dir 外並建 symlink；symlink 目標用 `str.startswith`（`/home3/x` 放行 `/home3/xevil`）；`/api/scan_fastq` 可列舉任意目錄 | `safe_join()`（resolve 後必須在根目錄內；拒 `..`、絕對路徑、NUL、逃逸 symlink）；`safe_target()` 以路徑**元件**比對並禁止 `/etc /proc /sys /root /usr …`；symlink 建在 resolve 後路徑；`scan_fastq` 只允許 `raw_dir` 的父目錄與 `PIPELINE_FASTQ_ROOTS`，其餘 403 且不回顯 | `test_safe_*`、`test_run_local_blocks_*`、`test_scan_fastq_*` |
| 3 | WSGI 與反向代理 | 只有 `app.run()`；`ProxyFix(x_host=1)` 信任 `X-Forwarded-Host` | `deploy/gunicorn.conf.py`、`deploy/circdex-web.service`、`deploy/nginx-circdex.conf`＋`circdex_proxy.inc`（80→443、HSTS、/login 限速、請求大小上限、覆寫 XFF）。`ProxyFix` 僅在 `PIPELINE_TRUSTED_PROXY=1` 啟用，且**永不信任 X-Forwarded-Host**。新增 `start_background_workers()`，因 gunicorn 不會執行 `__main__`（否則 queue worker 與 DB 初始化都不會啟動） | `test_forwarded_host_not_trusted_*`；部署檔需在 server 實測 |
| 3b | **Host header 汙染（計畫外新增）** | magic link 以 `request.url_root` 組成 → 攻擊者 POST `/login` 帶偽造 Host，受害者收到指向攻擊者網域的登入連結（token 外洩＝帳號接管） | 新增 `PIPELINE_PUBLIC_URL`；設定後所有 email 連結一律使用它。未設定時退回舊行為並於啟動印警告 | `test_magic_link_ignores_spoofed_host` |
| 4 | Job token 強度 | `generate_job_id` 用 `random.choices` 取 4 碼（約 170 萬種）；它是 `/report/<id>`、`/download/<id>` 唯一的存取依據 | `secrets.token_urlsafe(9)`（72 bits）；舊 job id 仍可查詢。magic-link token 原本就是 `token_urlsafe(32)`（OK）；確認 log 不含 token（含 console fallback，只印主旨）；`Referrer-Policy: no-referrer` | `test_job_id_strength_*`、`test_magic_link_token_is_strong_and_not_logged` |
| 5 | CSRF 與安全標頭 | CSRF 已由 `CSRFProtect` 全域啟用（保留）；無任何安全標頭；cookie 未設 Secure | `X-Content-Type-Options`、`X-Frame-Options: DENY`、`Referrer-Policy`、`Permissions-Policy`、HTTPS 時 HSTS、CSP（見下）；`SESSION_COOKIE_SECURE`：public URL 為 https 時自動開啟（可用 `PIPELINE_COOKIE_SECURE=0/1` 覆寫；**不對純 HTTP 內網強制開啟**，否則 Secure cookie 不會被送回、無人能登入）；`MAX_CONTENT_LENGTH`=5 MB | `test_csrf_enforced_*`、`test_security_headers_present`、`test_generated_report_is_sandboxed` |
| 6 | 反射型 XSS／錯誤訊息回顯 | `status.html:272` `{{ ('"'+job_id+'"') \| safe }}`，`job_id` 來自 URL → 可跳出 JS 字串（含 404 分支）；`status_job` 404、`scan_fastq`、`serve_qc` 回顯輸入或伺服器路徑 | 改 `{{ job_id \| tojson }}`；job id 先白名單驗證，404 不帶任何輸入；所有錯誤文字固定；`serve_qc` 不再洩漏路徑；`cross_dataset.html` 內嵌 JSON 以 `_json_for_script` 轉義 `< > &` | `test_status_page_does_not_reflect_job_id`（含 `";alert(1);"`）、`test_json_for_script_*` |
| 6b | 登入白名單（計畫外確認） | 任何 email 都可取得登入連結並提交 job | 實作 `PIPELINE_ALLOWED_EMAILS`：**有設定才強制**；未設定時維持舊行為並印警告；不在名單的 email 看到與名單內完全相同的回應（無帳號枚舉）；email 格式嚴格驗證（擋 CR/LF header 注入）；`notify_email` 表單欄位同樣驗證 | `test_allowlist_*`、`test_login_rejects_malformed_*`、`test_run_manual_bad_notify_email_*` |
| 6c | `/update` 型別轉換 | `int("x")` → HTTP 500；`de_method`、`tools` 未白名單 | 數值欄位 clamp＋預設值，方法／工具白名單 | `test_update_rejects_*` |
| 7 | 相依套件漏洞（pip-audit） | — | 見 §3 | — |

## 2. CSP 說明

一般頁面：`default-src 'self'; script-src 'self' 'unsafe-inline' https://cdn.plot.ly; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; font-src 'self' https://fonts.gstatic.com data:; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'; base-uri 'self'; form-action 'self'; object-src 'none'`。

- 仍含 `'unsafe-inline'`：模板大量使用 inline script／onclick。完整 nonce 化是較大重構，風險高，超出 L1，列為後續。
  因此 CSP 目前**不能**單獨阻擋 XSS，XSS 防線靠 #6 的輸出編碼；CSP 的價值是 frame 防護、`object-src`、`base-uri`、`form-action`、外部腳本來源限制。
- `/report/*`、`/qc/*`（自包含報告，內含 inline Plotly／Highcharts）：`sandbox allow-scripts allow-popups allow-downloads` — 腳本能執行，但沒有同源權限，**碰不到登入 session**。
- 冒煙測試：以內建 Chromium 載入 `/login`、`/`、`/queue`、`/status/<id>`、`/cross_dataset`（本機測試實例）→ **0 筆 CSP 違規**。
  `/cross_dataset` 的 1 個錯誤是雲端沙箱連不到 `cdn.plot.ly`（網路），與 CSP 無關。**未驗證真實報告頁**（雲端沒有任何報告）；請在 server 上開一份 `report.html` 確認互動圖表正常。

## 3. pip-audit

- repo 原本**沒有** requirements 檔；新增 `requirements-web.txt`（Web UI 專用）。
- **Python 3.13、最新版**（Flask 3.1.3、Werkzeug 3.1.9、Flask-WTF 1.3.0、Jinja2 3.1.6、PyYAML 6.0.3、pandas、gunicorn）：**No known vulnerabilities found**。
- **Python 3.7 天花板**（CLAUDE.md 記載 server conda env `ciriquant` 為 Python 3.7.12，Werkzeug ≥3 / Flask ≥3 無法安裝）：
  Flask 2.2.5（1 個公告，修在 3.1.3）、Werkzeug 2.2.3（9 個公告，修在 2.3.8～3.1.9）**無法在 3.7 修補**；
  gunicorn、Jinja2 在 3.7 上可升到 22.0.0／3.1.6 而不再有已知漏洞。
  **這是結論性風險，不是我能在 repo 修的**：我無法連 server 確認實際安裝版本。建議：
  Web UI 改用獨立的 Python ≥3.10 環境（只裝 `requirements-web.txt` 與 gunicorn），不要放在 `ciriquant` 環境。
  這需要你決定並由有權限的人操作（依專案規則，我不動 `ciriquant` 環境）。

## 4. 部署所需環境變數（不要寫進 repo）

| 變數 | 用途 |
|------|------|
| `PIPELINE_PUBLIC_URL` | 例 `https://circdex.cgm.ntu.edu.tw`；email 連結與 Secure cookie 判斷的依據 |
| `PIPELINE_ALLOWED_EMAILS` | 逗號分隔登入白名單；**未設定＝任何 email 可登入** |
| `PIPELINE_TRUSTED_PROXY=1` | 置於 nginx 後面時設定 |
| `PIPELINE_SECRET_KEY` | 建議由 systemd EnvironmentFile 提供；未設定時沿用 `jobs/secret_key.txt` |
| `PIPELINE_FASTQ_ROOTS` | 額外允許 `/api/scan_fastq`／本地 FASTQ 的根目錄（冒號分隔） |
| `PIPELINE_COOKIE_SECURE` | `0`／`1` 覆寫自動判斷 |
| `PIPELINE_DEV_PRINT_LINK=1` | **僅本機開發**：沒有郵件服務時把 magic link 印到 console。公開部署不可開（會把登入 token 寫進 log）|

`.gitignore` 新增 `jobs/`、`.env`、`.env.*`（初版 `*.db`、`secret_key*` 過寬，已在 repo 衛生稽核中收窄）（任務 1 發現 `jobs/` 原本未被忽略；跑舊版測試時就在 repo 內產生過 `jobs/secret_key.txt`，已刪除，內容為測試用隨機值）。

## 5. 未做／殘留風險（請知悉）

- 沒有速率限制於應用層（只在 nginx 範例設定）；沒有帳號鎖定；magic link 仍可被信箱持有者轉寄（by design）。
- `'unsafe-inline'` 仍在 CSP（見 §2）。
- `gunicorn` 必須 `workers=1`（queue worker 在行程內）；建議後續把 queue worker 拆成獨立 service。
- 認證為「登入即可」，沒有每位使用者只能看自己 job 的授權檢查（`/report/<id>`、`/download/<id>` 依賴 job id 不可猜）。72-bit id 已使暴力猜測不可行，但 ASVS L1 以上建議加入擁有者檢查。
- queue 資料庫與 job registry 仍為本機檔案；未加密。
- 部署檔（nginx／systemd／gunicorn）**未在任何機器上測試過**，只是範例。

## 6. 後續補充（同一分支）

- `generate_report.py`：`study_title`、condition／label、sample id、patient id、project id 以 `html.escape` 輸出（`study_title` 來自 GEO 外部文字，原本直接插入 HTML）。範圍僅限報告的 Samples 區塊與標題；報告其餘區塊（基因名稱等來自 TSV 的欄位）**未逐一稽核**，建議另案處理。
- `PIPELINE_DEV_PRINT_LINK`（見 §4）。
