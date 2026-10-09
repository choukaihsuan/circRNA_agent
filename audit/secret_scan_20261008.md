# Git 歷史機密外洩稽核（唯讀）

- 日期：2026-10-08（雲端 session）
- 範圍：`git log --all`，**51 個 commit**，分支 `main`、`claude/beautiful-tesla-f98688`（含 origin/*）；
  0 個 tag、0 個 stash、`git fsck --lost-found` 無懸空物件。歷史中曾出現過的路徑共 72 個。
- 方法：
  1. gitleaks v8.18.4（`detect --log-opts="--all" --redact`）
  2. 全歷史 patch（36,987 行）的手動 regex：SECRET_KEY／token／API key／SMTP／Bearer／
     私鑰標頭／Slack webhook／AWS／GitHub token 前綴／內網 IP／主機名稱／email／32+ 字元高熵字串
  3. 路徑清單檢查：`.env`、`*.pem`、`*.key`、`id_rsa`、`*.db`、`jobs/`、`credential*` 等
- 未對任何伺服器連線；未輪換任何憑證；未改寫歷史。

## 結論

**憑證類：零命中。** gitleaks 0 findings；手動掃描沒有任何硬編的 SECRET_KEY、token、API key、
密碼、webhook URL、私鑰，歷史中也從未出現 `.env`、`jobs/secret_key.txt`、`jobs/auth_tokens.db`、
`*.pem`、`*.key` 等檔案。所有憑證點都是從環境變數或 server 本機檔案讀取。

**非憑證類的資訊揭露：有，列於下表（低風險，但 repo 是 public）。**

## 命中清單

| # | 類型 | 內容（已遮蔽） | 首次出現 | 現行 HEAD 還在？ | 風險 |
|---|------|----------------|----------|------------------|------|
| 1 | 內網 IP | `172.16.0.1…`（HPC server，RFC1918 私有位址）| 946b587（2026-07-01，歷史最早 commit）起，多處 | 在（`CLAUDE.md`，含本次新增規則）| 低：私有位址無法從外網路由，但洩漏內網拓撲與 SSH 目標 |
| 2 | 內網 IP | `172.16.0.1…`（備用機）| 同上 | 在（`CLAUDE.md`）| 低 |
| 3 | server 使用者名稱與絕對路徑 | 使用者 `choukaihsuan`、`/home3/choukaihsuan/…`（105 行新增）、`/home/choukaihsuan/miniconda3/envs/…` | 946b587 | 在（`scripts/web_ui.py`、`benchmark/*`、`containers/build_and_deploy.sh`、`CLAUDE.md`）| 低～中：給出有效的 SSH 使用者名稱（搭配 #1 可做針對性爆破）|
| 4 | 本機路徑 | `/mnt/c/Users/User/…`（20 行）| 946b587 | 在（`CLAUDE.md`）| 低 |
| 5 | 個人 email | `chou.k…@gmail.com`（`NOTIFY_EMAIL_TO` 範例、git author）| 946b587 | 在 | 低：作者 email 本來就在 commit metadata（50 個 commit）|
| 6 | 佔位用密碼字串 | `NOTIFY_EMAIL_PASS="xxxx xxxx xxxx xxxx"`、`NOTIFY_SLACK_WEBHOOK="https://hooks.slack.com/services/..."` | 946b587 | 在（`CLAUDE.md`）| 無：確認為佔位文字，不是真值 |
| 7 | 預設值字串 | `RESEND_FROM` 預設 `onboarding@resend.dev`；`smtp.gmail.com:587` | 946b587 | 在 | 無：公開服務位址 |
| 8 | 服務位址 | `http://172.16.0.1…:5000` 寫在 `notify.py` 的 HTML 信件模板（歷史中，行 29586 的 patch 位置）| 946b587 | 不在 `scripts/` HEAD（`git grep` 於 scripts 無命中）| 低，僅歷史 |

「現行還在」＝存在於 HEAD；#8 是只存在於歷史的項目。

## 「必須輪換」清單

**憑證輪換：無。** 歷史中沒有任何曾經 commit 的憑證，因此不需要因為本稽核而輪換。

> 這個結論的前提是 Flask `SECRET_KEY`、Resend API key、Gmail App Password、SMTP 帳密只存在於
> server 的環境變數與 `jobs/secret_key.txt`，從未進入 git。以上掃描證實了 repo 這端沒有；
> 無法證實 server 端沒有其他外洩管道（例如 shell history、log），那不在本次範圍。

## 建議（不屬於輪換，需由你決定）

1. **`.gitignore` 沒有涵蓋 `jobs/`。** `web_ui.py` 會在 `jobs/secret_key.txt` 寫入 Flask `SECRET_KEY`、
   在 `jobs/auth_tokens.db` 存 magic-link token。目前歷史乾淨，但只要有人在 server 上 `git add -A`
   就會把它們推上 public repo。建議加 `jobs/`、`*.db`、`.env`、`secret_key*` 到 `.gitignore`
   （我會在任務 4 的 `.gitignore` 稽核中處理，除非你想先單獨處理）。
2. **內網資訊是否要從 CLAUDE.md 與 `scripts/` 拿掉？** #1–#4 都是低風險，且歷史中已無法收回
   （改寫歷史需 force push，這次不做）。若你在意，可在 HEAD 把 IP 與使用者名稱改成佔位符，
   但歷史版本仍可讀。是否處理由你決定。
3. 建議之後在 CI 加一個 gitleaks 步驟（`gitleaks detect --log-opts="--all"`），成本很低。
4. 若 server 的 SSH 對外網可達，建議確認僅允許金鑰登入；#3 讓使用者名稱公開。
