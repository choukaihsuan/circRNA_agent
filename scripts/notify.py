#!/usr/bin/env python3
"""
circRNA Agent 通知模組
支援：SMTP Email + Slack Webhook

用法：
    python scripts/notify.py --event success --report results/report.html
    python scripts/notify.py --event failure --rule dcc --log logs/pipeline_run.log
    python scripts/notify.py --event start --project GSE113230

環境變數（加到 ~/.bashrc）：
    NOTIFY_EMAIL_FROM    Gmail 寄件地址
    NOTIFY_EMAIL_PASS    Gmail App Password（非登入密碼）
    NOTIFY_EMAIL_TO      收件地址
    NOTIFY_SMTP_HOST     (選填，預設 smtp.gmail.com)
    NOTIFY_SMTP_PORT     (選填，預設 587)
    NOTIFY_SLACK_WEBHOOK Slack Incoming Webhook URL（選填）

注意：LINE Notify 已於 2025-03-31 停止服務。
"""

from __future__ import annotations

import argparse
import json
import os
import smtplib
import urllib.request
from datetime import datetime
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path

# ── 設定（環境變數 > config.yaml fallback）────────────────────────────────────
SMTP_HOST     = os.getenv("NOTIFY_SMTP_HOST", "smtp.gmail.com")
SMTP_PORT     = int(os.getenv("NOTIFY_SMTP_PORT", "587"))
SMTP_USER     = os.getenv("NOTIFY_EMAIL_FROM", "")
SMTP_PASS     = os.getenv("NOTIFY_EMAIL_PASS", "")
SLACK_WEBHOOK = os.getenv("NOTIFY_SLACK_WEBHOOK", "")

def _email_to_from_config() -> str:
    """Read notify.email_to from config.yaml as fallback."""
    try:
        import yaml
        cfg_path = Path(__file__).parent.parent / "config.yaml"
        if cfg_path.exists():
            cfg = yaml.safe_load(cfg_path.read_text())
            return cfg.get("notify", {}).get("email_to", "")
    except Exception:
        pass
    return ""

EMAIL_TO = os.getenv("NOTIFY_EMAIL_TO", "") or _email_to_from_config()

_ATTACH_LIMIT_MB = 20  # Gmail attachment size limit


# ── Transport helpers ─────────────────────────────────────────────────────────

def send_email(subject: str, body: str, attachment_path: str = "") -> None:
    if not all([SMTP_USER, SMTP_PASS, EMAIL_TO]):
        print("[notify] Email 設定不完整（NOTIFY_EMAIL_FROM / NOTIFY_EMAIL_PASS / NOTIFY_EMAIL_TO），跳過")
        return

    msg = MIMEMultipart()
    msg["From"]    = SMTP_USER
    msg["To"]      = EMAIL_TO
    msg["Subject"] = subject
    msg.attach(MIMEText(body, "html", "utf-8"))

    if attachment_path and Path(attachment_path).exists():
        size_mb = Path(attachment_path).stat().st_size / 1e6
        if size_mb <= _ATTACH_LIMIT_MB:
            with open(attachment_path, "rb") as f:
                part = MIMEBase("application", "octet-stream")
                part.set_payload(f.read())
            encoders.encode_base64(part)
            part.add_header(
                "Content-Disposition",
                f"attachment; filename={Path(attachment_path).name}",
            )
            msg.attach(part)
        else:
            print(f"[notify] 報告 {size_mb:.1f} MB > {_ATTACH_LIMIT_MB} MB，略過附件")

    try:
        with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as s:
            s.ehlo()
            s.starttls()
            s.login(SMTP_USER, SMTP_PASS)
            s.send_message(msg)
        print(f"[notify] Email → {EMAIL_TO}")
    except Exception as exc:
        print(f"[notify] Email 失敗: {exc}")


def send_slack(text: str) -> None:
    if not SLACK_WEBHOOK:
        return
    payload = json.dumps({"text": text}).encode()
    req = urllib.request.Request(
        SLACK_WEBHOOK,
        data=payload,
        headers={"Content-Type": "application/json"},
    )
    try:
        urllib.request.urlopen(req, timeout=10)
        print("[notify] Slack 通知已發送")
    except Exception as exc:
        print(f"[notify] Slack 失敗: {exc}")


# ── Public API ────────────────────────────────────────────────────────────────

def notify_start(project_id: str) -> None:
    now  = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"🚀 [{project_id}] Pipeline 啟動（{now}）"
    send_slack(text)
    print(f"[notify] start: {text}")


def notify_success(project_id: str, report_path: str = "",
                   stats: dict | None = None) -> None:
    now   = datetime.now().strftime("%Y-%m-%d %H:%M")
    stats = stats or {}

    subject = f"[circRNA Agent] {project_id} 分析完成 ✅"
    body = f"""
<h2 style="color:#16a34a">✅ circRNA 分析完成</h2>
<table style="border-collapse:collapse;font-family:sans-serif">
  <tr><td style="padding:4px 12px 4px 0"><b>專案</b></td>
      <td>{project_id}</td></tr>
  <tr><td style="padding:4px 12px 4px 0"><b>完成時間</b></td>
      <td>{now}</td></tr>
  <tr><td style="padding:4px 12px 4px 0"><b>總 circRNA 數</b></td>
      <td>{stats.get("total_circ", "N/A")}</td></tr>
  <tr><td style="padding:4px 12px 4px 0"><b>顯著 DECs</b></td>
      <td>{stats.get("n_sig", "N/A")}</td></tr>
  <tr><td style="padding:4px 12px 4px 0"><b>上調</b></td>
      <td>{stats.get("n_up", "N/A")}</td></tr>
  <tr><td style="padding:4px 12px 4px 0"><b>下調</b></td>
      <td>{stats.get("n_down", "N/A")}</td></tr>
</table>
{f'<p>報告已附加，或至 server 查看：<br><code>{report_path}</code></p>' if report_path else ''}
"""
    slack_text = (
        f"✅ *{project_id}* 分析完成（{now}）\n"
        f"顯著 circRNA: {stats.get('n_sig','N/A')} 個 "
        f"（↑{stats.get('n_up','?')} ↓{stats.get('n_down','?')}）"
    )

    send_email(subject, body, report_path)
    send_slack(slack_text)


def notify_failure(project_id: str, failed_rule: str = "unknown",
                   log_path: str = "") -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    log_tail = ""
    if log_path and Path(log_path).exists():
        lines    = Path(log_path).read_text(errors="replace").splitlines()
        log_tail = "\n".join(lines[-50:])

    subject = f"[circRNA Agent] {project_id} 分析失敗 ❌ ({failed_rule})"
    body = f"""
<h2 style="color:#dc2626">❌ Pipeline 失敗</h2>
<p><b>專案：</b>{project_id}<br>
<b>失敗 rule：</b>{failed_rule}<br>
<b>時間：</b>{now}</p>
<h3>Log（最後 50 行）</h3>
<pre style="background:#f5f5f5;padding:10px;font-size:12px">{log_tail}</pre>
"""
    slack_text = (
        f"❌ *{project_id}* 失敗（{now}）\n"
        f"Rule: `{failed_rule}`"
    )

    send_email(subject, body)
    send_slack(slack_text)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="circRNA Agent 通知工具"
    )
    parser.add_argument("--event",   required=True,
                        choices=["success", "failure", "start"],
                        help="通知類型")
    parser.add_argument("--project", default="circRNA",
                        help="專案 ID（如 GSE113230）")
    parser.add_argument("--report",  default="",
                        help="report.html 路徑（success 用）")
    parser.add_argument("--rule",    default="unknown",
                        help="失敗的 rule 名稱（failure 用）")
    parser.add_argument("--log",     default="",
                        help="Pipeline log 檔路徑（failure 用）")
    parser.add_argument("--stats",   default="{}",
                        help='統計 JSON，如 \'{"n_sig":10,"n_up":6,"n_down":4}\'')
    args = parser.parse_args()

    try:
        stats = json.loads(args.stats)
    except json.JSONDecodeError:
        stats = {}

    if args.event == "start":
        notify_start(args.project)
    elif args.event == "success":
        notify_success(args.project, args.report, stats)
    elif args.event == "failure":
        notify_failure(args.project, args.rule, args.log)


if __name__ == "__main__":
    main()
