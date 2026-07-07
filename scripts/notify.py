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
        print("[notify] Email not configured (NOTIFY_EMAIL_FROM / NOTIFY_EMAIL_PASS / NOTIFY_EMAIL_TO), skipping")
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
            print(f"[notify] Report {size_mb:.1f} MB > {_ATTACH_LIMIT_MB} MB, skipping attachment")

    try:
        with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as s:
            s.ehlo()
            s.starttls()
            s.login(SMTP_USER, SMTP_PASS)
            s.send_message(msg)
        print(f"[notify] Email → {EMAIL_TO}")
    except Exception as exc:
        print(f"[notify] Email failed: {exc}")


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
        print("[notify] Slack notification sent")
    except Exception as exc:
        print(f"[notify] Slack failed: {exc}")


# ── CircDEX Email Brand Header ────────────────────────────────────────────────

def _cd_email_header() -> str:
    return """
<div style="background:#0F2137;padding:20px 28px 16px;border-radius:8px 8px 0 0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif">
  <div style="display:inline-block">
    <span style="font-size:30px;font-weight:300;color:rgba(255,255,255,.85);letter-spacing:-.5px">Circ</span><span
          style="font-size:30px;font-weight:800;letter-spacing:.07em;color:#00B4C6">DEX</span>
  </div>
  <div style="font-size:14px;color:rgba(255,255,255,.55);margin-top:4px">
    From reads to circRNA biomarkers
  </div>
  <div style="margin-top:8px">
    <span style="display:inline-block;font-size:10px;font-weight:700;letter-spacing:.05em;text-transform:uppercase;color:#00B4C6;background:rgba(0,180,198,.15);border:1px solid rgba(0,180,198,.35);border-radius:4px;padding:3px 8px;margin-right:5px">Dual-tool consensus</span><span
          style="display:inline-block;font-size:10px;font-weight:700;letter-spacing:.05em;text-transform:uppercase;color:#00B4C6;background:rgba(0,180,198,.15);border:1px solid rgba(0,180,198,.35);border-radius:4px;padding:3px 8px;margin-right:5px">Differential expression</span><span
          style="display:inline-block;font-size:10px;font-weight:700;letter-spacing:.05em;text-transform:uppercase;color:#00B4C6;background:rgba(0,180,198,.15);border:1px solid rgba(0,180,198,.35);border-radius:4px;padding:3px 8px">6D biomarker ranking</span>
  </div>
</div>
<div style="height:1px;background:linear-gradient(90deg,transparent,#00B4C6,transparent);opacity:.5"></div>
"""


# ── Public API ────────────────────────────────────────────────────────────────

def notify_start(project_id: str) -> None:
    now  = datetime.now().strftime("%Y-%m-%d %H:%M")
    text = f"🚀 [{project_id}] Pipeline started ({now})"
    send_slack(text)
    print(f"[notify] start: {text}")


def notify_success(project_id: str, report_path: str = "",
                   stats: dict | None = None) -> None:
    now   = datetime.now().strftime("%Y-%m-%d %H:%M")
    stats = stats or {}

    subject = f"[CircDEX] {project_id} Analysis Complete ✅"
    body = f"""
{_cd_email_header()}
<div style="padding:24px 28px 28px;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif">
<h2 style="color:#16a34a;margin:0 0 16px">✅ Analysis Complete</h2>
<table style="border-collapse:collapse">
  <tr><td style="padding:4px 14px 4px 0;color:#64748b;font-size:13px"><b>Project</b></td>
      <td style="font-size:13px">{project_id}</td></tr>
  <tr><td style="padding:4px 14px 4px 0;color:#64748b;font-size:13px"><b>Completed at</b></td>
      <td style="font-size:13px">{now}</td></tr>
  <tr><td style="padding:4px 14px 4px 0;color:#64748b;font-size:13px"><b>Total circRNAs</b></td>
      <td style="font-size:13px">{stats.get("total_circ", "N/A")}</td></tr>
  <tr><td style="padding:4px 14px 4px 0;color:#64748b;font-size:13px"><b>Significant DECs</b></td>
      <td style="font-size:13px">{stats.get("n_sig", "N/A")}</td></tr>
  <tr><td style="padding:4px 14px 4px 0;color:#64748b;font-size:13px"><b>Up-regulated</b></td>
      <td style="font-size:13px">{stats.get("n_up", "N/A")}</td></tr>
  <tr><td style="padding:4px 14px 4px 0;color:#64748b;font-size:13px"><b>Down-regulated</b></td>
      <td style="font-size:13px">{stats.get("n_down", "N/A")}</td></tr>
</table>
{f'<p style="margin-top:16px;font-size:13px">Report attached, or view on server:<br><code style="background:#f1f5f9;padding:2px 6px;border-radius:4px">{report_path}</code></p>' if report_path else ''}
</div>
"""
    slack_text = (
        f"✅ *{project_id}* analysis complete ({now})\n"
        f"Significant circRNAs: {stats.get('n_sig','N/A')} "
        f"(↑{stats.get('n_up','?')} ↓{stats.get('n_down','?')})"
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

    subject = f"[CircDEX] {project_id} Analysis Failed ❌ ({failed_rule})"
    body = f"""
{_cd_email_header()}
<div style="padding:24px 28px 28px;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif">
<h2 style="color:#dc2626;margin:0 0 12px">❌ Pipeline Failed</h2>
<p style="font-size:13px"><b>Project:</b> {project_id}<br>
<b>Failed rule:</b> {failed_rule}<br>
<b>Time:</b> {now}</p>
<h3 style="margin:16px 0 8px;font-size:13px">Log (last 50 lines)</h3>
<pre style="background:#f5f5f5;padding:12px;font-size:11px;border-radius:6px;overflow-x:auto">{log_tail}</pre>
</div>
"""
    slack_text = (
        f"❌ *{project_id}* failed ({now})\n"
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
