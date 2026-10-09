"""gunicorn config for the CircDEX Web UI (scripts/web_ui.py).

Run (from the repo root, inside the conda env that has Flask):

    gunicorn -c deploy/gunicorn.conf.py --chdir scripts web_ui:app

IMPORTANT: keep ``workers = 1``. web_ui.py runs an in-process FIFO queue worker that launches
Snakemake; more than one worker process would start several queue workers and could launch the
same job twice. Concurrency comes from ``threads``. (Splitting the queue worker into its own
service is the clean long-term fix; see audit/security_hardening_20261008.md.)
"""
import os

bind = os.environ.get("CIRCDEX_BIND", "127.0.0.1:5000")   # nginx terminates TLS in front
workers = 1
threads = 4
worker_class = "gthread"
timeout = 120
graceful_timeout = 30
keepalive = 5
limit_request_line = 4094
limit_request_fields = 50
limit_request_field_size = 8190
forwarded_allow_ips = "127.0.0.1"          # only trust X-Forwarded-* from the local nginx
secure_scheme_headers = {"X-FORWARDED-PROTO": "https"}
accesslog = "-"
errorlog = "-"
# Never log query strings (they can carry job ids); keep method, path-only, status.
access_log_format = '%(h)s "%(m)s %(U)s" %(s)s %(b)s %(D)s'


def post_worker_init(worker):
    """Start the queue worker thread inside the (single) worker process."""
    import web_ui
    web_ui.start_background_workers()
