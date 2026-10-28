# gunicorn config for web_ui.py.
#
# Usage: gunicorn -c scripts/gunicorn_conf.py web_ui:app
#
# workers = 1 is a HARD REQUIREMENT, not a tuning knob. web_ui.py's Queue
# Worker is an in-process thread (see start_background_workers() in
# web_ui.py) that polls a SQLite job queue and runs Snakemake for whatever
# job it finds pending. Each gunicorn worker is a separate OS process; more
# than one worker means more than one Queue Worker thread polling the same
# queue independently, with no coordination between them -- the same job
# could be picked up and run twice, as two real Snakemake invocations
# against the same results_dir. If more concurrency is ever needed, use
# `threads` (below), never `workers`.
workers = 1
threads = 4
bind = "127.0.0.1:5000"
timeout = 120


def on_starting(server):
    """Runs once, in the master process, before any worker is forked."""
    import sys
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    from web_ui import run_preflight
    run_preflight()


def post_worker_init(worker):
    """Runs once per worker process, after it's forked.

    With workers = 1 this runs exactly once for the service's whole
    lifetime (per gunicorn start/reload) -- the single place
    start_background_workers() is invoked under gunicorn.
    """
    import sys
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    from web_ui import start_background_workers
    start_background_workers()
