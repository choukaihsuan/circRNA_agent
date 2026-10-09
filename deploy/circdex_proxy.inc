# /etc/nginx/circdex_proxy.inc  — shared proxy settings
proxy_pass         http://127.0.0.1:5000;
proxy_http_version 1.1;
proxy_set_header   Host              $host;            # the REAL host; web_ui ignores X-Forwarded-Host
proxy_set_header   X-Real-IP         $remote_addr;
proxy_set_header   X-Forwarded-For   $remote_addr;     # overwrite, never append client-supplied values
proxy_set_header   X-Forwarded-Proto $scheme;
proxy_set_header   X-Forwarded-Host  "";
proxy_read_timeout 120s;
