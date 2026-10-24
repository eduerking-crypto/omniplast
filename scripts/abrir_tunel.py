"""Lanza el tunel publico de OMNIPLAST (quick tunnel, sin cuenta) de forma detached."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
log_path = ROOT / "data" / "tunnel_log.txt"

for pid in subprocess.run(
        ["powershell", "-NoProfile", "-Command",
         "Get-Process cloudflared -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Id"],
        capture_output=True, text=True).stdout.split():
    subprocess.run(["powershell", "-NoProfile", "-Command",
                    f"Stop-Process -Id {pid} -Force -ErrorAction SilentlyContinue"])

log = open(log_path, "w", encoding="utf-8", buffering=1)
DETACHED = 0x00000008
NEW_GROUP = 0x00000200
proc = subprocess.Popen(
    [str(ROOT / "tools" / "cloudflared.exe"), "tunnel",
     "--url", "http://localhost:8787", "--no-autoupdate"],
    stdout=log, stderr=subprocess.STDOUT, cwd=str(ROOT),
    creationflags=DETACHED | NEW_GROUP, close_fds=True)
print("tunel lanzado, PID:", proc.pid)
sys.exit(0)
