@echo off
title OMNIPLAST en la nube - Link publico via Cloudflare Tunnel
cd /d "%~dp0"
echo ============================================================
echo   OMNIPLAST: preparando tu link publico (Cloudflare)
echo   La app debe quedar corriendo en http://localhost:8787
echo ============================================================
powershell -NoProfile -Command "try{Invoke-WebRequest -Uri http://127.0.0.1:8787 -UseBasicParsing -TimeoutSec 5 | Out-Null; Write-Host 'App ya activa'}catch{Write-Host 'Iniciando app...'; Start-Process python -ArgumentList '-m','streamlit','run','app\Home.py','--server.port','8787','--server.headless','true' -WindowStyle Hidden; Start-Sleep 16}"
if not exist tools\cloudflared.exe (
    echo Descargando cloudflared...
    powershell -NoProfile -Command "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; New-Item -ItemType Directory -Path tools -Force | Out-Null; Invoke-WebRequest -Uri 'https://github.com/cloudflare/cloudflared/releases/latest/download/cloudflared-windows-amd64.exe' -OutFile 'tools\cloudflared.exe' -UseBasicParsing"
)
python -X utf8 scripts\abrir_tunel.py
echo Esperando el enlace publico (20 s)...
timeout /t 20 /nobreak >nul
powershell -NoProfile -Command "$m = (Get-Content data\tunnel_log.txt -ErrorAction SilentlyContinue | Out-String) | Select-String -Pattern 'https://[a-z0-9-]+\.trycloudflare\.com'; if($m){ $l = $m.Matches[0].Value; Set-Content -Path data\link_publico.txt -Value $l -Encoding UTF8; Write-Host (''); Write-Host ('   TU LINK PUBLICO: ' + $l) -ForegroundColor Green; Start-Process $l } else { Write-Host 'Aun no hay URL: revisa data\link_publico.txt en unos segundos' }"
echo.
echo NOTAS:
echo  - El link es aleatorio y cambia al reiniciar este script.
echo  - La app corre en TU PC: mantenlo encendido para que el link funcione.
echo  - Para un link fijo con tu dominio: tunel con nombre en dash.cloudflare.com
echo    (Zero Trust ^> Tunnels) usando el mismo tools\cloudflared.exe.
pause
