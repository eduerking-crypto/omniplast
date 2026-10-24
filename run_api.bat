@echo off
title OMNIPLAST API REST (:8010)
cd /d "%~dp0"
python -m uvicorn api.main:app --host 127.0.0.1 --port 8010
pause
