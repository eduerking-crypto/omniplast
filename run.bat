@echo off
title OMNIPLAST - Plataforma de Analisis Integral de Microplasticos
cd /d "%~dp0"
echo ============================================================
echo   OMNIPLAST v1.0.0 - Analisis doctoral de microplasticos
echo   Abriendo la plataforma en http://localhost:8787
echo   (puerto propio: 8501 suele estar ocupado por otras apps)
echo   Para detener: cierre esta ventana o presione Ctrl+C
echo ============================================================
python -m streamlit run app\Home.py --server.port 8787 --browser.gatherUsageStats false
pause
