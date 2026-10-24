import asyncio
import json
from pathlib import Path

from playwright.async_api import async_playwright

BASE = "http://127.0.0.1:7965"
OUT = Path(__file__).parents[1] / "data" / "qa"
OUT.mkdir(parents=True, exist_ok=True)
CHROME = Path(r"C:\Users\reyes\AppData\Local\ms-playwright\chromium-1243\chrome-win64\chrome.exe")


async def main():
    reporte = {}
    js_errors = []
    async with async_playwright() as pw:
        browser = await pw.chromium.launch(
            executable_path=str(CHROME) if CHROME.exists() else None,
            headless=True, args=["--disable-gpu"])
        page = await browser.new_page(viewport={"width": 1440, "height": 900})
        page.on("pageerror", lambda e: js_errors.append(str(e)[:300]))
        await page.goto(BASE, wait_until="networkidle", timeout=45000)
        await page.wait_for_timeout(4000)
        enlaces = await page.query_selector_all('section[data-testid="stSidebar"] a')
        nombres = []
        for e in enlaces:
            nombres.append((await e.inner_text()).strip().replace("\n", " "))
        print("Links encontrados:", len(nombres))
        print(nombres)
        for i, e in enumerate(enlaces):
            nombre = nombres[i] if i < len(nombres) else f"pagina_{i}"
            try:
                await e.click(timeout=15000)
            except Exception as err:
                reporte[nombre] = {"error_click": str(err)[:200]}
                continue
            await page.wait_for_timeout(6500)
            exc = await page.query_selector('[data-testid="stException"]')
            texto_exc = ""
            if exc:
                texto_exc = (await exc.inner_text())[:400]
            body = await page.inner_text("body")
            reporte[nombre] = {
                "excepcion": texto_exc,
                "url": page.url.replace(BASE, ""),
                "tiene_hero_o_header": ("OMNIPLAST" in body) or ("Panel" in body) or len(body) > 800,
                "n_caracteres": len(body),
            }
            safe = "".join(c if c.isalnum() else "_" for c in nombre)[:40]
            await page.screenshot(path=str(OUT / f"v2_{i:02d}_{safe}.png"))
        await browser.close()
    (OUT / "reporte_qa_v2.json").write_text(json.dumps(reporte, indent=2, ensure_ascii=False),
                                            encoding="utf-8")
    fallas = {k: v for k, v in reporte.items() if v.get("excepcion")}
    cortas = {k: v for k, v in reporte.items() if v.get("n_caracteres", 0) < 300}
    print("PAGINAS:", len(reporte), "| con excepcion:", len(fallas), "| casi vacias:", len(cortas))
    for k, v in fallas.items():
        print(f"\n[{k}] {v['excepcion'][:280]}")
    for k, v in cortas.items():
        print(f"[VACIA?] {k}: {v}")
    print("JS errors:", len(js_errors))
    for e in js_errors[:6]:
        print("  -", e[:180])


asyncio.run(main())
