import asyncio
import json
from pathlib import Path

from playwright.async_api import async_playwright

CHROME = Path(r"C:\Users\reyes\AppData\Local\ms-playwright\chromium-1243\chrome-win64\chrome.exe")


async def main():
    async with async_playwright() as pw:
        browser = await pw.chromium.launch(
            executable_path=str(CHROME) if CHROME.exists() else None,
            headless=True, args=["--disable-gpu"])
        page = await browser.new_page(viewport={"width": 1440, "height": 900})
        await page.goto("http://127.0.0.1:8787", wait_until="networkidle", timeout=45000)
        await page.wait_for_timeout(4000)
        body = await page.inner_text("body")
        es_omniplast = "OMNIPLAST" in body and "microplásticos" in body.lower()
        enlaces = await page.query_selector_all('section[data-testid="stSidebar"] a')
        nombres = [((await e.inner_text()).strip()) for e in enlaces]
        print(json.dumps({"es_omniplast": es_omniplast,
                          "paginas_menu": len(nombres),
                          "menu": nombres[:6],
                          "hero": body[:90].replace("\n", " ")},
                         ensure_ascii=False, indent=2))
        await browser.close()

asyncio.run(main())
