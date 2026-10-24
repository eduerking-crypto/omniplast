import asyncio
import json
from pathlib import Path

from playwright.async_api import async_playwright

BASE = "http://127.0.0.1:7965"
CHROME = Path(r"C:\Users\reyes\AppData\Local\ms-playwright\chromium-1243\chrome-win64\chrome.exe")

PROBE = """
() => {
  const sb = document.querySelector('section[data-testid="stSidebar"]');
  const hero = document.querySelector('.omp-hero');
  const h1 = document.querySelector('h1');
  const sheet = document.getElementById('omniplast-theme');
  return {
    sidebar_bg: sb ? getComputedStyle(sb).backgroundColor : 'x',
    hero_bg: hero ? getComputedStyle(hero).backgroundImage.slice(0, 60) : 'x',
    h1_font: h1 ? getComputedStyle(h1).fontFamily.slice(0, 40) : 'x',
    style_omp_presente: !!sheet,
    n_reglas: sheet && sheet.sheet ? sheet.sheet.cssRules.length : 0
  };
}
"""


async def main():
    async with async_playwright() as pw:
        browser = await pw.chromium.launch(
            executable_path=str(CHROME) if CHROME.exists() else None,
            headless=True, args=["--disable-gpu"])
        page = await browser.new_page(viewport={"width": 1440, "height": 900})
        await page.goto(BASE, wait_until="networkidle", timeout=45000)
        await page.wait_for_timeout(4500)
        r = await page.evaluate(PROBE)
        print(json.dumps(r, indent=2))
        await page.screenshot(path="data/qa/tema_home.png")
        await browser.close()

asyncio.run(main())
