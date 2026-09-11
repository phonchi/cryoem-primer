"""Inspect the built SPA pages; writes acceptance evidence, never source files."""
import json
from pathlib import Path
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
HTML = ROOT/'book/_build/html'
OUT = ROOT/'notes/spa_revision_20260911/browser'
PAGES = ['00_intro','05_background','05_image_formation','05_statistical_inference',
         '06_workflow','06_alignment_classification','06_dimension_reduction',
         '06_reconstruction_validation','06_heterogeneity','06_resolution_validation',
         '08_resources','appendix_conventions']


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    report=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True)
        try:
            for width,height in [(1440,1000),(390,844)]:
                page=browser.new_page(viewport={'width':width,'height':height})
                for stem in PAGES:
                    errors=[]
                    def on_error(e): errors.append(str(e))
                    page.on('pageerror',on_error)
                    page.goto((HTML/f'{stem}.html').as_uri(),wait_until='networkidle',timeout=60000)
                    if page.locator('article .math').count():
                        page.wait_for_function('Boolean(window.MathJax && MathJax.startup && MathJax.startup.promise)',timeout=30000)
                        page.evaluate('async () => { await MathJax.startup.promise; }')
                    page.locator('article').wait_for()
                    details=page.locator('article details')
                    if details.count():
                        item=details.first
                        item.locator('summary').click()
                        assert item.evaluate('(e)=>e.open'),stem
                        item.locator('summary').click()
                        assert not item.evaluate('(e)=>e.open'),stem
                    result=page.evaluate('''() => ({
                       overflow: document.documentElement.scrollWidth > innerWidth + 2,
                       mathErrors: [...document.querySelectorAll('mjx-merror')].map(e=>e.textContent),
                       mathCount: document.querySelectorAll('article mjx-container').length,
                       brokenImages: [...document.querySelectorAll('article img')].filter(e=>!e.complete || e.naturalWidth===0).map(e=>e.src),
                       next: document.querySelector('a.next-page')?.getAttribute('href') || null,
                       detailsCount: document.querySelectorAll('article details').length
                    })''')
                    result.update(page=stem,width=width,errors=errors)
                    report.append(result)
                    if stem in ['00_intro','05_statistical_inference','06_dimension_reduction','06_resolution_validation']:
                        page.screenshot(path=str(OUT/f'{stem}-{width}.png'),full_page=True)
                    page.remove_listener('pageerror',on_error)
                # Exercise one actual cross-chapter link in addition to checking all local URLs statically.
                page.goto((HTML/'05_statistical_inference.html').as_uri(),wait_until='networkidle')
                target=page.locator('article a[href="06_workflow.html"]').first
                target.click()
                assert page.url.endswith('/06_workflow.html')
                page.close()
        finally:
            browser.close()
    (OUT/'acceptance.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    failures=[r for r in report if r['overflow'] or r['mathErrors'] or r['brokenImages'] or r['errors']]
    print(json.dumps({'pages':len(PAGES),'viewports':[1440,390],'observations':len(report),'failures':failures},ensure_ascii=False,indent=2))
    return bool(failures)

if __name__=='__main__': raise SystemExit(main())
