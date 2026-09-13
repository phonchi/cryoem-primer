"""Check the rendered original-notebook chapters on desktop and mobile."""
from pathlib import Path
import json
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'logs/image_revision_20260913/browser'
def main():
    errors=[]; rows=[]
    with sync_playwright() as p:
        browser=p.chromium.launch()
        page=browser.new_page()
        page.on('pageerror',lambda e:errors.append(str(e)))
        for width in (1440,390):
            page.set_viewport_size({'width':width,'height':1000})
            for chapter in ('01_image_basics','02_filter_segment','03_fourier','04_wavelet'):
                page.goto((ROOT/f'book/_build/html/{chapter}.html').as_uri(),wait_until='networkidle')
                if page.locator('.math').count():
                    page.wait_for_function('!window.MathJax || !MathJax.startup || !!document.querySelector("mjx-container")')
                broken=page.locator('img').evaluate_all('(xs)=>xs.filter(x=>!x.complete || x.naturalWidth===0).map(x=>x.src)')
                assert not broken,broken
                assert not page.evaluate('document.documentElement.scrollWidth>innerWidth+2'),(chapter,width)
                frames=[]
                for el in page.locator('iframe[srcdoc]').all():
                    el.scroll_into_view_if_needed();f=el.element_handle().content_frame();f.wait_for_load_state()
                    title=el.get_attribute('title');frames.append(title)
                    if title.startswith('影像處理互動'):
                        f.wait_for_function('window.ImageLesson && ImageLesson.last')
                        assert not f.evaluate('document.documentElement.scrollWidth>innerWidth+2'),title
                        if 'filter_sweep' in title:
                            for mode in ('low','high'):
                                f.locator('#'+mode).click()
                                for index in (0,11,23):
                                    f.locator('#cutoff').fill(str(index));f.locator('#cutoff').dispatch_event('input')
                                    state=f.evaluate('ImageLesson.last');assert state['cutoff']==index+1
                                    assert state['mode']==('lowpass' if mode=='low' else 'highpass')
                                    f.wait_for_function('document.querySelector("#filtered").complete')
                            el.screenshot(path=str(OUT/f'filter-sweep-{width}.png'))
                    if width==390: el.screenshot(path=str(OUT/f'{chapter}-{len(frames)}-mobile-site.png'))
                page.evaluate('scrollTo(0,0)');page.screenshot(path=str(OUT/f'{chapter}-{width}-site.png'))
                rows.append({'chapter':chapter,'width':width,'broken_images':broken,'frames':frames})
        browser.close()
    result={'pages':rows,'errors':errors};(OUT/'site_acceptance.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));assert not errors,errors
    print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
