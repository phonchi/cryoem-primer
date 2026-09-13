"""Browser + numerical acceptance of the original-example interactions."""
from pathlib import Path
import json
import numpy as np
from scipy.signal import convolve2d
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'logs/image_revision_20260913/browser';OUT.mkdir(parents=True,exist_ok=True)

def main():
    report=[];errors=[]
    references=np.load(OUT/"geometry_reference.npz")
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True)
        page=browser.new_page(viewport={'width':1280,'height':1000})
        page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto((ROOT/'notes/image_revision_20260913/interactive_preview.html').as_uri(),wait_until='networkidle')
        def get_frame(kind):
            element=page.locator(f'iframe[title="影像處理互動：{kind}"]')
            element.scroll_into_view_if_needed()
            frame=element.element_handle().content_frame()
            frame.wait_for_function('window.ImageLesson && ImageLesson.last')
            return frame,element
        g,ge=get_frame('geometry')
        for mode in ['rigid','origin','center','center50','similarity','shear','projective','text']:
            g.locator(f'[data-mode="{mode}"]').click()
            info=g.evaluate('ImageLesson.last');data=g.evaluate('ImageLesson.data')
            image=data['text_image' if mode=='text' else 'rgb' if mode.startswith('center') else 'image']
            a=np.asarray(image['pixels'],dtype=np.uint8).reshape((image['height'],image['width'])+(() if image['channels']==1 else (image['channels'],)))
            actual=np.asarray(g.evaluate('Array.from(document.querySelector("#geo-output").getContext("2d").getImageData(0,0,ImageLesson.last.width,ImageLesson.last.height).data)'),dtype=float).reshape(info['height'],info['width'],4)
            expected=references[mode]
            np.testing.assert_allclose(info['H'],references['H_'+mode],atol=1e-9)
            if expected.ndim==2:expected=np.repeat(expected[:,:,None],3,axis=2)
            if expected.shape[2]==3:expected=np.concatenate([expected,np.ones((*expected.shape[:2],1))],axis=2)
            # Compare displayed RGBA after canvas alpha quantization. Fully transparent RGB is immaterial.
            exp=np.round(expected*255)
            alpha_error=float(np.max(np.abs(exp[:,:,3]-actual[:,:,3])))
            displayed_expected=expected[:,:,:3]*expected[:,:,3:4]*255+255*(1-expected[:,:,3:4])
            displayed_actual=actual[:,:,:3]*(actual[:,:,3:4]/255)+255-actual[:,:,3:4]
            display_error=float(np.max(np.abs(displayed_expected-displayed_actual)))
            row={'kind':'geometry','preset':mode,'alpha_error':alpha_error,'display_error':display_error,'shape':[info['height'],info['width']]}
            report.append(row)
            assert alpha_error<=1 and display_error<=2,row
        g.locator('[data-mode="projective"]').click();g.locator('#matrix-edit').fill('1 0 0\n0 1 0\n0 0 0');g.locator('#apply-matrix').click();assert g.locator('#geometry-error').inner_text()
        g.locator('#reset').click();assert g.locator('#geometry-error').inner_text()==''
        g.locator('[data-mode="rigid"]').click();ge.screenshot(path=str(OUT/'geometry-desktop.png'))
        b,be=get_frame('fourier_basis')
        for kx,ky,a,kind in [(5,2,1,'cosine'),(-32,17,.4,'sine'),(0,0,.5,'cosine'),(0,0,.8,'sine')]:
            b.evaluate('([x,y,a,w])=>{for(const [id,v] of [["kx",x],["ky",y],["amplitude",a]])document.getElementById(id).value=v;document.getElementById("wave").value=w;document.getElementById("wave").dispatchEvent(new Event("change"));}',[kx,ky,a,kind])
            actual=np.array(b.evaluate('ImageLesson.last.values')).reshape(128,128)
            H=np.zeros((128,128),complex)
            if kx==ky==0:H[0,0]=a if kind=='cosine' else 0
            else:
                weight=a if kind=='cosine' else -1j*a;H[ky,kx]=weight;H[-ky,-kx]=np.conj(weight)
            err=float(np.max(np.abs(actual-np.fft.ifft2(H).real)));report.append({'kind':'basis','parameters':[kx,ky,a,kind],'error':err});assert err<1e-12
        b.locator('#reset').click();be.screenshot(path=str(OUT/'fourier-desktop.png'))
        c,ce=get_frame('convolution1d');d=c.evaluate('ImageLesson.data')
        for key,k in d['kernels'].items():
            c.locator('#kernel').select_option(key);actual=c.evaluate('ImageLesson.last.values');expected=np.convolve(d['signal'],k,'valid');err=float(np.max(np.abs(actual-expected)));report.append({'kind':'convolution1d','kernel':key,'error':err});assert err<1e-12
        c.locator('#play').click();c.wait_for_timeout(650);assert c.evaluate('ImageLesson.last.index')>0;c.locator('#play').click();c.locator('#reset').click();assert c.evaluate('ImageLesson.last.index')==0
        c2,c2e=get_frame('convolution2d');d=c2.evaluate('ImageLesson.data');actual=c2.evaluate('ImageLesson.last.values');err=float(np.max(np.abs(actual-convolve2d(d['image'],d['kernel'],mode='same'))));report.append({'kind':'convolution2d','error':err});assert err<1e-12
        c2.locator('#position').fill('16');c2.locator('#position').dispatch_event('input');assert np.isclose(c2.evaluate('ImageLesson.last.value'),4/9)
        c2e.screenshot(path=str(OUT/'convolution-desktop.png'))
        page.set_viewport_size({'width':390,'height':844})
        for kind in ['geometry','fourier_basis','convolution1d','convolution2d']:
            f,e=get_frame(kind)
            assert not f.evaluate('document.documentElement.scrollWidth>innerWidth+2'),kind
            e.screenshot(path=str(OUT/f'{kind}-mobile.png'))
        browser.close()
    result={'checks':report,'page_errors':errors};(OUT/'labs_acceptance.json').write_text(json.dumps(result,indent=2)+'\n');assert not errors,errors
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
