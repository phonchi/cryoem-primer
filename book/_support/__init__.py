"""Shared display helpers; computations remain in the teaching notebooks."""
from pathlib import Path
import base64,html,io,json,uuid
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from IPython.display import HTML,display
from IPython import get_ipython

_shell=get_ipython()
if _shell is not None:
    _shell.run_line_magic("matplotlib", "inline")

BOOK=Path(__file__).resolve().parents[1]

def image_path(name):
    path=BOOK/'images/original_notebooks'/name
    if not path.is_file():
        raise FileNotFoundError(f'Original notebook asset missing: {name}. Run scripts/fetch_original_assets.py first.')
    return path

def show_images(*images,titles=None,cmap=None,figsize=None,shared=False):
    titles=titles or ['']*len(images)
    fig,axes=plt.subplots(1,len(images),figsize=figsize or (4*len(images),4),squeeze=False)
    ranges={}
    if shared:
        ranges={'vmin':min(float(np.min(x)) for x in images),'vmax':max(float(np.max(x)) for x in images)}
    for ax,img,title in zip(axes.ravel(),images,titles):
        ax.imshow(img,cmap=cmap or ('gray' if np.ndim(img)==2 else None),**ranges)
        ax.set_title(title);ax.axis('off')
    fig.tight_layout();plt.show()
    return fig

def _image(array,signed=False,raw=False):
    a=np.asarray(array)
    low,high=float(np.nanmin(a)),float(np.nanmax(a))
    if raw:
        rendered=np.clip(a/(255 if a.dtype==np.uint8 else 1),0,1)
    elif a.ndim==2:
        if signed:
            high=max(abs(low),abs(high),1e-15);low=-high
        rendered=np.clip((a-low)/(high-low) if high>low else np.zeros_like(a),0,1)
    else:
        rendered=np.clip(a/(255 if a.dtype==np.uint8 else 1),0,1)
    picture=Image.fromarray(np.round(rendered*255).astype(np.uint8))
    stream=io.BytesIO();picture.save(stream,format='PNG')
    result={'src':'data:image/png;base64,'+base64.b64encode(stream.getvalue()).decode(),
            'width':a.shape[1],'height':a.shape[0],'min':low,'max':high}
    if raw:
        result.update(pixels=a.reshape(-1).tolist(),channels=1 if a.ndim==2 else a.shape[2],scale=255 if a.dtype==np.uint8 else 1)
    return result

def lab(kind,**data):
    """Render the same client-only interaction in notebooks and static HTML."""
    data=dict(data)
    if kind=='geometry':
        for key in ('image','rgb','text_image'):
            if key in data:data[key]=_image(data[key],raw=True)
    elif kind=='filter_sweep':
        for key in ('original','spectrum'):
            if key in data:data[key]=_image(data[key])
        for key in ('lowpass','highpass'):
            if key in data:data[key]=[_image(x,signed=key=='highpass') for x in data[key]]
    def encode(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    payload=json.dumps({'kind':kind,**data},default=encode,ensure_ascii=False).replace('<','\\u003c')
    folder=BOOK/'_support/image-lab'
    css=(folder/'lab.css').read_text();js=(folder/'lab.js').read_text()
    uid='image-lab-'+uuid.uuid4().hex[:12]
    page=f'''<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><style>{css}</style><body><main id="lab" aria-label="影像處理互動示範"></main><script id="lesson-data" type="application/json">{payload}</script><script>{js}</script><script>const frameId={json.dumps(uid)}; new ResizeObserver(()=>parent.postMessage({{imageLab:frameId,height:document.body.scrollHeight}},'*')).observe(document.body);</script></body></html>'''
    host="""<script>(function(){if(window.imageLabResizeReady)return;window.imageLabResizeReady=true;window.addEventListener('message',function(e){var d=e.data;if(!d||typeof d.imageLab!=='string')return;var f=document.getElementById(d.imageLab);if(f&&f.contentWindow===e.source&&Number.isFinite(d.height))f.style.height=Math.min(2400,Math.max(300,d.height+12))+'px';});})();</script>"""
    display(HTML(f'<iframe id="{uid}" title="影像處理互動：{html.escape(kind)}" class="image-lab-frame" srcdoc="{html.escape(page,quote=True)}" sandbox="allow-scripts" loading="lazy" style="width:100%;height:880px;border:0;border-radius:12px"></iframe>'+host))

def diagram(name):
    """Display an original, accessible SVG teaching schematic."""
    page=(BOOK/'_support/image-lab/diagrams'/f'{name}.html').read_text()
    uid='image-diagram-'+uuid.uuid4().hex[:12]
    observer=f"<script>new ResizeObserver(()=>parent.postMessage({{imageLab:{json.dumps(uid)},height:document.body.scrollHeight}},'*')).observe(document.body);</script>"
    page=page.replace('</body>',observer+'</body>')
    host="""<script>(function(){if(window.imageLabResizeReady)return;window.imageLabResizeReady=true;window.addEventListener('message',function(e){var d=e.data;if(!d||typeof d.imageLab!=='string')return;var f=document.getElementById(d.imageLab);if(f&&f.contentWindow===e.source&&Number.isFinite(d.height))f.style.height=Math.min(2400,Math.max(300,d.height+12))+'px';});})();</script>"""
    display(HTML(f'<iframe id="{uid}" title="教學關係圖：{html.escape(name)}" srcdoc="{html.escape(page,quote=True)}" sandbox="allow-scripts" loading="lazy" style="width:100%;height:760px;border:0"></iframe>'+host))
