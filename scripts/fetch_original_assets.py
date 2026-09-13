"""Fetch exactly the input/illustration assets referenced by original notebooks."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import base64,hashlib,io,json,re
from urllib.parse import urlsplit
import requests
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT.parent
DEST=ROOT/'book/images/original_notebooks'
REPORT=ROOT/'notes/image_revision_20260913/assets.json'

def inventory():
    items={}
    for nbpath in sorted(SOURCE.glob('0*_zh.ipynb')):
        chapter=nbpath.name[:2]
        nb=json.loads(nbpath.read_text())
        for index,cell in enumerate(nb['cells']):
            text=''.join(cell['source'])
            pairs=[]
            if cell['cell_type']=='code':
                pairs += [(name,url) for url,name in re.findall(r"download_from_pokemondb\(\s*['\"]([^'\"]+)['\"]\s*,\s*['\"]([^'\"]+)['\"]",text)]
                pairs += [(urlsplit(url).path.rsplit('/',1)[-1],url) for url in re.findall(r'!wget\s+(https?://\S+)',text)]
            else:
                matches=[(m.start(),m.group(1)) for m in re.finditer(r'<img\b[^>]*?src=[\"\']([^\"\']+)',text)]
                matches += [(m.start(),m.group(1)) for m in re.finditer(r'!\[[^\]]*\]\(([^\s)]+)\)',text)]
                pairs += [(f'ch{chapter}-cell{index}-{j}.png',url) for j,(_,url) in enumerate(sorted(matches),1)]
            for name,url in pairs:
                if name in items and items[name]['url']!=url: raise ValueError(f'conflicting original asset {name}')
                row=items.setdefault(name,{'name':name,'url':url,'sources':[]})
                row['sources'].append({'notebook':nbpath.name,'cell':index})
    return list(items.values())

def fetch(row):
    url=row['url'];path=DEST/row['name']
    if url.startswith('data:'):
        blob=base64.b64decode(url.split(',',1)[1]); public_url='embedded original notebook image'
    else:
        headers={"User-Agent":"cryoem-primer/1.0 (educational notebook restoration; https://phonchi.github.io/cryoem-primer/)"} if "wikimedia.org" in urlsplit(url).netloc else {}
        with requests.get(url,stream=True,timeout=40,headers=headers) as response:
            response.raise_for_status();response.raw.decode_content=True;blob=response.raw.read()
        public_url=url
    with Image.open(io.BytesIO(blob)) as im:
        fmt=im.format;size=im.size;mode=im.mode
        im.verify()
    path.write_bytes(blob)
    return {**row,'url':public_url,'path':str(path.relative_to(ROOT)),
            'sha256':hashlib.sha256(blob).hexdigest(),'format':fmt,'size':size,'mode':mode,'bytes':len(blob)}

def main():
    DEST.mkdir(parents=True,exist_ok=True)
    rows=inventory();results=[];failures=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        jobs=[(row,pool.submit(fetch,row)) for row in rows]
        for row,job in jobs:
            try:
                result=job.result();results.append(result);print('OK',result['name'],result['size'],result['format'],flush=True)
            except Exception as e:
                failures.append({'name':row['name'],'url':row['url'][:200],'error':str(e)});print('FAILED',row['name'],str(e),flush=True)
    REPORT.write_text(json.dumps({'assets':results,'failures':failures},ensure_ascii=False,indent=2)+'\n')
    if failures:raise SystemExit(1)
    print('all original assets:',len(results))
if __name__=='__main__':main()
