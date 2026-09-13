"""Execute only restored image chapters and save their notebook outputs."""
from pathlib import Path
import argparse,hashlib,json,time
import jupytext,nbformat
from nbclient import NotebookClient
ROOT=Path(__file__).resolve().parents[1]
NAMES=['01_image_basics','02_filter_segment','03_fourier','04_wavelet']
def main():
 p=argparse.ArgumentParser();p.add_argument('chapters',nargs='*');args=p.parse_args()
 report=[];logdir=ROOT/'logs/image_revision_20260913';logdir.mkdir(parents=True,exist_ok=True)
 for stem in args.chapters or NAMES:
  assert stem in NAMES
  source=ROOT/'book'/f'{stem}.py';nb=jupytext.read(source);start=time.time()
  print('START',stem,flush=True)
  try:
   NotebookClient(nb,timeout=300,kernel_name='cryoem-book',resources={'metadata':{'path':str(ROOT/'book')}},record_timing=False).execute()
  except Exception:
   nbformat.write(nb,logdir/f'failed_{stem}.ipynb');raise
  jupytext.write(nb,ROOT/'book'/f'{stem}.ipynb',fmt='ipynb')
  row={'chapter':stem,'seconds':round(time.time()-start,2),'cells':len(nb.cells),'code_cells':sum(c.cell_type=='code' for c in nb.cells),'errors':sum(o.output_type=='error' for c in nb.cells for o in c.get('outputs',[])),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
  report.append(row);print('DONE',json.dumps(row),flush=True)
  (logdir/'execution.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
