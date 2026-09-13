"""Independent skimage references for the eight original geometry presets."""
from pathlib import Path
import numpy as np
from PIL import Image
from skimage import data
from skimage.transform import EuclideanTransform,SimilarityTransform,AffineTransform,ProjectiveTransform,warp,rotate
ROOT=Path(__file__).resolve().parents[1]
a=np.array(Image.open(ROOT/'book/images/original_notebooks/venusaur-f.png').convert('RGBA'))
a[40,75]=[0,0,0,255];a[40,76]=[0,255,0,255];rgb=a[:,:,:3]
transforms={'rigid':EuclideanTransform(rotation=np.pi/12,translation=(30,-20)),
'origin':EuclideanTransform(rotation=np.pi/12),
'similarity':SimilarityTransform(scale=.5,rotation=np.pi/12,translation=(30,-20)),
'shear':AffineTransform(shear=np.pi/6),
'projective':ProjectiveTransform(matrix=np.array([[1,-.5,40],[.1,.9,20],[.0015,.0015,1]]))}
out={}
for name,t in transforms.items():out[name]=warp(a,t.inverse,order=1);out['H_'+name]=t.params
out['center']=rotate(rgb,15,mode='constant',cval=255,preserve_range=True,order=1)/255
out['center50']=rotate(rgb,-50,mode='constant',cval=1,order=1)
for name,degree in [('center',15),('center50',-50)]:
 c=(np.array(a.shape[:2][::-1])-1)/2
 T=lambda x,y:np.array([[1,0,x],[0,1,y],[0,0,1.]])
 out['H_'+name]=T(*c)@EuclideanTransform(rotation=-np.deg2rad(degree)).params@T(*(-c))
src=np.array([[0,0],[0,50],[300,50],[300,0]])
dst=np.array([[155,15],[65,40],[260,130],[360,95]])
t=ProjectiveTransform.from_estimate(src,dst);out['text']=warp(data.text(),t,output_shape=(50,300),order=1);out['H_text']=t.params
folder=ROOT/'logs/image_revision_20260913/browser';folder.mkdir(parents=True,exist_ok=True)
np.savez_compressed(folder/'geometry_reference.npz',**out)
print('Saved eight original preset references')
