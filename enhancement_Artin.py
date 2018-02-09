import nibabel as nib
# import Image
import matplotlib.pylab as plt


from PIL import ImageEnhance , Image , ImageFilter
import numpy as np


def enhancing(im):
    im = Image.fromarray(im)

    im = im.convert('L')
    enh = ImageEnhance.Contrast(im)
    im = enh.enhance(1)
    return im


Directory = '/media/data1/artin/vimp2/'
name = 'wmnpre.nii.gz'

im = nib.load(Directory+name)
imD = im.get_data()
sz = imD.shape

imEnhanced = np.zeros(sz)
for i in range(sz[2]):
    imEnhanced[:,:,i] = enhancing(imD[:,:,i])


sliceNum = 40
print type(imEnhanced[:,:,sliceNum])
fig , axes = plt.subplots(1,2 , figsize=(10,5))
axes[0].imshow(imEnhanced[:,:,sliceNum],cmap='gray',aspect='auto')
axes[1].imshow(imD[:,:,sliceNum],cmap='gray',aspect='auto')
plt.show()
