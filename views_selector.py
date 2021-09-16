
import numpy as np
import matplotlib.pyplot as plt

array = list(range(60*8))


views = []

for i in array:
    views.append(1)

# turn of views
views_off = np.concatenate( [ np.arange(0,15) ,   np.arange(26,75) ,
                              np.arange(86,135), np.arange(146,195) ,
                              np.arange(206,255), np.arange(266,315),
                              np.arange(326,375), np.arange(386,435), np.arange(446,480) ]).tolist()

## turn of views NOT including 170 nm segment 
##views_off = np.concatenate( [ np.arange(0,35) ,   np.arange(46,95) ,
##                              np.arange(106,155), np.arange(166,215) ,
##                              np.arange(226,275), np.arange(286,335),
##                              np.arange(346,395), np.arange(406,455), np.arange(466,480) ]).tolist()
##

print(views_off)
for view in views_off:
    views[view] = 0




ar = np.reshape(views,(8,60))


plt.figure()
plt.imshow(ar)
plt.show()





# plot with BF
