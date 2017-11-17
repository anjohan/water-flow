a=7.16
antx=30
anty=30
x=$(py "$antx*$a")
y=$(py "$anty*$a")
nederst=0
hardbunntykkelse=$(py "2*$a")
mykbunntykkelse=$(py "6*$a")
xmid=$(py "$x/2")
ymid=$(py "$y/2")
midtenavbunn=$(py "${nederst} + ${hardbunntykkelse}")
toppavbunn=$(py "${midtenavbunn} + ${mykbunntykkelse}")
sylinderhoyde=$(py "4*$a")
toppavsylinder=$(py "${toppavbunn} + ${sylinderhoyde} + $R")
sylinderradius=$(py "7*$a")

T=300
