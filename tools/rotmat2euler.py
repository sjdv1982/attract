from math import *
def rotmat2euler(rotmat):
  phi = atan2(rotmat[1][2],rotmat[0][2])
  ssi = acos(rotmat[2][2])
  rot = atan2(-rotmat[2][1],-rotmat[2][0])       

  if fabs(rotmat[2][2]) >= 0.9999: #gimbal lock
    phi = 0
    if fabs(rotmat[0][0]) >= 0.9999:
      if rotmat[0][0] < 0: 
        rot = pi	
      else:
        rot = 0      
      if rotmat[2][2] < 0: 
        ssi = pi	
      else:
        ssi = 0
    else:
      if rotmat[2][2] < 0: 
        ssi = pi	
        rot = -acos(-rotmat[0][0])
      else:
        ssi = 0
        rot = acos(rotmat[0][0])
    if (rotmat[0][1] < 0): rot *= -1

  return phi,ssi,rot
