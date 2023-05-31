# Copyright 2007-2009, Sjoerd de Vries
# This file is part of the Spyder module: "models3d" 
# For licensing information, see LICENSE.txt 

def get_opposite(x):
  if x == 1: return 2
  if x == 2: return 1
  if x == 3: return 4
  if x == 4: return 3
  if x == 5: return 6
  if x == 6: return 5

def move(posx,posy,posz,cdir):
  if cdir == 1: return posx+1,posy,posz
  if cdir == 2: return posx-1,posy,posz
  if cdir == 3: return posx,posy+1,posz
  if cdir == 4: return posx,posy-1,posz
  if cdir == 5: return posx,posy,posz+1
  if cdir == 6: return posx,posy,posz-1

def turn(dir1,dir2,dir3,t):
  if t == 0: return(dir2,dir1,get_opposite(dir3))
  if t == 1: return(dir3,dir1,dir2)
  if t == 2: return(get_opposite(dir2),dir1,dir3)
  if t == 3: return(get_opposite(dir3),dir1,get_opposite(dir2))

