#!/usr/bin/env python

class atom_data:
  def __init__(self):
    ntype=0 #types of atoms
    natoms=0 #numbers of atoms
    symbol=None
    atomSpe=None
    position=None
    constraint=None

class dos_data:
  def __init__(self):
    efermi=0
    ef=None
    dos0=None
    perelement=None
    taskname=""
    ncol=0
    par_element=[-1,-1,-1,-1] #partial orbital
    par_orbital=None

