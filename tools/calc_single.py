# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 12:02:37 2021

@author: Zsuzsanna Koczor-Benda, UCL
"""

import numpy as np
import math

def single_rot_IR(d,e,a=0,b=0,c=0):
    Rz=np.array([[math.cos(a), -math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]])
    Rx=np.array([[1,0,0],[0,math.cos(b), -math.sin(b)],[0,math.sin(b),math.cos(b)]])
    Rz2=np.array([[math.cos(c), -math.sin(c),0],[math.sin(c),math.cos(c),0],[0,0,1]])
    R=np.matmul(Rz2,np.matmul(Rx,Rz))
    q=np.matmul(np.transpose(R),e)
    ir=math.pow(np.matmul(q,d),2)
    return ir

def single_rot_R(p,e_in,e_out,a=0,b=0,c=0):
    Rz=np.array([[math.cos(a), -math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]])
    Rx=np.array([[1,0,0],[0,math.cos(b), -math.sin(b)],[0,math.sin(b),math.cos(b)]])
    Rz2=np.array([[math.cos(c), -math.sin(c),0],[math.sin(c),math.cos(c),0],[0,0,1]])
    R=np.matmul(Rz2,np.matmul(Rx,Rz))
    q1=np.matmul(np.transpose(R),e_in)
    q2=np.matmul(np.transpose(R),e_out)
    raman=math.pow(np.matmul(q1,np.matmul(p,np.transpose(q2))),2)
    return raman

def single_polar_IR(d,r,e):
    q=np.matmul(np.transpose(r),e)
    ir=math.pow(np.matmul(q,d),2)
    return ir

def single_polar_R(p,r,e_in,e_out):
    q1=np.matmul(np.transpose(r),e_in)
    q2=np.matmul(np.transpose(r),e_out)
    raman=math.pow(np.matmul(q1,np.matmul(p,np.transpose(q2))),2)
    return raman