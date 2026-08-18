#!/usr/bin/env python3
# Cell-center (x,y,z) for an AMR cell (i,j,k) at a given level, per the amrcomp_impact
# geometry convention in simulation.f90 (xlo/xhi/ylo/yhi fixed, zlo/zhi derived for nz=1,
# ref ratio 2 per level in x/y, ref ratio 1 in z).

import argparse

def cell_center(i,j,k,level,xlo,xhi,ylo,yhi,zlo,zhi,nx,ny,nz,ref_xy=2,ref_z=1):
   dx=(xhi-xlo)/(nx*ref_xy**level)
   dy=(yhi-ylo)/(ny*ref_xy**level)
   dz=(zhi-zlo)/(nz*ref_z**level)
   x=xlo+(i+0.5)*dx
   y=ylo+(j+0.5)*dy
   z=zlo+(k+0.5)*dz
   return x,y,z

def main():
   p=argparse.ArgumentParser(description='AMR cell indices (i,j,k) -> physical cell-center coordinates')
   p.add_argument('i',type=int); p.add_argument('j',type=int); p.add_argument('k',type=int)
   p.add_argument('--level',type=int,default=6,help='AMR level (default: finest, 6)')
   p.add_argument('--xlo',type=float,default=0.0)
   p.add_argument('--xhi',type=float,default=20.0,help='Domain length')
   p.add_argument('--ylo',type=float,default=-10.0)
   p.add_argument('--yhi',type=float,default=10.0)
   p.add_argument('--nx',type=int,default=32)
   p.add_argument('--ny',type=int,default=32)
   p.add_argument('--nz',type=int,default=1)
   p.add_argument('--maxlvl',type=int,default=6,help='Max AMR level, used to auto-derive zlo/zhi when nz=1')
   p.add_argument('--zlo',type=float,default=None)
   p.add_argument('--zhi',type=float,default=None)
   p.add_argument('--ref-xy',type=int,default=2,help='Refinement ratio per level in x,y')
   p.add_argument('--ref-z',type=int,default=1,help='Refinement ratio per level in z')
   args=p.parse_args()

   zlo,zhi=args.zlo,args.zhi
   if zlo is None or zhi is None:
      if args.nz!=1:
         raise SystemExit('--zlo/--zhi required when nz != 1')
      half=0.5*(args.yhi-args.ylo)/(args.ny*2**args.maxlvl)
      zlo,zhi=-half,half

   x,y,z=cell_center(args.i,args.j,args.k,args.level,args.xlo,args.xhi,args.ylo,args.yhi,zlo,zhi,
                      args.nx,args.ny,args.nz,args.ref_xy,args.ref_z)
   print(f'x={x!r}')
   print(f'y={y!r}')
   print(f'z={z!r}')

if __name__=='__main__':
   main()
