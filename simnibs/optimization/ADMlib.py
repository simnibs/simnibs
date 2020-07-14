#!/usr/local/bin/python
import os
import numpy as np
import fmm3dpy as fmm
import time
import numpy.linalg as la
#*********************************************************************
#  This library compute E-fields induced during TMS via reciprocity and auxiliary
#  dipoles. The method is implemented for magnetic dipole sources.
#  companion pre-print: https://www.biorxiv.org/content/10.1101/2020.05.27.120022v2
#  created by Luis Gomez (July 2020)!
#  Last modified by Luis Gomez (July 2020)
#*********************************************************************
def recipcode(rv,jv,rs,ks,A):
	#this function computes E-fields via reciprocity no auxiliary dipoles
	#rv is 3 by ntetra and has mesh tetrahedron centroid positions
	#jv is 3 by ntetra and has total conduction current at each tetrahedron
	#rs is 3 by ncoil and has the coil dipole positions centered about the origin
	#ks is 3 by ncoil and has the coil dipole weights
	#A is 4 by 4 by number of coilpositions and each 4 by 4 matrix is a translation of the coil to a point above the scalp

	#paramameter:
	prec=10**(-3); #this number determines the accuracy of H-primary evaluation higher accuracy=slower time increases ~log(accuracy)

	#reads in sizes of each array
	npos=A.shape[2]; #number of scalp positions
	ncoil=rs.shape[1]; 
	nsource=rv.shape[1];
	ntarget=ncoil*npos;
	
	#generate copies of coil
	rp=np.transpose(np.c_[np.transpose(rs) ,np.ones([ncoil])]); #pads coil positions to 4 by ncoil
	robs=np.array(np.zeros([3,ntarget]));
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		robs[:,st:en]=np.matmul(A[0:3,:,i],rp);

	start = time.time()
	print("Computing H-primary");
	Hprimary=computeHprimary(rv,jv,robs,prec);
	end = time.time()
	print("H-primary time: ",end-start);	
	kp=ks;
	Etotal=np.zeros(npos);
	print(Hprimary)
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		kp=np.matmul(A[0:3,0:3,i],ks);
		Etotal[i]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	return Etotal;

def ADM(rv,jv,rs,ks,A,coildir):
	#this function computes E-fields via reciprocity with auxiliary dipoles
	#rv is 3 by ntetra and has mesh tetrahedron centroid positions
	#jv is 3 by ntetra and has total conduction current at each tetrahedron
	#rs is 3 by ncoil and has the coil dipole positions centered about the origin
	#ks is 3 by ncoil and has the coil dipole weights
	#A is 4 by 4 by number of coilpositions and each 4 by 4 matrix is a translation of the coil to a point above the scalp
	#coildir is a 3 by number of coil orientations that gives the y-direction orientation


	#paramameter:
	prec=10**(-3); #this number determines the accuracy of H-primary evaluation higher accuracy=slower time increases ~log(accuracy)
	N=[17,17,2];
	#generate auxiliary dipoles
	Nj=coildir.shape[1];
	raux,kaux=resamplecoil(rs,ks,N,Nj,coildir);
	#reads in sizes of each array
	npos=A.shape[2]; #number of scalp positions
	ncoil=raux.shape[1]; 
	nsource=rv.shape[1];
	ntarget=ncoil*npos;
	
	#generate copies of coil
	rp=np.transpose(np.c_[np.transpose(raux) ,np.ones([ncoil])]); #pads coil positions to 4 by ncoil
	robs=np.array(np.zeros([3,ntarget]));
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		robs[:,st:en]=np.matmul(A[0:3,:,i],rp);

	start = time.time()
	print("Computing H-primary");
	Hprimary=computeHprimary(rv,jv,robs,prec);
	end = time.time()
	print("H-primary time: ",end-start);	
	kp=kaux[:,:,0];
	
	Etotal=np.zeros([Nj,npos]);
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		for j in range(0,Nj):
			kp=np.matmul(A[0:3,0:3,i],kaux[:,:,j]);#for flat coils you only need the 3rd element
			Etotal[j,i]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	return Etotal;

def recipcodemag(rv,jvx,jvy,jvz,rs,ks,A):
	#this function computes E-field unidirectional approximation of the magnitude via reciprocity no auxiliary dipoles
	#rv is 3 by ntetra and has mesh tetrahedron centroid positions
	#jv is 3 by ntetra and has total conduction current at each tetrahedron
	#rs is 3 by ncoil and has the coil dipole positions centered about the origin
	#ks is 3 by ncoil and has the coil dipole weights
	#A is 4 by 4 by number of coilpositions and each 4 by 4 matrix is a translation of the coil to a point above the scalp

	#paramameter:
	prec=10**(-3); #this number determines the accuracy of H-primary evaluation higher accuracy=slower time increases ~log(accuracy)

	#reads in sizes of each array
	npos=A.shape[2]; #number of scalp positions
	ncoil=rs.shape[1]; 
	nsource=rv.shape[1];
	ntarget=ncoil*npos;
	
	#generate copies of coil
	rp=np.transpose(np.c_[np.transpose(rs) ,np.ones([ncoil])]); #pads coil positions to 4 by ncoil
	robs=np.array(np.zeros([3,ntarget]));
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		robs[:,st:en]=np.matmul(A[0:3,:,i],rp);
	Etotal=np.zeros([npos,3]);
	kp=ks;
	start = time.time()
	print("Computing H-primary");
	Hprimary=computeHprimary(rv,jvx,robs,prec);
	end = time.time()	
	print("H-primary time: ",end-start);	
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		kp=np.matmul(A[0:3,0:3,i],ks);
		Etotal[i,0]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	start = time.time()
	print("Computing H-primary");
	Hprimary=computeHprimary(rv,jvy,robs,prec);
	end = time.time()	
	print("H-primary time: ",end-start);	
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		kp=np.matmul(A[0:3,0:3,i],ks);
		Etotal[i,1]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	start = time.time()
	print("Computing H-primary");
	Hprimary=computeHprimary(rv,jvz,robs,prec);
	end = time.time()	
	print("H-primary time: ",end-start);	
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		kp=np.matmul(A[0:3,0:3,i],ks);
		Etotal[i,2]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	Etotal=np.sqrt(Etotal[:,0]**2+Etotal[:,1]**2+Etotal[:,2]**2);
	return Etotal;

def ADMmag(rv,jvx,jvy,jvz,rs,ks,A,coildir):
	#this function computes E-field unidirectional approximation of the magnitude via reciprocity with auxiliary dipoles
	#rv is 3 by ntetra and has mesh tetrahedron centroid positions
	#jv is 3 by ntetra and has total conduction current at each tetrahedron
	#rs is 3 by ncoil and has the coil dipole positions centered about the origin
	#ks is 3 by ncoil and has the coil dipole weights
	#A is 4 by 4 by number of coilpositions and each 4 by 4 matrix is a translation of the coil to a point above the scalp
	#coildir is a 3 by number of coil orientations that gives the y-direction orientation


	#paramameter:
	prec=10**(-3); #this number determines the accuracy of H-primary evaluation higher accuracy=slower time increases ~log(accuracy)
	N=[17,17,2];
	#generate auxiliary dipoles
	Nj=coildir.shape[1];
	raux,kaux=resamplecoil(rs,ks,N,Nj,coildir);
	#reads in sizes of each array
	npos=A.shape[2]; #number of scalp positions
	ncoil=raux.shape[1]; 
	nsource=rv.shape[1];
	ntarget=ncoil*npos;
	
	#generate copies of coil
	rp=np.transpose(np.c_[np.transpose(raux) ,np.ones([ncoil])]); #pads coil positions to 4 by ncoil
	robs=np.array(np.zeros([3,ntarget]));
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		robs[:,st:en]=np.matmul(A[0:3,:,i],rp);

	start = time.time()
	Etotal=np.zeros([Nj,npos,3]);
	print("Computing H-primary");
	Hprimary=computeHprimary(rv,jvx,robs,prec);
	end = time.time()
	print("H-primary time: ",end-start);	
	kp=kaux[:,:,0];
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		for j in range(0,Nj):
			kp=np.matmul(A[0:3,0:3,i],kaux[:,:,j]);#for flat coils you only need the 3rd element
			Etotal[j,i,0]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	Hprimary=computeHprimary(rv,jvy,robs,prec);
	end = time.time()
	print("H-primary time: ",end-start);	
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		for j in range(0,Nj):
			kp=np.matmul(A[0:3,0:3,i],kaux[:,:,j]);#for flat coils you only need the 3rd element
			Etotal[j,i,1]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	Hprimary=computeHprimary(rv,jvz,robs,prec);
	end = time.time()
	print("H-primary time: ",end-start);	
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		for j in range(0,Nj):
			kp=np.matmul(A[0:3,0:3,i],kaux[:,:,j]);#for flat coils you only need the 3rd element
			Etotal[j,i,2]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	Etotal=np.sqrt(Etotal[:,:,0]**2+Etotal[:,:,1]**2+Etotal[:,:,2]**2);
	return Etotal;
	

def computeHprimary(rs,js,robs,prec):
	#this function computes H-fields via FMM3D library
	#convention is 3 by number of points
	#prec determines the accuracy of the multipole expansion
	#Note: for magnetic dipoles electromagnetic duality implies that
	#if we pass magnetic dipoles weights as js we get negative E-primary.
	# As such, this function is used to compute E-primary due to magnetic currents also.
	outex=fmm.Output();
	Hprimary=np.zeros([3,robs.shape[1]]);
	muover4pi=-1e-7;
	js1=js[0,:];
	print(rs)
	print(robs)
	out=fmm.lfmm3d(eps=prec,sources=rs,targets=robs,charges=js1,pgt=2);
	print("Run 1", out.gradtarg)
	Hprimary[1,:]=-out.gradtarg[2,:];
	Hprimary[2,:]=out.gradtarg[1,:];
	js1=js[1,:];
	out=fmm.lfmm3d(eps=prec,sources=rs,targets=robs,charges=js1,pgt=2);
	print("Run 2", out.gradtarg)
	Hprimary[0,:]=out.gradtarg[2,:];
	Hprimary[2,:]=Hprimary[2,:]-out.gradtarg[0,:];
	js1=js[2,:];
	out=fmm.lfmm3d(eps=prec,sources=rs,targets=robs,charges=js1,pgt=2);
	print("Run 3", out.gradtarg)
	Hprimary[0,:]=Hprimary[0,:]-out.gradtarg[1,:];
	Hprimary[1,:]=Hprimary[1,:]+out.gradtarg[0,:];
	Hprimary *= muover4pi
	return Hprimary;

def resamplecoil(rs,ks,N,Nj,coildir):
	#rs is 3 by ncoil and has the coil dipole positions centered about the origin
	#ks is 3 by ncoil and has the coil dipole weights
	#N is 3 by 1 and has the number of auxiliary dipoles along each dimension
	#Nj is an integer number of orientations 
	#coildir is 3 by Nj and has the y orientation of the coil
	#create copies of coil with different orientations
	rs2=np.zeros([3,rs.shape[1],Nj]);
	ks2=np.zeros([3,ks.shape[1],Nj]);
	#x=y cross z = y[1] x -y[0] y
	for i in range(0,Nj):
		rs2[2,:,i]=rs[2,:];
		rs2[0,:,i]= rs[0,:]*coildir[1,i]+rs[1,:]*coildir[0,i];
		rs2[1,:,i]=-rs[0,:]*coildir[0,i]+rs[1,:]*coildir[1,i];
		ks2[2,:,i]=ks[2,:];
		ks2[0,:,i]= ks[0,:]*coildir[1,i]+ks[1,:]*coildir[0,i];#unnecessary ks is z oriented
		ks2[1,:,i]=-ks[0,:]*coildir[0,i]+ks[1,:]*coildir[1,i];#unnecessary ks is z oriented
	#find range for interpolation
	Xm=min(rs2[0,:,:].flatten());
	Xp=max(rs2[0,:,:].flatten());
	Ym=min(rs2[1,:,:].flatten());
	Yp=max(rs2[1,:,:].flatten());
	Zm=min(rs2[2,:,:].flatten());
	Zp=max(rs2[2,:,:].flatten());
		#this assumes the original coil has 3 layers and the quadrature is gaussian
		#will need to be changed for other coils
	Zd=.6127016653792583*0.5*(Zp-Zm);
	
	Zm=Zm-Zd;
	Zp=Zp+Zd;

	#get gauss quadrature nodes for interpolation from txt file upto N=50
	interpnodes = np.loadtxt(os.path.join(os.dirname(__file__), 'gauss50.txt'))
	#translate quadrature rules to coil box
	XX=0.5*(interpnodes[int(N[0]*(N[0]-1)/2):int(N[0]*(N[0]+1)/2)]+1.0)*(Xp-Xm)+Xm;
	YY=0.5*(interpnodes[int(N[1]*(N[1]-1)/2):int(N[1]*(N[1]+1)/2)]+1.0)*(Yp-Ym)+Ym;
	ZZ=0.5*(interpnodes[int(N[2]*(N[2]-1)/2):int(N[2]*(N[2]+1)/2)]+1.0)*(Zp-Zm)+Zm;
	#generate auxiliary grid
	Y,Z,X=np.meshgrid(YY,ZZ,XX);
	raux=np.zeros([3,X.size]);
	kaux=np.zeros([3,X.size,Nj]);
	raux[0,:]=X.flatten();
	raux[1,:]=Y.flatten();
	raux[2,:]=Z.flatten();
	del X;
	del Y;
	del Z;
	#generate auxiliary dipole weights
	for kk in range(0,Nj):
		Lx=lagrange(rs2[0,:,kk],XX);
		Ly=lagrange(rs2[1,:,kk],YY);
		Lz=lagrange(rs2[2,:,kk],ZZ);
		for k in range(0,N[2]):
			for j in range(0,N[1]):
				for i in range(0,N[0]):
					L=Lx[i,:]*Ly[j,:]*Lz[k,:];
					kaux[0,int(i+(j+k*N[1])*N[0]),kk]=np.inner(ks2[0,:,kk],L);
					kaux[1,int(i+(j+k*N[1])*N[0]),kk]=np.inner(ks2[1,:,kk],L);
					kaux[2,int(i+(j+k*N[1])*N[0]),kk]=np.inner(ks2[2,:,kk],L);
	return raux,kaux
def lagrange(x,pointx):
	n=pointx.size;
	L=np.ones([n,x.size]);
	for i in range(0,n):
		for j in range(0,n):
			if i != j:
				L[i,:]=np.multiply(L[i,:],(x-pointx[j])/(pointx[i]-pointx[j]));
	return L;

def calculate_Hprim(source_positions, source_currents, query_positions):
    H = np.zeros_like(query_positions)
    for i in range(len(query_positions)):
        r = query_positions[i] - source_positions
        dist = np.linalg.norm(r, axis=1)
        H[i] = 1e-7 * np.sum(np.cross(source_currents, r)/dist[:, None]**3, axis=0)

    return H

def test_computeHprimary():
    np.random.seed(1)
    current_elm_positions = np.random.rand(100, 3)
    currents = np.random.rand(len(current_elm_positions), 3)
    observation_pos = np.random.rand(10, 3)

    fmm_Hprimary = computeHprimary(
        current_elm_positions.T,
        currents.T,
        observation_pos.T,
        1e-5
    ).T

    Hprim = calculate_Hprim(
        current_elm_positions,
        currents,
        observation_pos
    )

    print(fmm_Hprimary)
    print(Hprim)
    assert np.allclose(fmm_Hprimary, Hprim)


if __name__ == '__main__':
    test_computeHprimary()
