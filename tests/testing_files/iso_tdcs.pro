// simnibs automatically generated file for getdp simulation
// file created with getdp.py on 20150309_173537


WM = 1;
WM_surface = 1001;
GM = 2;
GM_surface = 1002;
CSF = 3;
CSF_surface = 1003;
Bone = 4;
Bone_surface = 1004;
Scalp = 5;
Scalp_surface = 1005;
Spongy_bone = 6;
Spongy_bone_surface = 1006;
Compact_bone = 7;
Compact_bone_surface = 1007;
Eye_balls = 8;
Eye_balls_surface = 1008;
Eye_region = 9;
Eye_region_surface = 1009;
el1 = 100;
el1_surface = 1100;
el2 = 101;
el2_surface = 1101;
gel_sponge1 = 500;
gel_sponge1_surface = 1500;
gel_sponge2 = 501;
gel_sponge2_surface = 1501;


Group {
	WM = Region[WM];
	WM_surface = Region[WM_surface];
	GM = Region[GM];
	GM_surface = Region[GM_surface];
	CSF = Region[CSF];
	CSF_surface = Region[CSF_surface];
	Bone = Region[Bone];
	Bone_surface = Region[Bone_surface];
	Scalp = Region[Scalp];
	Scalp_surface = Region[Scalp_surface];
	Spongy_bone = Region[Spongy_bone];
	Spongy_bone_surface = Region[Spongy_bone_surface];
	Compact_bone = Region[Compact_bone];
	Compact_bone_surface = Region[Compact_bone_surface];
	Eye_balls = Region[Eye_balls];
	Eye_balls_surface = Region[Eye_balls_surface];
	Eye_region = Region[Eye_region];
	Eye_region_surface = Region[Eye_region_surface];
	el1 = Region[el1];
	el1_surface = Region[el1_surface];
	el2 = Region[el2];
	el2_surface = Region[el2_surface];
	gel_sponge1 = Region[gel_sponge1];
	gel_sponge1_surface = Region[gel_sponge1_surface];
	gel_sponge2 = Region[gel_sponge2];
	gel_sponge2_surface = Region[gel_sponge2_surface];
	Omega = Region[{WM,GM,CSF,Bone,Scalp,Spongy_bone,Compact_bone,Eye_balls,Eye_region,el1,el2,gel_sponge1,gel_sponge2}];
	Gama = Region[{WM_surface, GM_surface, CSF_surface, Bone_surface, Scalp_surface, Spongy_bone_surface, Compact_bone_surface, Eye_balls_surface, Eye_region_surface, el1_surface, el2_surface, gel_sponge1_surface, gel_sponge2_surface}];
}


Function {
	sigma[WM] = 0.126;
	sigma[GM] = 0.275;
	sigma[CSF] = 1.654;
	sigma[Bone] = 0.01;
	sigma[Scalp] = 0.465;
	sigma[Spongy_bone] = 0.025;
	sigma[Compact_bone] = 0.008;
	sigma[Eye_balls] = 0.5;
	sigma[Eye_region] = 0.25;
	sigma[el1] = 0.1;
	sigma[el2] = 0.1;
	sigma[gel_sponge1] = 1.0;
	sigma[gel_sponge2] = 1.0;
}


Constraint {{
	Name ElectricScalarPotential;
	Type Assign;
	Case {
		{Region Region[{el1_surface}]; Value 1.;}
		{Region Region[{el2_surface}]; Value -1.;}
	}
}}


Jacobian {
	{ Name Volume;  Case {{ Region Omega; Jacobian Vol;  }} }
}


Integration {{
	Name GradGrad;
	Case {{
		Type Gauss;
		Case {
			{GeoElement Tetrahedron; NumberOfPoints 1;}
		}
	}}
}}


FunctionSpace {{
	Name Hgrad_vf_Ele;
	Type Form0;
	BasisFunction {{
		Name sn;
		NameOfCoef vn;
		Function BF_Node;
		Support Region[{Omega}];
		Entity NodesOf[All];
	}}
	Constraint {{
		NameOfCoef vn;
		EntityType NodesOf;
		NameOfConstraint ElectricScalarPotential;
	}}
}}


Formulation {{
	Name Electrostatics_Formulation;
	Type FemEquation;
	Quantity {{
		Name v;
		Type Local;
		NameOfSpace Hgrad_vf_Ele;
	}}
	Equation {
		Galerkin {
			[sigma[] * Dof{Grad v}, {Grad v}];
			In Omega;
			Jacobian Volume;
			Integration GradGrad;
		}
	}
}}


Resolution {{
	Name QS;
	System {{
		Name Electrostatic_System;
		NameOfFormulation Electrostatics_Formulation;
	}}
	Operation {
		Generate Electrostatic_System;
		Solve Electrostatic_System;
		SaveSolution Electrostatic_System;
	}
}}


